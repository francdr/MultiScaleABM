#include <Rcpp.h>
#include <vector>
#include <thread>
#include <future>
#include <random>
#include "header.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List run_abm(List params) {

  // ===========================================================================
  // Parse Parameters from Input List
  // ===========================================================================

  int n_agents            = params["n_agents"];
  int n_days              = params["n_days"];
  NumericMatrix contact_matrix = params["contact_matrix"];
  NumericVector age_probs = params["age_groups"];
  NumericVector k         = params["k"];
  NumericVector gamma     = params["gamma"];
  NumericVector sigma     = params["sigma"];
  double beta             = params["beta"];
  double dt               = params["dt"];
  int VL_days             = params["VL_days"];
  double V0               = params["V0"];
  double dV0              = params["dV0"];
  int n_initial_infected  = params["n_initial_infected"];
  double recovery_threshold = params["recovery_threshold"];

  // ---------------------------------------------------------------------------
  // Optional Flags with Defaults
  // ---------------------------------------------------------------------------

  bool use_parallel = params.containsElementNamed("parallel") ?
                      as<bool>(params["parallel"]) : false;

  int n_threads = (use_parallel && params.containsElementNamed("n_threads")) ?
                  as<int>(params["n_threads"]) :
                  std::thread::hardware_concurrency();

  int seed = params.containsElementNamed("seed") ?
             as<int>(params["seed"]) :
             std::random_device{}();

  bool VERBOSE = params.containsElementNamed("verbose") ?
                 as<bool>(params["verbose"]) : false;

  // ---------------------------------------------------------------------------
  // Verbose Mode Output
  // ---------------------------------------------------------------------------

  if (VERBOSE)
  {
    Rcpp::Rcout << "Using seed: " << seed << std::endl;
    Rcpp::Rcout << "Parallel: " << (use_parallel ? "Yes" : "No")
                << " | Threads: " << n_threads << std::endl;
  }

  // ===========================================================================
  // Set the random seed for reproducibility
  // ===========================================================================

  std::mt19937 rng(seed);
  std::uniform_real_distribution<> runif(0.0, 1.0);

  // ===========================================================================
  // Check
  // ===========================================================================

  check_initial_infections(n_agents, n_initial_infected);
  check_contact_matrix(contact_matrix, age_probs);

  // ===========================================================================
  // Define age profile
  // ===========================================================================

  std::discrete_distribution<int> age_dist(age_probs.begin(), age_probs.end());

  std::vector<int> sampled_age_groups(n_agents);
  for (int i = 0; i < n_agents; ++i)
  {
    sampled_age_groups[i] = age_dist(rng);
  }

  // ---------------------------------------------------------------------------
  // Initialize agents with sampled age groups
  // ---------------------------------------------------------------------------

  std::vector<Agent> agents;
  for (int i = 0; i < n_agents; ++i)
  {
    agents.emplace_back(i, sampled_age_groups[i]);
  }

  // ---------------------------------------------------------------------------
  // Count agents by age group
  // ---------------------------------------------------------------------------

  std::vector<int> n_agents_by_age(age_probs.size(), 0);
  for (int i = 0; i < n_agents; ++i)
  {
    int ag = agents[i].age_group;
    n_agents_by_age[ag]++;
  }

 // Print proportions of agents by age group based on n_agents_by_age
  if (VERBOSE)
  {
    Rcpp::Rcout << "Proportions of agents by age group:" << std::endl;

    for (size_t i = 0; i < n_agents_by_age.size(); ++i)
    {
      double proportion = static_cast<double>(n_agents_by_age[i]) / n_agents;
      Rcpp::Rcout << "Age group " << i << ": " << proportion << std::endl;
    }
  }

  // ===========================================================================
  // Initialize Agents with Initial Infections
  // ===========================================================================

  if (VERBOSE)
  {
    Rcpp::Rcout << "Number of infected agents at day 0: "
                << n_initial_infected << std::endl;
  }

  // Generate a list of all agent indices: 0, 1, ..., n_agents - 1
  std::vector<int> all_indices(n_agents);
  std::iota(all_indices.begin(), all_indices.end(), 0);

  // Randomly sample agent indices to infect initially
  std::vector<int> initial_cases;
  std::sample( all_indices.begin(),
               all_indices.end(),
               std::back_inserter(initial_cases),
               n_initial_infected,
               rng);

  // Infect selected agents
  for (int idx : initial_cases)
  {
    int ag = agents[idx].age_group;

    agents[idx].infect( 0, // infection starts at day 0
                        simulate_viral_load(VL_days, dt, V0, dV0, k[ag], gamma[ag], sigma[ag], rng));
  }

  // ===========================================================================
  // Initialize Output Vectors
  // ===========================================================================

  IntegerVector daily_infected(n_days);
  IntegerVector daily_recovered(n_days);
  IntegerVector prevalence_INFECTIOUS(n_days);
  IntegerVector prevalence_SUSCEPTIBLE(n_days);
  IntegerVector prevalence_RECOVERED(n_days);

  // Initial conditions at time t = 0
  daily_infected[0]         = n_initial_infected;
  prevalence_INFECTIOUS[0]  = n_initial_infected;
  prevalence_SUSCEPTIBLE[0] = n_agents - n_initial_infected;
  prevalence_RECOVERED[0]   = 0;

  // ============================================================================
  // Begin Time Loop
  // ============================================================================

  for (int day = 1; day < n_days; ++day)
  {

    int current_prevalence_INFECTIOUS  = 0;
    int current_prevalence_SUSCEPTIBLE = 0;
    int current_prevalence_RECOVERED   = 0;

    // =========================================================================
    // Sequential execution
    // =========================================================================

    if( use_parallel == FALSE)
    {

      // -----------------------------------------------------------------------
      // Compute Group-Level Infectiousness
      // -----------------------------------------------------------------------

      int n_age_groups = age_probs.size();
      std::vector<double> group_infectiousness(n_age_groups, 0.0);

      for (int i = 0; i < n_agents; ++i)
      {
        int ag = agents[i].age_group;

        if (agents[i].status == INFECTIOUS)
        {
          int t = day - agents[i].day_of_infection;

          if (t >= 0 && t < static_cast<int>(agents[i].viral_load_weights.size()))
          {
            group_infectiousness[ag] += agents[i].infectiousness(t);
          }
        }
      }

      // -----------------------------------------------------------------------
      // Calculate the FOI
      // -----------------------------------------------------------------------

      std::vector<double> lambda_by_age(n_age_groups, 0.0);
      std::vector<double> pi_by_age(n_age_groups, 0.0);

      for (int i = 0; i < n_age_groups; ++i)
      {
        for (int j = 0; j < n_age_groups; ++j)
        {
          if (n_agents_by_age[j] > 0)
          {
            lambda_by_age[i] += beta * contact_matrix(i, j) *
            group_infectiousness[j] / n_agents_by_age[j];
          }
        }

        // Calculate the probability of infection for each age group

        pi_by_age[i] = 1.0 - std::exp(-lambda_by_age[i]);
      }

      // -----------------------------------------------------------------------
      // Update the agents' epidemiological states
      // -----------------------------------------------------------------------

      for (int i = 0; i < n_agents; ++i)
      {
        switch (agents[i].status)
        {
          case SUSCEPTIBLE:
          {
            int ag = agents[i].age_group;

            if ( runif(rng) < pi_by_age[ag])
            {
              // Simulate viral load trajectory for new infection
              auto vload = simulate_viral_load(
                            VL_days, dt, V0, dV0,
                            k[ag], gamma[ag], sigma[ag], rng);

              agents[i].infect(day, vload);
              daily_infected[day] ++;
              current_prevalence_INFECTIOUS++;
            }
            else
            {
              current_prevalence_SUSCEPTIBLE++;
            }
            break;
          }

          case INFECTIOUS:
          {
            // Check for recovery based on viral load and threshold
            if ( agents[i].has_recovered(day, recovery_threshold) )
            {
              agents[i].recover();
              daily_recovered[day]++;
              current_prevalence_RECOVERED++;
            }
            else
            {
              current_prevalence_INFECTIOUS++;
            }
            break;
          }

          case RECOVERED:
          {
            current_prevalence_RECOVERED++;
            break;
          }
        } // end switch
      } // end agent loop
    }

    // =========================================================================
    // Parallel execution
    // =========================================================================

    else
    {
      int n_age_groups = age_probs.size();

      // Define population chunk size

      int chunk_size = n_agents / n_threads;

      // -----------------------------------------------------------------------
      // Compute Group-Level Infectiousness
      // -----------------------------------------------------------------------

      std::vector<std::future<void>> futures;
      std::vector<std::vector<double>> local_group_inf(n_threads, std::vector<double>(n_age_groups, 0.0));

      for (unsigned int t = 0; t < n_threads; ++t)
      {
        int start = t * chunk_size;
        int end;

        // Ensure the last thread takes any remaining agents
        if (t == n_threads - 1)
        {
          end = n_agents;  // Last thread takes the remainder
        }
        else
        {
          end = (t + 1) * chunk_size;
        }

        futures.push_back(std::async(std::launch::async, [&, t, start, end]()
        {
          for (int i = start; i < end; ++i)
          {
            int ag = agents[i].age_group;

            if (agents[i].status == INFECTIOUS)
            {
              int t_inf = day - agents[i].day_of_infection;

              if (t_inf >= 0 && t_inf < static_cast<int>(agents[i].viral_load_weights.size()))
              {
                local_group_inf[t][ag] += agents[i].infectiousness(t_inf);
              }
            }
          }
        }));
      }

      // Wait for all threads to complete before continuing

      for (auto& f : futures) f.get();

      // Combine results from all threads

      std::vector<double> group_infectiousness(n_age_groups, 0.0);

      for (unsigned int t = 0; t < n_threads; ++t)
      {
        for (int ag = 0; ag < n_age_groups; ++ag)
        {
          group_infectiousness[ag] += local_group_inf[t][ag];
        }
      }

      // -----------------------------------------------------------------------
      // Calculate the FOI
      // -----------------------------------------------------------------------

      std::vector<double> lambda_by_age(n_age_groups, 0.0);
      std::vector<double> pi_by_age(n_age_groups, 0.0);

      for (int i = 0; i < n_age_groups; ++i)
      {
        for (int j = 0; j < n_age_groups; ++j)
        {
          if (n_agents_by_age[j] > 0)
          {
            lambda_by_age[i] += beta * contact_matrix(i, j) *
              group_infectiousness[j] / n_agents_by_age[j];
          }
        }

        // Calculate the probability of infection for each age group

        pi_by_age[i] = 1.0 - std::exp(-lambda_by_age[i]);
      }

      // -----------------------------------------------------------------------
      // Update the agents' epidemiological states
      // -----------------------------------------------------------------------

      std::vector<std::future<void>> state_futures;

      // Local counters for this thread
      std::vector<int> local_infected(n_threads, 0);
      std::vector<int> local_prevalence_INFECTIOUS(n_threads, 0);
      std::vector<int> local_prevalence_SUSCEPTIBLE(n_threads, 0);
      std::vector<int> local_prevalence_RECOVERED(n_threads, 0);
      std::vector<int> local_recovered(n_threads, 0);

      for (unsigned int t = 0; t < n_threads; ++t)
      {
        int start = t * chunk_size;
        int end;

        // Ensure the last thread takes any remaining agents

        if (t == n_threads - 1)
        {
          end = n_agents;  // Last thread takes the remainder
        } else {
          end = (t + 1) * chunk_size;
        }

        state_futures.push_back(std::async(std::launch::async, [&, t, start, end]()
        {
          // Use Cantor Pairing function to generate a unique seed for each thread
          // int unique_seed = (t + day) * (t + day + 1) / 2 + day;
          int unique_seed = seed + t*100000 + day;
          std::mt19937 rng_local(unique_seed);
          // rng_local.discard(start);      // skip to start of this thread's chunk
          std::uniform_real_distribution<> runif(0.0, 1.0);

          for (int i = start; i < end; ++i)
          {
            switch (agents[i].status)
            {
              case SUSCEPTIBLE:
              {
                int ag = agents[i].age_group;

                if (runif(rng_local) < pi_by_age[ag])
                {
                  // Simulate viral load trajectory for new infection
                  auto vload = simulate_viral_load(
                                VL_days, dt, V0, dV0,
                                k[ag], gamma[ag], sigma[ag], rng_local);
                  agents[i].infect(day, vload);
                  local_infected[t]++;
                  local_prevalence_INFECTIOUS[t]++;
                }
                else
                {
                  local_prevalence_SUSCEPTIBLE[t]++;
                }
                break;
              }

              case INFECTIOUS:
              {
                // Check for recovery based on viral load and threshold
                if ( agents[i].has_recovered(day, recovery_threshold) )
                {
                  agents[i].recover();
                  local_recovered[t]++;
                  local_prevalence_RECOVERED[t]++;
                }
                else
                {
                  local_prevalence_INFECTIOUS[t]++;
                }
                break;
              }

              case RECOVERED:
                local_prevalence_RECOVERED[t]++;
                break;
            }
          }
        }));
      }

      // Wait for all threads to complete before continuing
      for (auto& f : state_futures) f.get();

      // Combine results from all threads

      for (unsigned int t = 0; t < n_threads; ++t)
      {
          daily_infected[day]  += local_infected[t];
          daily_recovered[day] += local_recovered[t];
          current_prevalence_INFECTIOUS  += local_prevalence_INFECTIOUS[t];
          current_prevalence_SUSCEPTIBLE += local_prevalence_SUSCEPTIBLE[t];
          current_prevalence_RECOVERED   += local_prevalence_RECOVERED[t];
      }
    }

    prevalence_INFECTIOUS[day] = current_prevalence_INFECTIOUS;
    prevalence_SUSCEPTIBLE[day] = current_prevalence_SUSCEPTIBLE;
    prevalence_RECOVERED[day]   = current_prevalence_RECOVERED;
  }

  return List::create(
    Named("daily_infected") = daily_infected,
    Named("daily_recovered") = daily_recovered,
    Named("prevalence_INFECTIOUS") = prevalence_INFECTIOUS,
    Named("prevalence_SUSCEPTIBLE") = prevalence_SUSCEPTIBLE,
    Named("prevalence_RECOVERED") = prevalence_RECOVERED
  );
}
