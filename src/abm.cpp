#include <Rcpp.h>
#include "header.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List run_abm(List params) {

  // Extract parameters from the input list

  int n_agents                 = params["n_agents"];
  int n_days                   = params["n_days"];
  NumericMatrix contact_matrix = params["contact_matrix"];
  NumericVector age_probs      = params["age_groups"];
  NumericVector k              = params["k"];
  NumericVector gamma          = params["gamma"];
  NumericVector sigma          = params["sigma"];
  double beta                  = params["beta"];
  double dt                    = params["dt"];
  int VL_days                  = params["VL_days"];
  double V0                    = params["V0"];
  double dV0                   = params["dV0"];
  double recovery_rate         = params["recovery_rate"];
  int n_initial_infected       = params["n_initial_infected"];
  double prob_recovery         = (1-exp(-recovery_rate));
  // Some checks

  check_initial_infections(n_agents, n_initial_infected);
  check_contact_matrix(contact_matrix, age_probs);

  // Define age groups

  IntegerVector age_labels         = seq(0, age_probs.size() - 1);
  IntegerVector sampled_age_groups = Rcpp::sample(age_labels, n_agents, true, age_probs);

  // Define agents' age groups and initialize them

  std::vector<Agent> agents;
  for (int i = 0; i < n_agents; ++i)
  {
    agents.emplace_back(i, sampled_age_groups[i]);
  }

  // Set up the random number generator

  RNGScope scope;

  // Epidemic Initialization: Sample n_initial_infected individuals and infect them

  Rcpp::Rcout << "Number of infected agents at day 0: " << n_initial_infected << std::endl;

  IntegerVector initial_cases = Rcpp::sample(n_agents, n_initial_infected, false);

  for (int i = 0; i < n_initial_infected; ++i)
  {
    int idx = initial_cases[i];
    int ag = agents[idx].age_group;

    std::vector<double> vl = simulate_viral_load(VL_days, dt, V0, dV0, k[ag], gamma[ag], sigma[ag]);
    agents[idx].infect(0, vl);
  }

  // Initialize vectors to store daily infected, recovered, and prevalence

  IntegerVector daily_infected(n_days);
  IntegerVector daily_recovered(n_days);
  IntegerVector prevalence(n_days);

  daily_infected[0] = n_initial_infected;
  prevalence[0]     = n_initial_infected;

  // Loop over time

  for (int day = 1; day < n_days; ++day)
  {

    int current_prevalence = 0;

    // Recovery loop

    for (int i = 0; i < n_agents; ++i)
    {
      if (agents[i].status == INFECTIOUS)
      {
        if ( R::runif(0, 1) < prob_recovery )
        {
          agents[i].recover();
          daily_recovered[day]++;
        }
        else
        {
          current_prevalence++;
        }
      }
    }

    // Compute nr. of people in each age group to speed up the infectious process

    int n_age_groups = age_probs.size();
    std::vector<double> group_infectiousness(n_age_groups, 0.0);
    std::vector<int> n_agents_by_age(n_age_groups, 0);

    // Aggregate infectiousness and count group sizes

    for (int i = 0; i < n_agents; ++i)
    {
      int ag = agents[i].age_group;
      n_agents_by_age[ag]++;

      if (agents[i].status == INFECTIOUS)
      {
        int t = day - agents[i].day_of_infection;
        if (t >= 0 && t < agents[i].viral_load_weights.size())
        {
          group_infectiousness[ag] += agents[i].infectiousness(t);
        }
      }
    }

    // Compute FOI (lambda) for each age group

    std::vector<double> lambda_by_age(n_age_groups, 0.0);

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
    }

    // Infect susceptible agents based on their age group's FOI

    for (int i = 0; i < n_agents; ++i)
    {
      if (agents[i].status != SUSCEPTIBLE) continue;

      int ag = agents[i].age_group;
      double lambda = lambda_by_age[ag];
      double pi = 1.0 - std::exp(-lambda);

      if ( R::runif(0, 1) < pi )
      {
        std::vector<double> vl = simulate_viral_load(
          VL_days, dt, V0, dV0, k[ag], gamma[ag], sigma[ag]
        );
        agents[i].infect(day, vl);
        daily_infected[day]++;
        current_prevalence++;
      }
    }

    prevalence[day] = current_prevalence;

  } // End loop over time

  // Output the results

  return List::create(
    Named("daily_infected") = daily_infected,
    Named("daily_recovered") = daily_recovered,
    Named("prevalence") = prevalence
  );
}
