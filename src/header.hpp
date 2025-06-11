#ifndef HEADER_HPP
#define HEADER_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <Rcpp.h>

enum Status { SUSCEPTIBLE, INFECTIOUS, RECOVERED };

class Agent {
public:
  int id;
  int age_group;
  Status status;
  int day_of_infection;
  std::vector<double> viral_load_weights;

  // Constructor for Agent class -----------------------------------------------

  Agent(int _id, int _age_group)
    : id(_id), age_group(_age_group), status(SUSCEPTIBLE), day_of_infection(-1) {}

  // Infect the agent and define a viral load profile --------------------------

  void infect(int day, const std::vector<double>& viral_load)
  {
    status = INFECTIOUS;
    day_of_infection = day;
    viral_load_weights = viral_load;
  }

  // Recover from the disease --------------------------------------------------

  void recover()
  {
    status = RECOVERED;
    day_of_infection = -1;
    viral_load_weights.clear();
  }

  // Define infectiousness at time t since infection, based on the viral load --

  double infectiousness(int t_since_infection) const
  {
    if ( t_since_infection < 0 || t_since_infection >= static_cast<int>(viral_load_weights.size()))
      return 0.0;
    return viral_load_weights[t_since_infection];
  }

  // Check if the agent has recovered based on the current day -----------------

  bool has_recovered(int current_day, double recovery_threshold ) const
  {
    int t = current_day - day_of_infection;

    if (t < 0)
    {
      Rcpp::Rcout << "Warning: has_recovered() called before infection.\n";
      return false;
    }

    // Two cases for recovering:
    // 1) Duration of infection has passed the length of the viral load profile
    // 2) Viral load has dropped below the recovery threshold after the peak

    // Case 1.

    if (t >= static_cast<int>(viral_load_weights.size()))
      return true;   // Infection has run its course --> recovered

    // Case 2.

    // Find the peak of the viral load

    auto it = std::max_element(viral_load_weights.begin(), viral_load_weights.end());
    int peak_day = std::distance(viral_load_weights.begin(), it);

    // If we are past the peak and VL has dropped below threshold

    if (t > peak_day && viral_load_weights[t] < recovery_threshold)
      return true;   // Infection has run its course --> recovered

    return false;
  }
};

// Forward declaration

std::vector<double> simulate_viral_load(int days, double dt,
                                        double V0, double dV0,
                                        double k, double gamma, double sigma,
                                        std::mt19937& rng);
std::vector<double> simulate_viral_load_v2(int days, double dt,
                                            double V0, double dV0,
                                            double k, double gamma, double sigma,
                                            std::mt19937& rng);
Rcpp::NumericVector test_viral_load_trajectory(int age_group,
                                               double dt,
                                               int VL_days,
                                               double I0,
                                               double baseline,
                                               Rcpp::NumericVector k,
                                               Rcpp::NumericVector gamma,
                                               Rcpp::NumericVector sigma);

// Validation functions
void check_initial_infections(int n_agents, int n_initial_infected);
void check_contact_matrix(const Rcpp::NumericMatrix& contact_matrix, const Rcpp::NumericVector& age_probs);

#endif
