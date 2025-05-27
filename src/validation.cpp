#include <Rcpp.h>
using namespace Rcpp;

// Check that n_initial_infected is <= n_agents
void check_initial_infections(int n_agents, int n_initial_infected) {
  if (n_initial_infected > n_agents) {
    Rcpp::stop("n_initial_infected cannot be larger than n_agents");
  }
}

// Check that contact matrix matches age group structure
void check_contact_matrix(const NumericMatrix& contact_matrix, const NumericVector& age_probs) {
  if (contact_matrix.nrow() != age_probs.size() || contact_matrix.ncol() != age_probs.size()) {
    Rcpp::stop("Contact matrix must be square and match number of age groups");
  }
}

