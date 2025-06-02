#include <Rcpp.h>
#include "header.hpp"

using namespace Rcpp;

// std::vector<double> simulate_viral_load(int days, double dt,
//                                         double V0, double dV0,
//                                         double k, double gamma, double sigma)
// {
//   std::vector<double> V(days), dV(days), Y(days), w(days);
//   V[0] = V0;
//   dV[0] = dV0;
//   Y[0] = std::exp(V0);
//
//   for (int t = 1; t < days; ++t)
//   {
//     double dW = R::rnorm(0, std::sqrt(dt));
//     double ddV = -k * V[t - 1] - gamma * dV[t - 1];
//     dV[t] = dV[t - 1] + ddV * dt;
//     V[t] = V[t - 1] + dV[t] * dt + sigma * dW; // stochastic part
//     Y[t] = Y[t - 1] + Y[t - 1] * dV[t] * dt + 0.5 * Y[t - 1] * sigma * sigma * dt; // drift part
//   }
//
//   double maxY = *std::max_element(Y.begin(), Y.end());
//   for (int t = 0; t < days; ++t) w[t] = Y[t] / maxY;
//
//   return w;
// }


// [[Rcpp::export]]
std::vector<double> simulate_viral_load(int days, double dt,
                                        double V0, double dV0,
                                        double k, double gamma, double sigma)
{
  int n_steps = static_cast<int>(days / dt);
  std::vector<double> V(n_steps), dV(n_steps), Y(n_steps);
  V[0] = V0;
  dV[0] = dV0;
  Y[0] = std::exp(V0);

  // Simulate viral dynamics
  for (int t = 1; t < n_steps; ++t)
  {
    double dW = R::rnorm(0.0, std::sqrt(dt));
    double ddV = -k * V[t - 1] - gamma * dV[t - 1];
    dV[t] = dV[t - 1] + ddV * dt ;
    V[t]  = V[t - 1] + dV[t] * dt+ sigma * dW ;
    //Y[t]  = Y[t - 1] + Y[t - 1] * dV[t] * dt + 0.5 * Y[t - 1] * sigma * sigma * dt;
    Y[t] = std::exp(V[t]);
  }

  // Compute daily averages
  std::vector<double> daily_Y(days, 0.0);
  int steps_per_day = static_cast<int>(1.0 / dt);

  for (int d = 0; d < days; ++d)
  {
    double sum = 0.0;
    int count = 0;

    // Average over time steps within day d
    for (int t = d * steps_per_day; t < (d + 1) * steps_per_day && t < n_steps; ++t)
    {
      sum += Y[t];
      count++;
    }

    if (count > 0)
    {
      daily_Y[d] = sum / count;
    }
    else
    {
      daily_Y[d] = 0.0;
    }
  }

  // Normalize the curve
  double maxY = *std::max_element(daily_Y.begin(), daily_Y.end());
  std::vector<double> w(days, 0.0);

  for (int d = 0; d < days; ++d)
  {
    if (maxY > 0.0)
    {
      w[d] = daily_Y[d] / maxY;
    }
    else
    {
      w[d] = 0.0;
    }
  }

  return w;
}

// [[Rcpp::export]]
Rcpp::NumericVector test_viral_load_trajectory(int age_group,
                                               double dt,
                                               int VL_days,
                                               double I0,
                                               double baseline,
                                               Rcpp::NumericVector k,
                                               Rcpp::NumericVector gamma,
                                               Rcpp::NumericVector sigma)
{

  std::vector<double> vl = simulate_viral_load(VL_days, dt, I0, baseline,
                                               k[age_group], gamma[age_group], sigma[age_group]);

  return Rcpp::wrap(vl);
}
