library(MultiScaleABM)
library(ggplot2)

# Define the parameters for the ABM
params <- list(
  n_agents = 1e4,
  n_days = 180,
  contact_matrix = matrix(c(0.1, 0.05, 0.01,
                            0.05, 0.1, 0.02,
                            0.01, 0.02, 0.1),
                          nrow = 3, byrow = TRUE),
  age_groups = c(0.2, 0.5, 0.3),
  k = c(0.05, 0.05, 0.05),
  gamma = c(0.3, 0.3, 0.3),
  sigma = c(0.25, 0.25, 0.25),
  beta = 10,
  dt = 0.0005,
  VL_days = 20,
  V0 = log(0.01),
  dV0 = 10,
  recovery_rate = 0.1,
  n_initial_infected = 10
)

# Check viral load trajectory. If spaghetti=T, plot N stochastic trajectories
plot_viral_load_trajectory(params = params, spaghetti = T, N=100)

# Run the agent-based model
out <- run_abm(params)

# Plot the results

d <- data.frame(  day        = 1:params$n_days,
                  infected   = out$daily_infected,
                  recovered  = out$daily_recovered,
                  prevalence = out$prevalence )

p <- ggplot(d, aes(x = day)) +
      geom_line(aes(y = infected, color = "Infected")) +
      geom_line(aes(y = recovered, color = "Recovered")) +
      labs( title = "Agent-Based Model Simulation Results",
            x = "Day",
            y = "Incidence") +
      scale_color_manual(values = c("Infected" = "red",
                                    "Recovered" = "green",
                                    "Prevalence" = "blue")) +
      theme_minimal()
print(p)

# Save the plot
# Create folder if it doesn't exist
#if (!dir.exists("plots")) {
#  dir.create("plots")
#}
#ggsave("plots/abm_simulation_results.png", plot = p, width = 10, height = 6)
