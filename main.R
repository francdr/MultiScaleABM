# ==============================================================================
# Agent-Based Epidemic Simulation using MultiScaleABM
# ==============================================================================

library(MultiScaleABM)
library(ggplot2)

# ------------------------------------------------------------------------------
# Load contact matrix and age distribution from CSV files
# ------------------------------------------------------------------------------

contact_matrix <- as.matrix(read.csv("data/contact_matrix.csv", header = FALSE))
age_distribution <- as.numeric(read.csv("data/age_distribution.csv", header = FALSE)[, 2])

stopifnot(abs( sum(age_distribution/sum(age_distribution)) - 1) < 1e-6)  # sanity check: must sum to 1

# ------------------------------------------------------------------------------
# Set Parameters
# ------------------------------------------------------------------------------

params <- list(
  # ABM parameters
  n_agents = 1e5,
  n_days = 180,
  contact_matrix = contact_matrix,
  age_groups = age_distribution,
  # Viral load parameters
  k = rep(0.06, length(age_distribution)),
  gamma = rep(0.3, length(age_distribution)),
  sigma = rep(0.36, length(age_distribution)),
  beta = 0.08,
  dt = 0.05,
  VL_days = 20,
  V0 = log(0.06),
  dV0 = 5,
  # Simulation parameters
  n_initial_infected = 10,
  recovery_threshold = 0.01,
  parallel = TRUE,
  n_threads = 12,
  verbose = T
)

# ------------------------------------------------------------------------------
# Plot Sample Viral Load Trajectories
# ------------------------------------------------------------------------------

plot_viral_load_trajectory(params = params, spaghetti = TRUE, N = 100)

# ------------------------------------------------------------------------------
# Run ABM and Plot SIR Dynamics
# ------------------------------------------------------------------------------

out <- run_abm(params)

sir_data <- data.frame(
  day = 1:params$n_days,
  infected = out$daily_infected,
  recovered = out$daily_recovered,
  S = out$prevalence_SUSCEPTIBLE,
  I = out$prevalence_INFECTIOUS,
  R = out$prevalence_RECOVERED
)
sir_data$Total <- rowSums(sir_data[, c("S", "I", "R")])

p_sir <- ggplot(sir_data, aes(x = day)) +
  geom_line(aes(y = S, color = "Susceptible"), size = 1.2) +
  geom_line(aes(y = I, color = "Infectious"), size = 1.2) +
  geom_line(aes(y = R, color = "Recovered"), size = 1.2) +
  labs(title = "",
       x = "Day",
       y = "Number of Individuals",
       color = "Compartment") +
  scale_color_manual(values = c(
    Susceptible = "#1b9e77",  # teal green
    Infectious  = "#d95f02",  # rich orange
    Recovered   = "#7570b3"   # muted purple
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_sir)
# ggsave("plots/abm_simulation_results.png", p_sir, width = 10, height = 6)

# ------------------------------------------------------------------------------
# Run Multiple Stochastic Realizations and Make a Spaghetti Plot
# ------------------------------------------------------------------------------

plot_spaghetti <- function( params, N = 10 )
{
  results <- data.frame()
  seed_sequence <- sample(1:100000, N)

  for (seed in seed_sequence) {
    params$seed <- seed
    out <- run_abm(params)
    tmp <- data.frame(
      day = 1:params$n_days,
      S = out$prevalence_SUSCEPTIBLE,
      I = out$prevalence_INFECTIOUS,
      R = out$prevalence_RECOVERED,
      seed = seed
    )
    results <- rbind(results, tmp)
  }

  ggplot(results, aes(x = day, group = seed)) +
    geom_line(aes(y = S, color = "Susceptible"), alpha = 0.4, size = 0.8) +
    geom_line(aes(y = I, color = "Infectious"),  alpha = 0.4, size = 0.8) +
    geom_line(aes(y = R, color = "Recovered"),   alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c(
      Susceptible = "#1b9e77",  # teal green
      Infectious  = "#d95f02",  # rich orange
      Recovered   = "#7570b3"   # muted purple
    )) +
    labs(
      title = paste("Spaghetti Plot of", N, "Simulations"),
      x = "Day",
      y = "Number of Individuals",
      color = "Compartment"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

p_spaghetti <- plot_spaghetti(params, N = 50)
print(p_spaghetti + xlim(0, 80) )

