# ==============================================================================
# Agent-Based Epidemic Simulation using MultiScaleABM
# ==============================================================================

library(MultiScaleABM)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# Load contact matrix and age distribution from CSV files
# ------------------------------------------------------------------------------

contact_matrix   <- as.matrix(read.csv("data/contact_matrix.csv",    header = FALSE))
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
  sigma = rep(0.2, length(age_distribution)),
  beta = 1.8e-5,
  dt = 0.05,
  VL_days = 20,
  V0 = log(0.06),
  dV0 = 5,
  # Simulation parameters
  n_initial_infected = 10,
  recovery_threshold = 20,#0.01,
  parallel = TRUE,
  n_threads = 12,
  verbose = T
)

# ------------------------------------------------------------------------------
# Plot Sample Viral Load Trajectories
# ------------------------------------------------------------------------------

p_vl <- plot_viral_load_trajectory(params = params, spaghetti = TRUE, N = 500)
p_vl <- p_vl + labs(title="", x = "Day", y="Viral Load (copies/mL)")
# Save plot as pdf
# ggsave("plots/viral_load_trajectories.pdf", p_vl, width = 10, height = 6)

#Version 2 with geom_ribbon
summary_df <- p_vl$data %>%
  group_by(time) %>%
  summarise(
    mean_vl = mean(viral_load),
    lower = quantile(viral_load, 0.025),
    upper = quantile(viral_load, 0.975)
  )

p_vl.2  <- ggplot(summary_df, aes(x = time, y = mean_vl)) +
            geom_line(color = "#d95f02", size = 1) +
            geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#d95f02", alpha = 0.2) +
            labs(
              title = "Mean Viral Load with 95% Confidence Interval",
              x = "Time (days)",
              y = "Viral Load"
            ) +
            theme_minimal(base_size = 14)

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
  geom_line(aes(y = S, color = "Susceptible"), linewidth = 1.2) +
  geom_line(aes(y = I, color = "Infectious"), linewidth = 1.2) +
  geom_line(aes(y = R, color = "Recovered"), linewidth = 1.2) +
  labs(title = "",
       x = "Day",
       y = "Number of Individuals",
       color = "") +
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
    geom_line(aes(y = S, color = "Susceptible"), alpha = 0.4, linewidth = 0.6) +
    geom_line(aes(y = I, color = "Infectious"),  alpha = 0.4, linewidth = 0.6) +
    geom_line(aes(y = R, color = "Recovered"),   alpha = 0.4, linewidth = 0.6) +
    scale_color_manual(values = c(
      Susceptible = "#1b9e77",  # teal green
      Infectious  = "#d95f02",  # rich orange
      Recovered   = "#7570b3"   # muted purple
    )) +
    labs(
      title = paste("Spaghetti Plot of", N, "Simulations"),
      x = "Day",
      y = "Number of Individuals",
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

p_spaghetti <- plot_spaghetti(params, N = 100)
print(p_spaghetti <-p_spaghetti + xlim(0, 80) + labs(title="") )
# Save high resolution plot
# ggsave("plots/abm_spaghetti_plot.png", p_spaghetti, width = 10, height = 6, dpi = 300)
# pdf version
# ggsave("plots/abm_spaghetti_plot.pdf", p_spaghetti, width = 10, height = 6)

# ------------------------------------------------------------------------------
# Combine and Save Plots
# ------------------------------------------------------------------------------

combined_plot <- (p_spaghetti + p_vl ) +
  plot_annotation(
    tag_levels = 'A',
    title = ""
  )

# Show the plot
print(combined_plot)

# Save the figure
ggsave("plots/MSABM_panel_figure.pdf", combined_plot, width = 12, height = 6)

# ------------------------------------------------------------------------------
# Test the model variability for different sigma values
# ------------------------------------------------------------------------------

sigma_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                  1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2)
results_sigma <- function(sigma_values) {
  results <- data.frame()

  for( sigma in sigma_values )
  {
    params$sigma <- rep(sigma, length(params$age_groups))
    for( i in seq(50))
    {
      params$seed <- sample(1:100000, 1)  # Random seed for each run
      out <- run_abm(params)
      tmp <- data.frame(  day = 1:params$n_days,
                          S   = out$prevalence_SUSCEPTIBLE,
                          I   = out$prevalence_INFECTIOUS,
                          R   = out$prevalence_RECOVERED,
                          seed = params$seed,
                          sigma = sigma
                        )
      results <- rbind(results, tmp)
    }
  }
  return(results)
}
df <- results_sigma(sigma_values)  # Run for the first sigma value to initialize¨
# Compute cumulative infections by sigma and seed
df2 <- df %>%
  group_by(sigma, seed) %>%
  summarise(I = sum(I), .groups = 'drop') %>%
  ungroup()

# Box plot in ggplot
ggplot(df2, aes(x = factor(sigma), y = I)) +
  geom_boxplot(fill = "#d95f02", alpha = 0.7) +
  labs(
    title = "Cumulative Infections by Sigma Value",
    x = "Sigma Value",
    y = "Cumulative Infections"
  ) +
  theme_minimal(base_size = 14) +
  ylim(900000, 1200000)   # Adjust y-axis limit as needed
