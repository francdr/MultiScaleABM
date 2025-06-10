# ==============================================================================
# Benchmarking Script for MultiScaleABM
# Description: Compares performance and epidemic dynamics with and without
#              parallelization, and measures scaling with population size.
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
  n_agents = 1e5,
  n_days = 180,
  contact_matrix = contact_matrix,
  age_groups = age_distribution,
  k = rep(0.06, length(age_distribution)),
  gamma = rep(0.3, length(age_distribution)),
  sigma = rep(0.36, length(age_distribution)),
  beta = 0.08,
  dt = 0.05,
  VL_days = 20,
  V0 = log(0.06),
  dV0 = 5,
  n_initial_infected = 10,
  recovery_threshold = 0.01,
  parallel = TRUE,
  n_threads = 12,
  verbose = T
)


# ------------------------------------------------------------------------------
# Compare runtime: parallel vs serial
# ------------------------------------------------------------------------------

cat("\n=== Runtime Comparison: Parallel vs Serial ===\n")
params$seed <- 42

params$parallel <- TRUE
start_time <- Sys.time()
d<-run_abm(params)
end_time <- Sys.time()
cat("Parallel time: ", end_time - start_time, "\n")

params$parallel <- FALSE
start_time <- Sys.time()
d<-run_abm(params)
end_time <- Sys.time()
cat("Serial time:   ", end_time - start_time, "\n")

# ------------------------------------------------------------------------------
# Compare epidemic curves over N stochastic simulations
# ------------------------------------------------------------------------------

compare_parallel_vs_serial <- function(params, N = 30, plot = TRUE) {
  run_once <- function(parallel, seed) {
    params$parallel <- parallel
    params$seed <- seed
    run_abm(params)
  }

  n_days <- params$n_days
  results_serial <- matrix(0, nrow = N, ncol = n_days)
  results_parallel <- matrix(0, nrow = N, ncol = n_days)

  for (i in 1:N) {
    seed <- 1000 + i
    results_serial[i, ]   <- run_once(parallel = FALSE, seed = seed)$daily_infected
    results_parallel[i, ] <- run_once(parallel = TRUE,  seed = seed)$daily_infected
  }

  mean_serial   <- colMeans(results_serial)
  mean_parallel <- colMeans(results_parallel)

  if (plot) {
    matplot(t(rbind(mean_serial, mean_parallel)),
            type = "l", lty = 1, col = c("blue", "red"),
            ylab = "Mean Daily Infected", xlab = "Day",
            main = paste("Serial vs Parallel (N =", N, ")"))
    legend("topright", legend = c("Serial", "Parallel"),
           col = c("blue", "red"), lty = 1, bty = "n")
  }

  return(list(
    mean_serial = mean_serial,
    mean_parallel = mean_parallel,
    raw_serial = results_serial,
    raw_parallel = results_parallel
  ))
}

cat("\n=== Comparing Epidemic Curves ===\n")
comparison <- compare_parallel_vs_serial(params, N = 50)

# ------------------------------------------------------------------------------
# Evaluate scaling: run time vs population size
# ------------------------------------------------------------------------------

run_time_vs_population <- function(pop_sizes, params, n_days = 180, PLOT = TRUE) {
  run_times <- numeric(length(pop_sizes))

  for (i in seq_along(pop_sizes)) {
    cat(sprintf("Simulating %d agents...\n", pop_sizes[i]))
    params$n_agents <- pop_sizes[i]
    start_time <- Sys.time()
    run_abm(params)
    end_time <- Sys.time()
    run_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cat(sprintf("  Time: %.2f seconds\n", run_times[i]))
  }

  results <- data.frame(
    population_size = pop_sizes,
    run_time = run_times
  )

  p <- ggplot(results, aes(x = population_size, y = run_time)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Run Time vs Population Size",
         x = "Population Size (log10)", y = "Run Time (seconds, log10)") +
    theme_minimal()

  if (PLOT) print(p)

  return(list(results = results, plot = p))
}

cat("\n=== Scaling Benchmark: Varying Population Size ===\n")
pop_sizes <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
scaling <- run_time_vs_population(pop_sizes, params, n_days = 180)

# ------------------------------------------------------------------------------
# Optional: Save results (uncomment to enable)
# ------------------------------------------------------------------------------

# dir.create("benchmark_results", showWarnings = FALSE)
# write.csv(scaling$results, "benchmark_results/runtime_vs_population.csv", row.names = FALSE)
# ggsave("benchmark_results/runtime_scaling_plot.png", scaling$plot, width = 8, height = 5)
