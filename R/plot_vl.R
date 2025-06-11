#' Plot viral load trajectory for one or multiple individuals
#'
#' @param age_group Index of the age group (starting at 0). If NULL, one is sampled.
#' @param params List of parameters (must include dt, VL_days, k, gamma, sigma, age_groups)
#' @param I0 Initial viral load (default 0.1)
#' @param baseline Baseline viral load (default 0.0)
#' @param spaghetti If TRUE, plots N sample trajectories
#' @param N Number of sample trajectories to plot if spaghetti = TRUE
#' @export
plot_viral_load_trajectory <- function(age_group = NULL, params, spaghetti = FALSE, N = 50) {
  if (is.null(age_group)) {
    age_group <- sample(seq_along(params$age_groups) - 1, 1)
    message("Sampled age group: ", age_group)
  }

  n_days <- params$VL_days
  time_points <- seq(0, by = params$dt, length.out = n_days)
  time_points <- 1:params$VL_days

  if (!spaghetti) {
    vl <- test_viral_load_trajectory(
      age_group = age_group,
      dt = params$dt,
      VL_days = n_days,
      I0 = params$V0,
      baseline = params$dV0,
      k = params$k,
      gamma = params$gamma,
      sigma = params$sigma
    )

    df <- data.frame(
      time = time_points,
      viral_load = vl
    )

    p <- ggplot(df, aes(x = time, y = viral_load)) +
      geom_line(color = "#1f77b4", linewidth = 1.2) +
      labs(
        title = paste("Viral Load Trajectory – Age group:", age_group),
        x = "Days",
        y = "Relative Infectiousness"
      ) +
      theme_minimal(base_size = 14)

  } else {
    mat <- matrix(NA, nrow = n_days, ncol = N)
    for (i in 1:N) {
      mat[, i] <- test_viral_load_trajectory(
        age_group = age_group,
        dt = params$dt,
        VL_days = n_days,
        I0 = params$V0,
        baseline = params$dV0,
        k = params$k,
        gamma = params$gamma,
        sigma = params$sigma
      )
    }

    df <- data.frame(
      time = rep(time_points, times = N),
      viral_load = as.vector(mat),
      traj = rep(1:N, each = n_days)
    )

    p <- ggplot(df, aes(x = time, y = viral_load, group = traj)) +
      geom_line(color = "black", alpha = 0.2) +
      labs(
        title = paste("Spaghetti Plot –", N, "Trajectories – Age group:", age_group),
        x = "Days",
        y = "Relative Infectiousness"
      ) +
      theme_minimal(base_size = 14)
  }

  return(p)
}
