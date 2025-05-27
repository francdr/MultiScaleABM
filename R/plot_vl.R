#' Plot viral load trajectory for one individual
#'
#' @param age_group Index of the age group (starting at 0). If NULL, one is sampled.
#' @param params List of parameters (must include dt, VL_days, k, gamma, sigma, age_groups)
#' @param I0 Initial viral load (default 0.1)
#' @param baseline Baseline viral load (default 0.0)
#' @export
plot_viral_load_trajectory <- function(age_group = NULL, params)
{
  if (is.null(age_group))
  {
    age_group <- sample(seq_along(params$age_groups) - 1, 1)
    message("Sampled age group: ", age_group)
  }

  vl <- test_viral_load_trajectory(
    age_group = age_group,
    dt = params$dt,
    VL_days = params$VL_days,
    I0 = params$V0,
    baseline = params$dV0,
    k = params$k,
    gamma = params$gamma,
    sigma = params$sigma
  )

  plot(seq_along(vl), vl, type = "l", lwd = 2, col = "blue",
       xlab = "Days", ylab = "Relative Infectiousness",
       main = paste("Viral Load Trajectory — Age group:", age_group))
}
