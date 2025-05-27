#' Run ABM Simulation
#' @useDynLib MultiScaleABM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @param params A list of simulation parameters
#' @export
run_abm <- function(params) {
  .Call('_MultiScaleABM_run_abm', params)
}
