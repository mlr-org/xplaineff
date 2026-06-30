#' Calculate ALE Heterogeneity
#'
#' @param Y (`list()` or `data.frame()`) \cr
#'   ALE effect data.
#'
#' @return (`numeric()`) \cr
#'   Heterogeneity value(s): vector per feature when Y is list, single value when Y is data.frame.
#' @keywords internal
calculate_ale_heterogeneity_cpp = function(Y) {
  # Handle both data.frame and list cases
  if (is.data.frame(Y)) {
    d_l = Y$d_l
    interval_index = Y$interval_index
    calculate_ale_heterogeneity_single_cpp(d_l, interval_index)
  } else if (is.list(Y)) {
    calculate_ale_heterogeneity_list_cpp(Y)
  } else {
    cli::cli_abort("Y must be either a data.frame or a list")
  }
}
