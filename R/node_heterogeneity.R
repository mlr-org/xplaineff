#' Compute heterogeneity per feature from effect matrices
#'
#' Given Y (list of numeric matrices): for each matrix, computes sum over
#' columns of (sum of squares - (sum)^2/n) via \code{node_heterogeneity_cpp}.
#' Returns numeric vector of length \code{length(Y)}.
#'
#' @param Y (`list()`) \cr
#'   List of numeric matrices (ICE effect per feature).
#'
#' @return (`numeric()`) \cr
#'   Heterogeneity per feature, length \code{length(Y)}.
#' @keywords internal
node_heterogeneity = function(Y) {
  y_mat = mlr3misc::map(Y, as.matrix)
  node_heterogeneity_cpp(y_mat)
}

#' Total local-effect sum of squares per feature
#'
#' For each (mean-centered) ICE matrix, the total sum of squares
#' \eqn{\sum_k \sum_i Y_{ik}^2}. By the variance decomposition this equals
#' \eqn{\mathcal{R}_j + \mathcal{B}_j}, i.e. the interaction (between-curve) part plus the
#' main-effect (within-curve) part. It is the denominator of the interaction fraction
#' \eqn{q_j = \mathcal{R}_j / (\mathcal{R}_j + \mathcal{B}_j + \delta)} used by selective early
#' stopping Method 3. Using the total directly (instead of computing \eqn{\mathcal{B}_j}
#' separately) keeps the identity exact also when columns carry differing numbers of
#' non-missing entries.
#'
#' @param Y (`list()`) \cr
#'   List of numeric matrices (ICE effect per feature), mean-centered per node.
#'
#' @return (`numeric()`) \cr
#'   Total sum of squares per feature, length \code{length(Y)}, named like \code{Y}.
#' @keywords internal
total_effect_sum_of_squares = function(Y) {
  vapply(Y, function(mat) sum(as.matrix(mat)^2, na.rm = TRUE), NA_real_)
}
