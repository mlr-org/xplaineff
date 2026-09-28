#' Common preprocessing for ALE and PD split data.
#'
#' Resolves feature_set/split_feature and ensures factors.
#'
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Data.
#' @param target_feature_name (`character(1)` or `NULL`) \cr
#'   Target column; \code{NULL} = all columns are features.
#' @param feature_set (`character()` or `NULL`) \cr
#'   Features; \code{NULL} = all.
#' @param split_feature (`character()` or `NULL`) \cr
#'   Split features; \code{NULL} = all.
#' @return (`list()`) \cr
#'   \code{data}, \code{feature_set}, \code{split_feature}.
#' @keywords internal
prepare_split_data_common = function(data, target_feature_name, feature_set, split_feature) {
  all_features = if (is.null(target_feature_name)) colnames(data) else setdiff(colnames(data), target_feature_name)
  feature_set = resolve_split_features(feature_set, all_features, "Features")
  split_feature = resolve_split_features(split_feature, all_features, "Split features")
  data = ensure_factors(data, union(feature_set, split_feature))
  list(data = data, feature_set = feature_set, split_feature = split_feature)
}

#' Resolve feature names against available columns.
#'
#' @param requested (`character()` or `NULL`) \cr
#'   Requested feature names; \code{NULL} = use all.
#' @param all_features (`character()`) \cr
#'   Available feature names.
#' @param err_label (`character(1)`) \cr
#'   Label for error message (e.g. "Features", "Split features").
#'
#' @return (`character()`) \cr
#'   Resolved feature names.
#' @keywords internal
resolve_split_features = function(requested, all_features, err_label) {
  if (is.null(requested)) return(all_features)
  miss = setdiff(requested, all_features)
  if (length(miss) > 0L) {
    cli::cli_abort(c(
      "{err_label} not found in data: {paste(miss, collapse = ', ')}.",
      i = "Available: {paste(all_features, collapse = ', ')}"
    ))
  }
  requested
}

#' Ensure character columns are factors.
#'
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Data to modify.
#' @param cols (`character()`) \cr
#'   Column names to ensure as factor.
#'
#' @return (`data.frame()` or `data.table()`) \cr
#'   Copy of \code{data} with the listed character columns converted to factor.
#' @keywords internal
ensure_factors = function(data, cols) {
  for (c in cols) {
    if (is.character(data[[c]])) data[[c]] = factor(data[[c]])
  }
  data
}

#' Subset columns from data.frame or data.table.
#'
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Data.
#' @param cols (`character()`) \cr
#'   Column names.
#'
#' @return (`data.frame()` or `data.table()`) \cr
#'   Subset of columns.
#' @keywords internal
take_cols = function(data, cols) {
  if (data.table::is.data.table(data)) data[, cols, with = FALSE] else data[, cols, drop = FALSE]
}
