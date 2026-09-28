#' Convert factor to numeric.
#'
#' Given f: if all level labels parse as numbers,
#' uses as.numeric(as.character(f)); otherwise uses as.numeric(f) (level indices).
#' Returns numeric vector.
#'
#' @param f (`factor()`) \cr
#'   Factor to convert.
#' @return (`numeric()`) \cr
#'   Numeric vector.
#' @keywords internal
factor_to_numeric = function(f) {
  checkmate::assert_factor(f, .var.name = "f")
  f = droplevels(f)
  lv = levels(f)
  is_all_numeric = all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", lv))
  if (is_all_numeric) {
    as.numeric(as.character(f))
  } else {
    as.numeric(f)
  }
}
