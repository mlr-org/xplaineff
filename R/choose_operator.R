#' Return comparison operator for child relative to parent.
#'
#' Given parent_node and current_node: if current is left child returns \code{<=} (numeric) or \code{=} (categorical);
#' if right returns \code{>} or \code{!=}.
#'
#' @param parent_node (`Node`) \cr
#'   Node with \code{children} and \code{split$value}.
#' @param current_node (`Node`) \cr
#'   One of the children.
#' @return (`character(1)`) \cr
#'   \code{"<="}, \code{">"}, \code{"="}, or \code{"!="}.
#' @keywords internal
choose_operator = function(parent_node, current_node) {
  op = NA_character_
  if (!is.null(parent_node) && !is.null(parent_node$children)) {
    if (!is.null(parent_node$children$left_child) &&
        identical(parent_node$children$left_child, current_node)) {
      op = if (is.numeric(parent_node$split$value)) "<=" else "="
    } else if (!is.null(parent_node$children$right_child) &&
        identical(parent_node$children$right_child, current_node)) {
      op = if (is.numeric(parent_node$split$value)) ">" else "!="
    }
  }
  op
}
