#' Node Transform ALE
#'
#' Subsets ALE effect data to the current node's row indices and recomputes
#' per-interval statistics. When \code{is_child} is \code{TRUE},
#' forces \code{d_l = 0} for any feature whose values are constant in this node
#' (single unique value).
#'
#' @param Y (`list()`) \cr
#'   ALE effect data per feature.
#' @param idx (`integer()`) \cr
#'   Sample indices in the current node.
#' @param is_child (`logical(1)`) \cr
#'   Whether the current node is a child node; \code{FALSE} skips constant-feature zeroing.
#'
#' @return (`list()`) \cr
#'   Transformed ALE effects per feature.
#'
node_transform_ale = function(Y, idx, is_child = FALSE) {
  idx = as.integer(idx)
  if (!is_child && is_full_ale_node(Y, idx)) {
    return(Y)
  }

  y_subset = lapply(names(Y), function(feat) {
    y_j = Y[[feat]]
    refresh_ale_interval_stats(subset_ale_rows(y_j, idx))
  })
  names(y_subset) = names(Y)
  if (is_child) {
    y_processed = lapply(names(y_subset), function(feat) {
      y_j = y_subset[[feat]]
      # Zero out d_l for constant features in this node (ALE undefined when all values equal)
      if (length(unique(y_j$feat_val)) == 1) {
        y_j$d_l = 0
        y_j = refresh_ale_interval_stats(y_j)
      }
      y_j
    })
    names(y_processed) = names(y_subset)
    y_processed
  } else {
    y_subset
  }
}

is_full_ale_node = function(Y, idx) {
  if (length(Y) == 0L) {
    return(FALSE)
  }
  n_rows = nrow(Y[[1L]])
  length(idx) == n_rows && identical(idx, seq_len(n_rows))
}

subset_ale_rows = function(y_j, idx) {
  row_id = y_j$row_id
  if (!anyNA(idx) && length(row_id) == nrow(y_j) && length(idx) > 0L &&
      min(idx) >= 1L && max(idx) <= nrow(y_j) &&
      identical(as.integer(row_id), seq_len(nrow(y_j)))) {
    return(y_j[idx])
  }
  pos = match(idx, row_id)
  y_j[pos[!is.na(pos)]]
}

refresh_ale_interval_stats = function(y_j) {
  if (nrow(y_j) == 0L) {
    return(y_j)
  }
  interval_index = as.integer(y_j$interval_index)
  d_l = y_j$d_l
  use_fast_stats = !anyNA(interval_index) && all(interval_index >= 1L) && all(is.finite(d_l))
  if (use_fast_stats) {
    agg = cpp_ale_interval_aggregate(d_l, interval_index)
    y_j[, `:=`(
      int_n = agg$interval_n[interval_index],
      int_s1 = agg$interval_s1[interval_index],
      int_s2 = agg$interval_s2[interval_index]
    )]
    return(y_j)
  }
  y_j[, `:=`(
    int_n  = .N,
    int_s1 = sum(d_l, na.rm = TRUE),
    int_s2 = sum(d_l^2, na.rm = TRUE)
  ), by = interval_index]
  y_j
}
