# Internal helpers for categorical split handling.

is_ale_ordered_categorical_split = function(is_categorical, strategy) {
  if (!isTRUE(is_categorical) || is.null(strategy) || !identical(strategy$name, "ale")) {
    return(FALSE)
  }
  categorical_split = tryCatch(strategy$categorical_split, error = function(e) NULL)
  is.null(categorical_split) || identical(categorical_split, "ordered_prefix")
}

format_split_level_set = function(levels) {
  if (is.null(levels) || !length(levels)) {
    return("{}")
  }
  paste0("{", paste(as.character(levels), collapse = ", "), "}")
}

format_categorical_split_condition = function(feature, levels) {
  paste0(feature, " in ", format_split_level_set(levels))
}

categorical_split_groups = function(x, split_levels) {
  checkmate::assert_factor(x, .var.name = "x")
  levels_x = levels(x)
  left_levels = unique(as.character(split_levels))
  if (!length(left_levels) || anyNA(left_levels) || any(!(left_levels %in% levels_x))) {
    cli::cli_abort("Categorical split levels must be non-missing levels of {.arg x}.")
  }
  list(left_levels = left_levels, right_levels = setdiff(levels_x, left_levels))
}

categorical_left_mask = function(x, split_levels) {
  mask = as.character(x) %in% as.character(split_levels)
  mask[is.na(x)] = NA
  mask
}

ordered_categorical_split_groups = function(x, split_value) {
  checkmate::assert_factor(x, .var.name = "x")
  levels_x = levels(x)
  split_level_id = match(as.character(split_value), levels_x)
  if (is.na(split_level_id)) {
    cli::cli_abort("Categorical split value {.val {split_value}} is not a level of {.arg x}.")
  }
  left_levels = levels_x[seq_len(split_level_id)]
  right_levels = if (split_level_id < length(levels_x)) {
    levels_x[(split_level_id + 1L):length(levels_x)]
  } else {
    character(0)
  }
  list(left_levels = left_levels, right_levels = right_levels, split_level_id = split_level_id)
}

one_vs_rest_categorical_split_groups = function(x, split_value) {
  categorical_split_groups(x, split_value)
}

ordered_categorical_left_mask = function(x, split_value) {
  groups = ordered_categorical_split_groups(x, split_value)
  as.integer(x) <= groups$split_level_id
}
