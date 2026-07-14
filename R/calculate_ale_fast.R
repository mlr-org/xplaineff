#' Calculate ALE via C++ kernels
#'
#' Fast ALE path that keeps model prediction in R, while using C++ kernels for
#' interval indexing and interval-wise aggregation.
#' Numeric and categorical branches stack lower/upper (or plus/minus) configurations
#' into \code{2n} rows and call the predictor once per feature, matching the batched
#' \code{predict_newdata_fast} pattern used in \code{xplaineff::ale}.
#'
#' @param model (`any`) \cr
#'   Fitted model with predict interface.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param feature_set (`character()`) \cr
#'   Features to compute ALE for.
#' @param target_feature_name (`character(1)`) \cr
#'   Target variable name.
#' @param n_intervals (`integer(1)`) \cr
#'   Equal-frequency intervals for numeric features.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   \code{function(model, data)} returning predictions; \code{NULL} = default.
#' @return (`list()`) \cr
#'   Named list of ALE data.tables, same schema as \code{calculate_ale()}.
#' @keywords internal
calculate_ale_fast = function(
  model, data, feature_set, target_feature_name, n_intervals = 10, predict_fun = NULL
) {
  X = if (data.table::is.data.table(data)) {
    data[, setdiff(colnames(data), target_feature_name), with = FALSE]
  } else {
    data[, setdiff(colnames(data), target_feature_name), drop = FALSE]
  }
  x_features_dt = data.table::as.data.table(X)

  if (is.null(predict_fun)) {
    predict_fun = default_predict_fun
  }
  n_intervals = as.integer(n_intervals)
  predictor = make_predictor(model = model, predict_fun = predict_fun)

  # Pre-allocate stacked once: avoids p x 2 x (n*p) copies inside the feature loop.
  # Each ale_feature call modifies only the current feature's column in-place, then restores it,
  # so background features always reflect observed X values (required by ALE).
  n_rows = nrow(X)
  stacked = data.table::rbindlist(list(X, X), use.names = TRUE)
  idx_lower = seq_len(n_rows)
  idx_upper = seq.int(n_rows + 1L, 2L * n_rows)

  mlr3misc::map(setNames(nm = feature_set), function(feat) {
    ale_feature(
      data = data,
      X = X,
      stacked = stacked,
      idx_lower = idx_lower,
      idx_upper = idx_upper,
      feature = feat,
      n_intervals = n_intervals,
      predictor = predictor
    )
  })
}

calculate_ale_fast_compact = function(
  model, data, feature_set, target_feature_name, n_intervals = 10, predict_fun = NULL
) {
  X = if (data.table::is.data.table(data)) {
    data[, setdiff(colnames(data), target_feature_name), with = FALSE]
  } else {
    data[, setdiff(colnames(data), target_feature_name), drop = FALSE]
  }
  x_features_dt = data.table::as.data.table(X)
  feature_cols = x_features_dt[, feature_set, with = FALSE]
  is_supported_col = vapply(feature_cols, function(x) is.numeric(x) || is.integer(x), logical(1L))
  if (!all(is_supported_col)) {
    return(NULL)
  }

  if (is.null(predict_fun)) {
    predict_fun = default_predict_fun
  }
  n_intervals = as.integer(n_intervals)
  predictor = make_predictor(model = model, predict_fun = predict_fun)

  n_rows = nrow(X)
  p = length(feature_set)
  stacked = data.table::rbindlist(list(X, X), use.names = TRUE)
  idx_lower = seq_len(n_rows)
  idx_upper = seq.int(n_rows + 1L, 2L * n_rows)

  d_l_mat = matrix(0.0, nrow = p, ncol = n_rows, dimnames = list(feature_set, NULL))
  interval_idx_mat = matrix(1L, nrow = p, ncol = n_rows, dimnames = list(feature_set, NULL))
  feature_value_mat = matrix(NA_real_, nrow = p, ncol = n_rows, dimnames = list(feature_set, NULL))

  for (j in seq_along(feature_set)) {
    feat = feature_set[[j]]
    x_num = data[[feat]]
    feature_value_mat[j, ] = as.numeric(x_num)
    if (length(unique(na.omit(x_num))) <= 1L) {
      next
    }

    prep = cpp_ale_numeric_prepare(x = as.numeric(x_num), n_intervals = n_intervals)
    if (isTRUE(prep$zero_effect)) {
      next
    }

    original = data.table::copy(stacked[[feat]])
    if (is.integer(original)) {
      data.table::set(stacked, j = feat, value = as.numeric(original))
    }
    data.table::set(stacked, i = idx_lower, j = feat, value = prep$x_left)
    data.table::set(stacked, i = idx_upper, j = feat, value = prep$x_right)
    preds_all = predictor$predict(stacked) + 0
    data.table::set(stacked, j = feat, value = original)

    d_l = preds_all[idx_upper] - preds_all[idx_lower]
    d_l[!is.finite(d_l)] = 0.0
    interval_index = as.integer(prep$interval_index)
    interval_index[is.na(interval_index) | interval_index < 1L] = 1L
    d_l_mat[j, ] = d_l
    interval_idx_mat[j, ] = interval_index
  }

  structure(
    list(
      feature_names = feature_set,
      d_l_mat = d_l_mat,
      interval_idx_mat = interval_idx_mat,
      feature_value_mat = feature_value_mat
    ),
    class = "xplaineff_ale_compact"
  )
}


#' Fast ALE for a single feature.
#'
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param X (`data.table()`) \cr
#'   Features (excl. target). Never modified.
#' @param stacked (`data.table()`) \cr
#'   Pre-allocated 2n-row matrix shared across features; modified in-place per call.
#' @param idx_lower (`integer()`) \cr
#'   Row indices for the lower/plus half of stacked (1..n).
#' @param idx_upper (`integer()`) \cr
#'   Row indices for the upper/minus half of stacked ((n+1)..2n).
#' @param feature (`character(1)`) \cr
#'   Feature name.
#' @param n_intervals (`integer(1)`) \cr
#'   Number of intervals.
#' @param predictor (`list()`) \cr
#'   Prediction function wrapper from \code{make_predictor}.
#'
#' @return (`data.table()`) \cr
#'   ALE data with \code{row_id}, \code{feat_val}, \code{d_l}, \code{interval_index}, etc.
#' @keywords internal
ale_feature = function(data, X, stacked, idx_lower, idx_upper, feature, n_intervals = 10L, predictor) {
  if (is.factor(data[[feature]])) {
    ale_categorical(data = data, X = X, stacked = stacked,
      idx_lower = idx_lower, idx_upper = idx_upper, feature = feature, predictor = predictor)
  } else {
    ale_numeric(data = data, X = X, stacked = stacked,
      idx_lower = idx_lower, idx_upper = idx_upper,
      feature = feature, n_intervals = n_intervals, predictor = predictor)
  }
}


ale_numeric = function(data, X, stacked, idx_lower, idx_upper, feature, n_intervals = 10L, predictor) {
  x_num = data[[feature]]
  if (length(unique(na.omit(x_num))) <= 1L) {
    return(ale_zero(feat_val = x_num))
  }

  prep = cpp_ale_numeric_prepare(x = x_num, n_intervals = n_intervals)
  if (isTRUE(prep$zero_effect)) {
    return(ale_zero(feat_val = x_num))
  }

  original = data.table::copy(stacked[[feature]])
  if (is.integer(original)) {
    data.table::set(stacked, j = feature, value = as.numeric(original))
  }
  data.table::set(stacked, i = idx_lower, j = feature, value = prep$x_left)
  data.table::set(stacked, i = idx_upper, j = feature, value = prep$x_right)
  preds_all = predictor$predict(stacked) + 0
  data.table::set(stacked, j = feature, value = original)

  data.table::as.data.table(cpp_ale_numeric_effect_table(
    feat_val = as.numeric(x_num),
    x_left = prep$x_left,
    x_right = prep$x_right,
    interval_index = prep$interval_index,
    preds_lower = preds_all[idx_lower],
    preds_upper = preds_all[idx_upper]
  ))
}

ale_categorical = function(data, X, stacked, idx_lower, idx_upper, feature, predictor) {
  x_cat = droplevels(data[[feature]])
  k = nlevels(x_cat)
  if (k <= 1L) {
    return(ale_zero(feat_val = x_cat, interval_index = as.integer(x_cat)))
  }

  levels_id = as.integer(x_cat)
  levels_orig = levels(x_cat)
  prep = cpp_ale_categorical_prepare(levels_id = levels_id, n_levels = k)
  original = data.table::copy(stacked[[feature]])
  data.table::set(stacked, i = idx_lower, j = feature,
    value = factor(levels_orig[prep$right_id], levels = levels_orig))
  data.table::set(stacked, i = idx_upper, j = feature,
    value = factor(levels_orig[prep$left_id], levels = levels_orig))
  preds_all = predictor$predict(stacked) + 0
  data.table::set(stacked, j = feature, value = original)

  out = data.table::as.data.table(cpp_ale_categorical_effect_table(
    feat_val = levels_id,
    x_left = prep$left_id,
    x_right = prep$right_id,
    interval_index = levels_id,
    y_hat_plus = preds_all[idx_lower],
    y_hat_neg = preds_all[idx_upper]
  ))
  out$feat_val = factor(levels_orig[out$feat_val], levels = levels_orig)
  out$x_left = factor(levels_orig[out$x_left], levels = levels_orig)
  out$x_right = factor(levels_orig[out$x_right], levels = levels_orig)
  out
}

ale_zero = function(feat_val, interval_index = 1L) {
  n_rows = length(feat_val)
  interval_index = rep(as.integer(interval_index), length.out = n_rows)
  data.table::data.table(
    row_id = seq_len(n_rows),
    feat_val = feat_val,
    x_left = feat_val,
    x_right = feat_val,
    d_l = rep(0, n_rows),
    interval_index = interval_index,
    int_n = rep(n_rows, n_rows),
    int_s1 = rep(0, n_rows),
    int_s2 = rep(0, n_rows)
  )
}

make_predictor = function(model, predict_fun) {
  if (identical(predict_fun, default_predict_fun)) {
    return(list(predict = function(newdata) default_predict_fun(model, newdata)))
  }
  list(predict = function(newdata) {
    extract_numeric_prediction(predict_fun(model, newdata), expected_n = nrow(newdata))
  })
}
