# Internal helpers and C++ symbols (not for direct use).

risk_from_stats = function(n, s1, s2) ifelse(n <= 1L, 0.0, s2 - (s1 * s1) / n)

active_effect_rel_tol = function(default = 1e-4) {
  tol = getOption("xplaineff.active_effect_rel_tol", default)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0 || tol > 1) {
    return(default)
  }
  tol
}

active_effect_indices = function(objective_value_j, rel_tol = active_effect_rel_tol()) {
  values = as.numeric(objective_value_j)
  if (length(values) == 0L || rel_tol <= 0) {
    return(seq_along(values))
  }
  finite_positive = is.finite(values) & values > 0
  total = sum(values[finite_positive], na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(seq_along(values))
  }

  active = finite_positive & values > total * rel_tol
  if (any(active)) {
    return(which(active))
  }

  finite_values = values
  finite_values[!is.finite(finite_values)] = -Inf
  which.max(finite_values)
}

active_effect_names = function(objective_value_j, effect_names = names(objective_value_j),
  rel_tol = active_effect_rel_tol()) {
  idx = active_effect_indices(objective_value_j = objective_value_j, rel_tol = rel_tol)
  if (is.null(effect_names) || length(effect_names) != length(objective_value_j)) {
    return(as.character(idx))
  }
  effect_names[idx]
}

subset_ale_compact_features = function(effect, features) {
  pos = match(features, effect$feature_names)
  if (anyNA(pos)) {
    cli::cli_abort("Compact ALE effect does not contain requested features.")
  }
  out = effect
  out$feature_names = effect$feature_names[pos]
  out$d_l_mat = effect$d_l_mat[pos, , drop = FALSE]
  out$interval_idx_mat = effect$interval_idx_mat[pos, , drop = FALSE]
  out$feature_value_mat = effect$feature_value_mat[pos, , drop = FALSE]
  if (!is.null(effect$feature_types)) {
    out$feature_types = effect$feature_types[pos]
  }
  if (!is.null(effect$feature_levels)) {
    out$feature_levels = effect$feature_levels[pos]
  }
  out
}

prune_effects_for_split_search = function(Y, objective_value_j, rel_tol = active_effect_rel_tol()) {
  effect_names = if (is_ale_compact(Y)) Y$feature_names else names(Y)
  if (!is.null(effect_names) && length(effect_names) == length(objective_value_j)) {
    names(objective_value_j) = effect_names
  }
  active_idx = active_effect_indices(objective_value_j = objective_value_j, rel_tol = rel_tol)
  active_objective_value_j = objective_value_j[active_idx]
  active_names = if (is.null(effect_names)) names(active_objective_value_j) else effect_names[active_idx]
  Y_active = if (is_ale_compact(Y)) {
    subset_ale_compact_features(Y, active_names)
  } else {
    Y[active_idx]
  }

  list(
    Y = Y_active,
    objective_value_j = active_objective_value_j,
    objective_value = sum(active_objective_value_j, na.rm = TRUE),
    active_features = active_names
  )
}

assert_ale_effect_list = function(Y, var_name = "Y") {
  required_cols = c("row_id", "interval_index", "d_l", "int_n", "int_s1", "int_s2")
  checkmate::assert_list(Y, min.len = 1, .var.name = var_name)
  checkmate::assert_true(all(mlr3misc::map_lgl(Y, is.data.frame)), .var.name = var_name)
  missing_cols = mlr3misc::map(Y, function(dt) setdiff(required_cols, colnames(dt)))
  bad = Filter(function(x) length(x) > 0, missing_cols)
  if (length(bad) > 0) {
    cli::cli_abort(
      "{.arg {var_name}} elements are missing required columns: {.val {unique(unlist(bad))}}."
    )
  }
}

default_predict_fun = function(model, data) {
  if (is.function(model)) {
    return(extract_numeric_prediction(model(data), expected_n = nrow(data)))
  }
  fast_pred = fast_predict_regression_model(model, data)
  if (!is.null(fast_pred)) {
    return(extract_numeric_prediction(fast_pred, expected_n = nrow(data)))
  }
  predict_newdata_fast_dispatch(model, data)
}

make_effect_predictor = function(model, predict_fun = NULL) {
  checkmate::assert_function(predict_fun, null.ok = TRUE, .var.name = "predict_fun")
  prefer_data_table = is.null(predict_fun) && has_predict_method(model, "predict_newdata_fast")
  predict = if (is.null(predict_fun)) {
    function(newdata) default_predict_fun(model, newdata)
  } else {
    function(newdata) {
      extract_numeric_prediction(predict_fun(model, newdata), expected_n = nrow(newdata))
    }
  }
  list(
    predict = predict,
    prefer_data_table = prefer_data_table
  )
}

has_predict_method = function(model, method) {
  candidate = tryCatch(model[[method]], error = function(e) NULL)
  !is.null(candidate) && is.function(candidate)
}

# Coerce mlr3 Prediction / matrix / vector outputs to a numeric vector (one value per row).
# For multiclass probability PD/ICE on a specific class, pass a custom `predict_fun` in PdStrategy.
# When `expected_n` is set, abort if length differs (including empty vs nonempty rows).
extract_numeric_prediction = function(pred, expected_n = NULL) {
  out = if (inherits(pred, "Prediction")) {
    if (!is.null(pred$response)) {
      as.numeric(pred$response)
    } else if (!is.null(pred$prob)) {
      if (ncol(pred$prob) > 1L) {
        cli::cli_warn(
          "Multiclass model detected: using first class probability column ({colnames(pred$prob)[1L]}) as prediction.
           Pass {.arg predict_fun} to override."
        )
      }
      as.numeric(pred$prob[, 1L])
    } else {
      cli::cli_abort(
        "{.cls Prediction} object has neither {.field response} nor {.field prob}; cannot extract numeric predictions."
      )
    }
  } else if (is.list(pred) && !is.null(pred$response)) {
    as.numeric(pred$response)
  } else if (is.list(pred) && !is.null(pred$predictions)) {
    as.numeric(pred$predictions)
  } else if (is.data.frame(pred)) {
    as.numeric(pred[[1L]])
  } else if (is.matrix(pred)) {
    as.numeric(pred[, 1L])
  } else {
    as.numeric(pred)
  }

  if (!is.null(expected_n)) {
    if (length(out) != expected_n) {
      cli::cli_abort(
        "Prediction length mismatch: got {length(out)} value{?s} but 
        expected {expected_n} (one per row of {.arg newdata})."
      )
    }
  }
  out
}

# Prediction dispatch priority:
# 1) mlr3 predict_newdata_fast method (requires data.table),
# 2) mlr3 predict_newdata method,
# 3) generic stats::predict fallback.
predict_newdata_fast_dispatch = function(model, newdata) {
  n_rows = nrow(newdata)
  if (has_predict_method(model, "predict_newdata_fast")) {
    newdata_fast = if (data.table::is.data.table(newdata)) newdata else data.table::as.data.table(newdata)
    pred = model$predict_newdata_fast(newdata_fast)
    return(extract_numeric_prediction(pred, expected_n = n_rows))
  }
  if (has_predict_method(model, "predict_newdata")) {
    pred = model$predict_newdata(newdata)
    return(extract_numeric_prediction(pred, expected_n = n_rows))
  }
  pred = tryCatch(
    stats::predict(model, newdata = newdata, type = "response"),
    error = function(e1) {
      tryCatch(
        stats::predict(model, data = newdata, type = "response"),
        error = function(e2) {
          cli::cli_abort(c(
            "Prediction failed for both {.code stats::predict(..., newdata =)} 
            and {.code stats::predict(..., data =)}.",
            "i" = paste0("newdata: ", conditionMessage(e1)),
            "i" = paste0("data: ", conditionMessage(e2))
          ))
        }
      )
    }
  )
  extract_numeric_prediction(pred, expected_n = n_rows)
}

fast_predict_regression_model = function(model, newdata) {
  pred = predict_ranger_regression_fast(model, newdata)
  if (!is.null(pred)) {
    return(pred)
  }
  pred = predict_rpart_regression_fast(model, newdata)
  if (!is.null(pred)) {
    return(pred)
  }
  pred = predict_xgboost_regression_fast(model, newdata)
  if (!is.null(pred)) {
    return(pred)
  }
  NULL
}

predict_ranger_regression_fast = function(model, newdata) {
  info = extract_ranger_regression_model(model)
  if (is.null(info)) {
    return(NULL)
  }
  newdata_selected = select_newdata_features(newdata, info$feature_names)
  if (is.null(newdata_selected)) {
    return(NULL)
  }
  predict_args = list(object = info$model, data = as.data.frame(newdata_selected))
  if (!is.null(info$num_threads)) {
    predict_args$num.threads = info$num_threads
  }
  pred = tryCatch(
    do.call(stats::predict, predict_args)$predictions,
    error = function(e) NULL
  )
  if (is.null(pred)) {
    return(NULL)
  }
  as.numeric(pred)
}

extract_ranger_regression_model = function(model) {
  if (inherits(model, "ranger") && identical(model$treetype, "Regression")) {
    feature_names = tryCatch(model$forest$independent.variable.names, error = function(e) NULL)
    return(list(model = model, feature_names = feature_names, num_threads = extract_ranger_num_threads(model)))
  }
  if (inherits(model, "LearnerRegr")) {
    ranger_model = extract_mlr3_native_model(model)
    if (inherits(ranger_model, "ranger") && identical(ranger_model$treetype, "Regression")) {
      return(list(
        model = ranger_model,
        feature_names = extract_mlr3_feature_names(model),
        num_threads = extract_ranger_num_threads(ranger_model, learner = model)
      ))
    }
  }
  NULL
}

extract_ranger_num_threads = function(model, learner = NULL) {
  num_threads = tryCatch(model$call$num.threads, error = function(e) NULL)
  if (is.numeric(num_threads) && length(num_threads) == 1L && is.finite(num_threads)) {
    return(as.integer(num_threads))
  }
  if (!is.null(learner)) {
    num_threads = tryCatch(learner$param_set$values[["num.threads"]], error = function(e) NULL)
    if (is.numeric(num_threads) && length(num_threads) == 1L && is.finite(num_threads)) {
      return(as.integer(num_threads))
    }
  }
  NULL
}

predict_rpart_regression_fast = function(model, newdata) {
  info = extract_rpart_regression_model(model)
  if (is.null(info)) {
    return(NULL)
  }
  newdata_selected = select_newdata_features(newdata, info$feature_names)
  if (is.null(newdata_selected)) {
    return(NULL)
  }
  pred = tryCatch(
    stats::predict(info$model, newdata = as.data.frame(newdata_selected), type = "vector"),
    error = function(e) NULL
  )
  if (is.null(pred)) {
    return(NULL)
  }
  as.numeric(pred)
}

extract_rpart_regression_model = function(model) {
  if (inherits(model, "rpart") && identical(model$method, "anova")) {
    return(list(model = model, feature_names = NULL))
  }
  if (inherits(model, "LearnerRegr")) {
    rpart_model = extract_mlr3_native_model(model)
    if (inherits(rpart_model, "rpart") && identical(rpart_model$method, "anova")) {
      return(list(model = rpart_model, feature_names = extract_mlr3_feature_names(model)))
    }
  }
  NULL
}

predict_xgboost_regression_fast = function(model, newdata) {
  info = extract_xgboost_regression_model(model)
  if (is.null(info)) {
    return(NULL)
  }
  x = numeric_matrix_for_prediction(newdata, info$feature_names)
  if (is.null(x)) {
    return(NULL)
  }
  pred = tryCatch(
    stats::predict(info$model, newdata = x),
    error = function(e) NULL
  )
  if (is.null(pred)) {
    return(NULL)
  }
  as.numeric(pred)
}

extract_xgboost_regression_model = function(model) {
  if (inherits(model, "LearnerRegr")) {
    booster = extract_mlr3_native_model(model)
    if (inherits(booster, "xgb.Booster")) {
      return(list(model = booster, feature_names = extract_mlr3_feature_names(model)))
    }
  }
  if (inherits(model, "xgb.Booster") && isTRUE(xgboost_booster_is_regression(model))) {
    return(list(model = model, feature_names = NULL))
  }
  NULL
}

xgboost_booster_is_regression = function(model) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    return(FALSE)
  }
  xgb_config = tryCatch(get("xgb.config", envir = asNamespace("xgboost"), mode = "function"), error = function(e) NULL)
  if (is.null(xgb_config)) {
    return(FALSE)
  }
  config = tryCatch(xgb_config(model), error = function(e) NULL)
  if (is.null(config)) {
    return(FALSE)
  }
  objective = config$learner$learner_train_param$objective
  if (is.null(objective)) {
    objective = config$learner$objective$name
  }
  is.character(objective) && length(objective) == 1L && grepl("^reg:", objective)
}

extract_mlr3_native_model = function(model) {
  native_model = tryCatch(model$native_model, error = function(e) NULL)
  if (!is.null(native_model)) {
    return(native_model)
  }
  learner_model = tryCatch(model$model, error = function(e) NULL)
  if (is.list(learner_model) && !is.null(learner_model$model)) {
    return(learner_model$model)
  }
  learner_model
}

extract_mlr3_feature_names = function(model) {
  state = tryCatch(model$state, error = function(e) NULL)
  feature_names = state$feature_names
  if (is.character(feature_names) && length(feature_names)) {
    feature_names
  } else {
    NULL
  }
}

select_newdata_features = function(newdata, feature_names = NULL) {
  if (is.null(feature_names)) {
    return(newdata)
  }
  if (identical(names(newdata), feature_names)) {
    return(newdata)
  }
  if (!all(feature_names %in% names(newdata))) {
    return(NULL)
  }
  if (data.table::is.data.table(newdata)) {
    newdata[, feature_names, with = FALSE]
  } else {
    newdata[, feature_names, drop = FALSE]
  }
}

numeric_matrix_for_prediction = function(newdata, feature_names = NULL) {
  selected = select_newdata_features(newdata, feature_names)
  if (is.null(selected)) {
    return(NULL)
  }
  supported = vapply(selected, function(x) is.numeric(x) || is.integer(x) || is.logical(x), logical(1L))
  if (!all(supported)) {
    return(NULL)
  }
  x = as.matrix(selected)
  storage.mode(x) = "double"
  x
}

# Wrap each line of a node label (split by \\n) for ggraph tree plots.
wrap_tree_label = function(text, width = 34L) {
  if (length(text) == 0L) {
    return("")
  }
  text = text[[1L]]
  if (is.na(text) || !nzchar(text)) {
    return(text)
  }
  lines = strsplit(text, "\n", fixed = TRUE)[[1L]]
  wrapped = vapply(lines, function(line) {
    if (!nzchar(line)) {
      return(line)
    }
    paste(strwrap(line, width = width, exdent = 2L), collapse = "\n")
  }, character(1L))
  paste(wrapped, collapse = "\n")
}

#' Internal C++ helpers and package symbols
#'
#' Functions and symbols used internally by the package. Not intended for direct use.
#'
#' @name xplaineff_internal
#' @aliases ale_sweep_cpp calculate_ale_matrix calculate_ale_heterogeneity_list_cpp
#'   calculate_ale_heterogeneity_single_cpp re_mean_center_ice_cpp
#'   search_best_split_cpp default_predict_fun make_effect_predictor predict_newdata_fast_dispatch
#'   fast_predict_regression_model predict_ranger_regression_fast extract_ranger_regression_model
#'   extract_ranger_num_threads predict_rpart_regression_fast
#'   extract_rpart_regression_model predict_xgboost_regression_fast
#'   extract_xgboost_regression_model xgboost_booster_is_regression
#'   extract_mlr3_native_model extract_mlr3_feature_names select_newdata_features numeric_matrix_for_prediction
#'   extract_numeric_prediction has_predict_method cpp_pd_stack_newdata risk_from_stats assert_ale_effect_list
#'   active_effect_rel_tol active_effect_indices active_effect_names subset_ale_compact_features
#'   prune_effects_for_split_search d_l interval_index level x x_grid x_left x_right y
#' @keywords internal
NULL
