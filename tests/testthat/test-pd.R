test_that("compute_ice cpp matches r (numeric focal feature)", {
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(2L)
  n = 45L
  d = data.frame(x1 = runif(n), x2 = runif(n), x3 = rnorm(n))
  fit = stats::lm(x3 ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = seq(0.05, 0.95, length.out = 11L)
  ice_r = xplaineff:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("compute_ice cpp matches r (factor focal feature)", {
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1L),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(3L)
  n = 40L
  d = data.frame(
    x1 = factor(sample(c("low", "high"), n, replace = TRUE)),
    x2 = runif(n),
    y = rnorm(n)
  )
  fit = stats::lm(y ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = levels(d$x1)
  ice_r = xplaineff:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("re_mean_center_ice_cpp centers grid columns and preserves metadata", {
  tryCatch(
    xplaineff:::re_mean_center_ice_cpp(
      Y = list(x = matrix(1, nrow = 1L, dimnames = list(NULL, "a"))),
      grid = list(x = "a"),
      idx = 1L
    ),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ re_mean_center_ice not loaded")
      }
      stop(e)
    }
  )

  mat = matrix(
    c(1, 2, 10, NA, 4, 20, 7, NA, 30),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("g1", "g2", "meta"))
  )

  result = xplaineff:::re_mean_center_ice_cpp(
    Y = list(x = mat),
    grid = list(x = c("g1", NA_character_, "g2")),
    idx = c(1L, 2L, 4L)
  )

  testthat::expect_named(result, "x")
  centered = result$x
  testthat::expect_equal(dim(centered), c(3L, 3L))
  testthat::expect_equal(colnames(centered), c("g1", "g2", "meta"))
  testthat::expect_equal(centered[1L, ], c(g1 = -0.5, g2 = 0.5, meta = NA_real_))
  testthat::expect_equal(centered[2L, ], c(g1 = NA_real_, g2 = 0, meta = NA_real_))
  testthat::expect_true(all(is.na(centered[3L, ])))
})

test_that("PdStrategy reuses full-root centered matrices", {
  mat = matrix(
    c(1, -1, 2, -2, 3, -3),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )
  Y = list(x = mat)
  grid = list(x = c("a", "b"))
  strategy = PdStrategy$new()

  result = strategy$node_transform(Y = Y, idx = 1:3, grid = grid)

  testthat::expect_identical(result, Y)
})

test_that("PdStrategy only skips full-grid centering for prepared centered PD matrices", {
  strategy = PdStrategy$new()
  raw_y = list(x = matrix(c(0, 2, 10, 14), nrow = 2L, byrow = TRUE, dimnames = list(NULL, c("a", "b"))))
  raw_result = strategy$node_transform(Y = raw_y, idx = 1L, grid = list(x = c("a", "b")))

  testthat::expect_equal(unname(raw_result$x), matrix(c(-1, 1), nrow = 1L))

  centered_y = list(
    x = matrix(c(-1, 0, 1, -2, 0, 2), nrow = 2L, byrow = TRUE, dimnames = list(NULL, c("a", "b", "c"))),
    z = matrix(c(-3, 3, -4, 4), nrow = 2L, byrow = TRUE, dimnames = list(NULL, c("l", "r")))
  )
  attr(centered_y, "xplaineff_pd_centered") = TRUE

  full_grid_result = strategy$node_transform(
    Y = centered_y,
    idx = c(2L, 1L),
    grid = list(x = c("a", "b", "c"), z = c("l", "r"))
  )
  testthat::expect_equal(full_grid_result$x, centered_y$x[c(2L, 1L), , drop = FALSE])
  testthat::expect_equal(full_grid_result$z, centered_y$z[c(2L, 1L), , drop = FALSE])

  restricted_grid_result = strategy$node_transform(
    Y = centered_y,
    idx = 1:2,
    grid = list(x = c("a", "b"), z = c("l", "r"))
  )
  testthat::expect_equal(
    unname(restricted_grid_result$x),
    matrix(c(-0.5, 0.5, NA, -1, 1, NA), nrow = 2L, byrow = TRUE)
  )
  testthat::expect_equal(restricted_grid_result$z, centered_y$z)
})

test_that("compute_ice cpp matches r (ranger native model, data= predict)", {
  testthat::skip_if_not_installed("ranger")
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(4L)
  n = 50L
  d = data.frame(x1 = runif(n), x2 = runif(n), y = rnorm(n))
  fit = ranger::ranger(
    data = d,
    dependent.variable.name = "y",
    num.trees = 50L,
    num.threads = 1L,
    seed = 4L
  )
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = seq(0.1, 0.9, length.out = 9L)
  ice_r = xplaineff:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("compute_ice cpp chunks ranger grid batches without changing predictions", {
  testthat::skip_if_not_installed("ranger")
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(41L)
  n = 60L
  d = data.frame(x1 = runif(n), x2 = rnorm(n), x3 = rnorm(n), y = rnorm(n))
  fit = ranger::ranger(
    y ~ .,
    data = d,
    num.trees = 30L,
    num.threads = 1L,
    seed = 41L
  )
  x_only = d[, c("x1", "x2", "x3"), drop = FALSE]
  grid = seq(0.1, 0.9, length.out = 9L)

  withr::local_options(list(xplaineff.pd.ranger_grid_chunk_size = NULL))
  unchunked = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  withr::local_options(list(xplaineff.pd.ranger_grid_chunk_size = 2L))
  chunked = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)

  testthat::expect_equal(chunked, unchunked, tolerance = 1e-10)
})

test_that("ranger fast PD path matches default prediction for native regression forests", {
  testthat::skip_if_not_installed("ranger")
  tryCatch(
    xplaineff:::ranger_pd_numeric_cpp(
      list(
        num.trees = 0,
        child.nodeIDs = list(),
        split.varIDs = list(),
        split.values = list()
      ),
      matrix(numeric(), nrow = 0L),
      integer(),
      list()
    ),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ ranger fast path not loaded")
      }
      stop(e)
    }
  )
  set.seed(8L)
  n = 80L
  dat = data.frame(x1 = runif(n), x2 = rnorm(n), x3 = runif(n))
  dat$y = sin(dat$x1) + dat$x2 * (dat$x3 > 0.5)
  fit = ranger::ranger(
    y ~ .,
    data = dat,
    num.trees = 30L,
    mtry = 3L,
    min.node.size = 2L,
    num.threads = 1L,
    seed = 8L
  )
  x_only = data.table::as.data.table(dat[, c("x1", "x2", "x3")])
  grids = list(x1 = seq(0.1, 0.9, length.out = 5L), x2 = seq(-1, 1, length.out = 4L))

  withr::local_options(list(xplaineff.pd.ranger_fast = TRUE))
  fast_info = xplaineff:::pd_ranger_fast_info(fit, x_only, names(grids), grids, predict_fun = NULL)
  fast = xplaineff:::calculate_pd_ranger_matrix(
    forest = fast_info$forest,
    x_features_dt = x_only,
    forest_features = fast_info$forest_features,
    feature_indices = fast_info$feature_indices,
    feature_set = names(grids),
    grids = grids
  )

  expected_x1 = xplaineff:::compute_ice_r(fit, x_only, "x1", grids$x1, predict_fun = NULL)
  expected_x2 = xplaineff:::compute_ice_r(fit, x_only, "x2", grids$x2, predict_fun = NULL)
  testthat::expect_equal(unname(fast$Y$x1), unname(expected_x1), tolerance = 1e-10)
  testthat::expect_equal(unname(fast$Y$x2), unname(expected_x2), tolerance = 1e-10)
})

test_that("ranger fast PD path matches default prediction for mlr3 ranger learners", {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("mlr3learners")
  testthat::skip_if_not_installed("ranger")
  tryCatch(
    xplaineff:::ranger_pd_numeric_cpp(
      list(
        num.trees = 0,
        child.nodeIDs = list(),
        split.varIDs = list(),
        split.values = list()
      ),
      matrix(numeric(), nrow = 0L),
      integer(),
      list()
    ),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ ranger fast path not loaded")
      }
      stop(e)
    }
  )
  set.seed(9L)
  n = 80L
  dat = data.frame(x1 = runif(n), x2 = rnorm(n), x3 = runif(n))
  dat$y = dat$x1 - 0.5 * dat$x2 + dat$x3
  task = mlr3::TaskRegr$new("pd_ranger_fast", backend = dat, target = "y")
  learner = mlr3::lrn("regr.ranger", num.trees = 30L, mtry = 3L, min.node.size = 2L, num.threads = 1L)
  learner$train(task)
  x_only = data.table::as.data.table(dat[, c("x1", "x2", "x3")])
  grids = list(x1 = seq(0.1, 0.9, length.out = 5L), x3 = seq(0.2, 0.8, length.out = 4L))

  withr::local_options(list(xplaineff.pd.ranger_fast = TRUE))
  fast_info = xplaineff:::pd_ranger_fast_info(learner, x_only, names(grids), grids, predict_fun = NULL)
  fast = xplaineff:::calculate_pd_ranger_matrix(
    forest = fast_info$forest,
    x_features_dt = x_only,
    forest_features = fast_info$forest_features,
    feature_indices = fast_info$feature_indices,
    feature_set = names(grids),
    grids = grids
  )

  expected_x1 = xplaineff:::compute_ice_r(learner, x_only, "x1", grids$x1, predict_fun = NULL)
  expected_x3 = xplaineff:::compute_ice_r(learner, x_only, "x3", grids$x3, predict_fun = NULL)
  testthat::expect_equal(unname(fast$Y$x1), unname(expected_x1), tolerance = 1e-10)
  testthat::expect_equal(unname(fast$Y$x3), unname(expected_x3), tolerance = 1e-10)
})

test_that("custom predict_fun PD cpp engine uses the same values as the R backend", {
  set.seed(10L)
  data = data.frame(x1 = runif(40L), x2 = rnorm(40L), x3 = runif(40L))
  data$y = data$x1 + data$x2
  pred_fun = function(model, newdata) {
    newdata$x1 - 0.5 * newdata$x2 + ifelse(newdata$x3 > 0.5, newdata$x1, 0)
  }

  cpp = xplaineff:::calculate_pd_matrix(
    model = "toy",
    data = data,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    predict_fun = pred_fun,
    n_grid = 5L,
    pd_engine = "cpp"
  )
  r = xplaineff:::calculate_pd_matrix(
    model = "toy",
    data = data,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    predict_fun = pred_fun,
    n_grid = 5L,
    pd_engine = "r"
  )

  testthat::expect_equal(cpp$Y, r$Y)
  testthat::expect_equal(cpp$grid, r$grid)
})

test_that("PdStrategy fit aborts when target column is missing from data", {
  d = data.frame(x1 = 1, y = 2)
  fit = stats::lm(y ~ x1, data = d)
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 1L, min_node_size = 1L)
  testthat::expect_error(
    tree$fit(data = d, target_feature_name = "yy", model = fit),
    regexp = "target_feature_name|yy|subset"
  )
})

test_that("AleStrategy fit aborts when target column is missing from data", {
  d = data.frame(x1 = rnorm(20), y = rnorm(20))
  fit = stats::lm(y ~ x1, data = d)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1L, min_node_size = 1L)
  testthat::expect_error(
    tree$fit(data = d, target_feature_name = "yy", model = fit, n_intervals = 3L),
    regexp = "target_feature_name|yy|subset"
  )
})

test_that("pd_feature_grid aborts for all-NA numeric predictor", {
  testthat::expect_error(
    xplaineff:::pd_feature_grid(c(NA_real_, NA_real_), 5L),
    regexp = "grid|quantiles|finite|NA"
  )
})

test_that("compute_ice cpp matches r for integer focal feature (C++ promotes stacked column to double)", {
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(5L)
  n = 30L
  d = data.frame(x1 = sample(1L:5L, n, replace = TRUE), x2 = runif(n), y = rnorm(n))
  fit = stats::lm(y ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = sort(unique(d$x1))[1:3]
  ice_r = xplaineff:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = xplaineff:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("compute_ice_r preserves fractional grid values for cached integer features", {
  data = data.table::data.table(x = 1L:4L, z = 0)
  grid = c(1, 2.5, 4)
  predict_fun = function(model, data) as.numeric(data$x)
  stacked_pd_cache = list(
    stacked = data.table::as.data.table(lapply(data, rep, times = length(grid))),
    max_g = length(grid),
    n_obs = nrow(data)
  )

  expect_warning({
    ice = xplaineff:::compute_ice_r(
      model = NULL,
      data = data,
      feature = "x",
      grid = grid,
      predict_fun = predict_fun,
      base_data_dt = data,
      stacked_pd_cache = stacked_pd_cache
    )
  },
    NA
  )

  testthat::expect_equal(ice, matrix(rep(grid, each = nrow(data)), nrow = nrow(data)))
  testthat::expect_equal(stacked_pd_cache$stacked$x, as.numeric(rep(data$x, times = length(grid))))
})

test_that("extract_numeric_prediction uses response then first prob column for mlr3 Prediction", {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("mlr3learners")
  testthat::skip_if_not_installed("ranger")
  task = mlr3::tsk("iris")
  lrn = mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = 30L, num.threads = 1L)
  lrn$train(task)
  pred = lrn$predict_newdata(iris[1:10, ])
  testthat::expect_equal(
    xplaineff:::extract_numeric_prediction(pred),
    as.numeric(pred$response)
  )
})

test_that("PdStrategy PD path can target one class prob via predict_fun", {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("mlr3learners")
  testthat::skip_if_not_installed("ranger")
  set.seed(1L)
  task = mlr3::tsk("iris")
  lrn = mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = 30L, num.threads = 1L)
  lrn$train(task)
  dat = iris[1:20, ]
  pf = function(model, newdata) {
    p = model$predict_newdata(newdata)$prob[, "versicolor"]
    as.numeric(p)
  }
  ice = xplaineff:::compute_ice_r(
    lrn, dat[, c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")],
    "Petal.Length", grid = c(1.4, 1.5), predict_fun = pf
  )
  testthat::expect_equal(nrow(ice), 20L)
  testthat::expect_equal(ncol(ice), 2L)
})

test_that("pd_feature_grid returns sorted unique levels for character x", {
  g = xplaineff:::pd_feature_grid(c("b", "a", "a", NA), n_grid = 5L)
  testthat::expect_equal(g, c("a", "b"))
})

test_that("prepare_split_data_pd infers feature_set from precomputed effect", {
  effect = list(results = list(
    x = data.frame(
      .id = rep(1:3, each = 2),
      .type = "ice",
      .feature = "x",
      .borders = rep(c(0, 1), times = 3),
      .value = seq_len(6)
    )
  ))
  data = data.frame(x = 1:3, z = 4:6, y = 7:9)

  prepared = xplaineff:::prepare_split_data_pd(
    effect = effect,
    data = data,
    target_feature_name = "y"
  )

  testthat::expect_equal(names(prepared$Y), "x")
  testthat::expect_equal(names(prepared$Z), c("x", "z"))
})

test_that("calculate_pd_matrix matches long-format PD preprocessing", {
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(6L)
  n = 35L
  dat = data.frame(x1 = runif(n), x2 = rnorm(n))
  dat$y = 1 + 2 * dat$x1 - dat$x2
  fit = stats::lm(y ~ x1 + x2, data = dat)

  effect_long = xplaineff:::calculate_pd(
    model = fit,
    data = dat,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    n_grid = 7L,
    pd_engine = "cpp"
  )
  effect_matrix = xplaineff:::calculate_pd_matrix(
    model = fit,
    data = dat,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    n_grid = 7L,
    pd_engine = "cpp"
  )

  long_prepared = xplaineff:::mean_center_ice(effect_long, feature_set = c("x1", "x2"))
  matrix_prepared = xplaineff:::mean_center_ice(effect_matrix, feature_set = c("x1", "x2"))

  testthat::expect_s3_class(effect_matrix, "xplaineff_pd_matrix")
  testthat::expect_equal(matrix_prepared$grid, long_prepared$grid)
  testthat::expect_equal(matrix_prepared$Y$x1, as.matrix(long_prepared$Y$x1), tolerance = 1e-10)
  testthat::expect_equal(matrix_prepared$Y$x2, as.matrix(long_prepared$Y$x2), tolerance = 1e-10)
})

test_that("PdStrategy model path caches matrix-native effects", {
  tryCatch(
    xplaineff:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(7L)
  n = 40L
  dat = data.frame(x1 = runif(n), x2 = rnorm(n))
  dat$y = dat$x1 + dat$x2
  fit = stats::lm(y ~ x1 + x2, data = dat)
  strat = PdStrategy$new()
  tree = GadgetTree$new(strategy = strat, n_split = 1L, min_node_size = 10L)
  tree$fit(data = dat, target_feature_name = "y", model = fit, n_grid = 5L, pd_engine = "cpp")

  testthat::expect_s3_class(strat$effect, "xplaineff_pd_matrix")
  testthat::expect_true(is.list(strat$effect$Y))
  testthat::expect_true(all(vapply(strat$effect$Y, is.matrix, TRUE)))
  testthat::expect_true(!is.null(tree$root))
})

test_that("calculate_y_range for PD omits raw target when mean_center is TRUE", {
  prepared_data = list(
    x1 = data.frame(`0` = 0, `1` = 0.1, node = 1L, check.names = FALSE)
  )
  dat = data.frame(y = c(100, 100))
  yr_mc = xplaineff:::calculate_y_range(prepared_data, dat, "y", mean_center = TRUE)
  yr_raw = xplaineff:::calculate_y_range(prepared_data, dat, "y", mean_center = FALSE)
  testthat::expect_true(yr_mc$ymax < 50)
  testthat::expect_true(yr_raw$ymax >= 100)
})
