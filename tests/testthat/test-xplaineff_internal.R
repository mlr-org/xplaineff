test_that("native ranger regression fast prediction matches ranger predict", {
  skip_if_not_installed("ranger")
  set.seed(30L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = sin(data$x1) + data$x2 * (data$x3 > 0.5)
  model = ranger::ranger(
    y ~ .,
    data = data,
    num.trees = 30L,
    mtry = 3L,
    min.node.size = 2L,
    num.threads = 1L,
    seed = 30L
  )
  x = data[, c("x1", "x2", "x3")]

  info = xplaineff:::extract_ranger_regression_model(model)
  expect_s3_class(info$model, "ranger")
  expected = as.numeric(predict(model, data = x, num.threads = 1L)$predictions)

  withr::local_options(list(xplaineff.ranger.num_threads = 1L))
  expect_equal(xplaineff:::predict_ranger_regression_fast(model, x), expected)
  expect_equal(xplaineff:::default_predict_fun(model, x), expected)
})

test_that("select_newdata_features keeps already aligned feature data", {
  x = data.frame(x1 = 1:3, x2 = 4:6)

  expect_identical(xplaineff:::select_newdata_features(x, names(x)), x)
  expect_identical(xplaineff:::select_newdata_features(x, rev(names(x))), x[rev(names(x))])
})

test_that("mlr3 regr.ranger fast prediction uses the trained ranger model", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")
  set.seed(30L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = sin(data$x1) + data$x2 * (data$x3 > 0.5)
  task = mlr3::TaskRegr$new("ranger_fast", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger", num.trees = 30L, mtry = 3L, min.node.size = 2L, num.threads = 1L)
  learner$train(task)
  x = data[, c("x1", "x2", "x3")]

  info = xplaineff:::extract_ranger_regression_model(learner)
  expect_s3_class(info$model, "ranger")
  expect_equal(info$feature_names, c("x1", "x2", "x3"))
  expected = as.numeric(learner$predict_newdata(x)$response)

  expect_equal(xplaineff:::predict_ranger_regression_fast(learner, x), expected)
  expect_equal(xplaineff:::default_predict_fun(learner, x), expected)
})

test_that("rpart regression fast prediction matches native predict", {
  skip_if_not_installed("rpart")
  set.seed(31L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = sin(data$x1) + data$x2 * (data$x3 > 0.5)
  model = rpart::rpart(
    y ~ .,
    data = data,
    control = rpart::rpart.control(minsplit = 5L, cp = 0.001)
  )
  x = data[, c("x1", "x2", "x3")]

  info = xplaineff:::extract_rpart_regression_model(model)
  expect_s3_class(info$model, "rpart")
  expected = as.numeric(predict(model, newdata = x, type = "vector"))

  expect_equal(xplaineff:::predict_rpart_regression_fast(model, x), expected)
  expect_equal(xplaineff:::default_predict_fun(model, x), expected)
})

test_that("mlr3 regr.rpart fast prediction bypasses the learner prediction object", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("rpart")
  set.seed(32L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = data$x1 - 0.5 * data$x2 + data$x3
  task = mlr3::TaskRegr$new("rpart_fast", backend = data, target = "y")
  learner = mlr3::lrn("regr.rpart", minsplit = 5L, cp = 0.001)
  learner$train(task)
  x = data[, c("x1", "x2", "x3")]

  info = xplaineff:::extract_rpart_regression_model(learner)
  expect_s3_class(info$model, "rpart")
  expect_equal(info$feature_names, c("x1", "x2", "x3"))
  expected = as.numeric(learner$predict_newdata(x)$response)

  expect_equal(xplaineff:::predict_rpart_regression_fast(learner, x), expected)
  expect_equal(xplaineff:::default_predict_fun(learner, x), expected)
})

test_that("native xgboost regression fast prediction matches booster predict", {
  skip_if_not_installed("xgboost")
  Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  set.seed(33L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = data$x1 - data$x2 + 0.25 * data$x3
  x = data[, c("x1", "x2", "x3")]
  x_mat = as.matrix(x)
  model = xgboost::xgboost(
    x = x_mat,
    y = data$y,
    nrounds = 8L,
    objective = "reg:squarederror",
    nthread = 1L,
    verbosity = 0L
  )

  info = xplaineff:::extract_xgboost_regression_model(model)
  expect_s3_class(info$model, "xgb.Booster")
  expected = as.numeric(predict(model, newdata = x_mat))

  expect_equal(xplaineff:::predict_xgboost_regression_fast(model, x), expected)
  expect_equal(xplaineff:::default_predict_fun(model, x), expected)
})

test_that("mlr3 regr.xgboost fast prediction uses the trained booster", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("xgboost")
  Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  set.seed(34L)
  data = data.frame(x1 = runif(80L), x2 = rnorm(80L), x3 = runif(80L))
  data$y = 0.5 * data$x1 - data$x2 + data$x3
  task = mlr3::TaskRegr$new("xgboost_fast", backend = data, target = "y")
  learner = mlr3::lrn("regr.xgboost", nrounds = 8L, nthread = 1L, verbosity = 0L)
  learner$train(task)
  x = data[, c("x1", "x2", "x3")]

  info = xplaineff:::extract_xgboost_regression_model(learner)
  expect_s3_class(info$model, "xgb.Booster")
  expect_equal(info$feature_names, c("x1", "x2", "x3"))
  expected = as.numeric(learner$predict_newdata(x)$response)

  expect_equal(xplaineff:::predict_xgboost_regression_fast(learner, x), expected)
  expect_equal(xplaineff:::default_predict_fun(learner, x), expected)
})

test_that("rpart fast prediction is used in the PD matrix path", {
  skip_if_not_installed("rpart")
  set.seed(35L)
  data = data.frame(x1 = runif(50L), x2 = rnorm(50L), x3 = runif(50L))
  data$y = sin(data$x1) + data$x2 * (data$x3 > 0.5)
  model = rpart::rpart(
    y ~ .,
    data = data,
    control = rpart::rpart.control(minsplit = 5L, cp = 0.001)
  )
  pred_fun = function(model, newdata) {
    as.numeric(predict(model, newdata = as.data.frame(newdata), type = "vector"))
  }

  fast = xplaineff:::calculate_pd_matrix(
    model = model,
    data = data,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    predict_fun = NULL,
    n_grid = 5L,
    pd_engine = "cpp"
  )
  expected = xplaineff:::calculate_pd_matrix(
    model = model,
    data = data,
    target_feature_name = "y",
    feature_set = c("x1", "x2"),
    predict_fun = pred_fun,
    n_grid = 5L,
    pd_engine = "cpp"
  )

  expect_equal(fast$Y$x1, expected$Y$x1)
  expect_equal(fast$Y$x2, expected$Y$x2)
})

test_that("mlr3 regr.xgboost fast prediction is used in the ALE path", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("xgboost")
  Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  set.seed(36L)
  data = data.frame(x1 = runif(50L), x2 = rnorm(50L), x3 = runif(50L))
  data$y = 0.5 * data$x1 - data$x2 + data$x3
  task = mlr3::TaskRegr$new("xgboost_fast_ale", backend = data, target = "y")
  learner = mlr3::lrn("regr.xgboost", nrounds = 8L, nthread = 1L, verbosity = 0L)
  learner$train(task)
  pred_fun = function(model, newdata) {
    as.numeric(model$predict_newdata(newdata)$response)
  }

  fast = xplaineff:::calculate_ale_fast(
    model = learner,
    data = data,
    feature_set = c("x1", "x2"),
    target_feature_name = "y",
    n_intervals = 5L,
    predict_fun = NULL
  )
  expected = xplaineff:::calculate_ale_fast(
    model = learner,
    data = data,
    feature_set = c("x1", "x2"),
    target_feature_name = "y",
    n_intervals = 5L,
    predict_fun = pred_fun
  )

  compare_cols = c("row_id", "feat_val", "x_left", "x_right", "d_l", "interval_index", "int_n", "int_s1", "int_s2")
  expect_equal(fast$x1[, compare_cols, with = FALSE], expected$x1[, compare_cols, with = FALSE])
  expect_equal(fast$x2[, compare_cols, with = FALSE], expected$x2[, compare_cols, with = FALSE])
})
