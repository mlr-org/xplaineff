test_that("GadgetTree can be created", {
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 4)

  expect_true(inherits(tree, "GadgetTree"))
  expect_equal(tree$n_split, 4)
  expect_true(inherits(tree$strategy, "PdStrategy"))
})

test_that("GadgetTree fit with ALE strategy works", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(123)
  n = 80
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 2, min_node_size = 20)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  expect_true(!is.null(tree$root))
  expect_true(inherits(tree$root, "R6"))
})

test_that("GadgetTree plot_tree_structure works after fit", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(456)
  n = 60
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1, min_node_size = 15)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  p = tree$plot_tree_structure()
  expect_true(inherits(p, "gg"))
})

test_that("GadgetTree extract_split_info returns data frame", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(789)
  n = 60
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 2, min_node_size = 15)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  split_info = tree$extract_split_info()
  expect_true(is.data.frame(split_info))
  expect_true(nrow(split_info) >= 1)
  expect_true("depth" %in% names(split_info))
  expect_true("id" %in% names(split_info))
  expect_identical(names(split_info)[1:2], c("id", "depth"))
  expect_false("time" %in% names(split_info))

  timed_split_info = tree$extract_split_info(include_timing = TRUE)
  expect_true(is.data.frame(timed_split_info))
  if (length(tree$split_benchmark) > 0L) {
    expect_true("time" %in% names(timed_split_info))
  }
})

test_that("GadgetTree plot (ALE) returns list of plots", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(101)
  n = 50
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1, min_node_size = 15)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  plot_result = tree$plot(data = data, target_feature_name = "y", show_plot = FALSE)
  expect_true(is.list(plot_result))
  # Structure: plot_result[[depth_name]] = list(Node_id = patchwork, ...); get first actual plot
  if (length(plot_result) > 0 && length(plot_result[[1]]) > 0) {
    first_plot = plot_result[[1]][[1]]
    expect_true(inherits(first_plot, "gg") || inherits(first_plot, "patchwork"))
  }
})
