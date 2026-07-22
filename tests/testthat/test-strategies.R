test_that("PdStrategy can be created", {
  strategy = PdStrategy$new()
  expect_true(inherits(strategy, "PdStrategy"))
  expect_equal(strategy$name, "pd")
  expect_equal(strategy$categorical_split, "one_vs_rest")
  expect_equal(strategy$max_exhaustive_levels, 12L)
  exhaustive_strategy = PdStrategy$new(categorical_split = "exhaustive", max_exhaustive_levels = 8L)
  expect_equal(exhaustive_strategy$categorical_split, "exhaustive")
  expect_equal(exhaustive_strategy$max_exhaustive_levels, 8L)
})

test_that("PdStrategy find_best_split returns expected structure", {
  tryCatch({
    xplaineff:::search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
  set.seed(1)
  n = 30
  Z = data.frame(x = runif(n), y = runif(n))
  Y = list(
    matrix(rnorm(n * 2), ncol = 2),
    matrix(rnorm(n * 2), ncol = 2)
  )
  strategy = PdStrategy$new()
  res = strategy$find_best_split(Z = Z, Y = Y, min_node_size = 5, n_quantiles = NULL)
  expect_true(is.data.frame(res))
  expect_true(nrow(res) >= 1)
  expect_true(all(c("split_feature", "is_categorical", "split_point",
        "split_objective", "best_split") %in% names(res)))
})

test_that("fit-time active-effect hidden option keeps PD grids full", {
  tryCatch({
    xplaineff:::search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
  n = 10L
  signal = c(rep(2, 5L), rep(-2, 5L))
  weak = seq(-0.02, 0.02, length.out = n)
  Y = list(
    signal = cbind(left = signal, right = -signal),
    weak = cbind(left = weak, right = -weak)
  )
  effect = structure(
    list(Y = Y, grid = list(signal = c("left", "right"), weak = c("left", "right"))),
    class = "xplaineff_pd_matrix"
  )
  data = data.frame(signal = seq_len(n), weak = seq_len(n), y = 0)

  withr::local_options(list(xplaineff.active_effect_rel_tol = 1e-4))
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 0L, min_node_size = 2L)
  tree$fit(data = data, target_feature_name = "y", effect = effect)
  expect_named(tree$root$objective$value_j, "signal")
  expect_named(tree$root$grid, c("signal", "weak"))

  withr::local_options(list(xplaineff.active_effect_rel_tol = 0))
  unpruned_tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 0L, min_node_size = 2L)
  unpruned_tree$fit(data = data, target_feature_name = "y", effect = effect)
  expect_named(unpruned_tree$root$objective$value_j, c("signal", "weak"))
})

test_that("PdStrategy fit can use exhaustive categorical level-set splits", {
  tryCatch({
    xplaineff:::search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
  group = factor(rep(c("a", "b", "c", "d"), each = 2L), levels = c("a", "b", "c", "d"))
  curve_left = c(5, 5, -5, -5)
  curve_right = -curve_left
  Y = matrix(
    c(rep(curve_left, 4L), rep(curve_right, 4L)),
    nrow = 8L,
    byrow = TRUE,
    dimnames = list(NULL, levels(group))
  )
  effect = structure(
    list(Y = list(group = Y), grid = list(group = levels(group))),
    class = "xplaineff_pd_matrix"
  )
  data = data.frame(group = group, y = 0)

  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 1L, min_node_size = 2L, impr_par = 0)
  tree$fit(
    data = data,
    target_feature_name = "y",
    effect = effect,
    categorical_split = "exhaustive",
    max_exhaustive_levels = 4L
  )

  expect_equal(tree$strategy$categorical_split, "exhaustive")
  expect_equal(tree$strategy$max_exhaustive_levels, 4L)
  expect_equal(tree$root$split$levels, c("a", "b"))
  expect_equal(tree$root$children$left_child$parent$split_condition, "group in {a, b}")
  expect_equal(tree$root$children$right_child$parent$split_condition, "group in {c, d}")
})

test_that("PdStrategy heterogeneity returns numeric vector", {
  Y = list(matrix(rnorm(20), ncol = 2), matrix(rnorm(20), ncol = 2))
  strategy = PdStrategy$new()
  h = strategy$heterogeneity(Y)
  expect_true(is.numeric(h))
  expect_length(h, 2)
  expect_true(all(!is.na(h)))
  expect_true(all(h >= 0))
})

test_that("PdStrategy child objectives reuse non-split feature objectives only", {
  tryCatch({
    xplaineff:::search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
  strategy = PdStrategy$new()
  Y = list(
    x = matrix(
      c(
        1, 2, 3, 10, 11, 12,
        2, 3, 4, 11, 12, 13,
        3, 4, 5, 12, 13, 14,
        4, 5, 6, 13, 14, 15,
        5, 6, 7, 14, 15, 16,
        6, 7, 8, 15, 16, 17
      ),
      nrow = 6L,
      byrow = TRUE,
      dimnames = list(NULL, as.character(1:6))
    ),
    f = matrix(
      c(
        1, -1,
        1, -1,
        1, -1,
        -1, 1,
        -1, 1,
        -1, 1
      ),
      nrow = 6L,
      byrow = TRUE,
      dimnames = list(NULL, c("a", "b"))
    )
  )
  grid_left = list(x = as.character(1:3), f = c("a", "b"))
  grid_right = list(x = as.character(4:6), f = c("a", "b"))
  idx_left = 1:3
  idx_right = 4:6
  raw = xplaineff:::search_best_split_cpp(
    Z = data.frame(x = 1:6),
    Y = Y,
    min_node_size = 2L
  )
  split_info = list(
    split_feature = "x",
    split_value = 3.5,
    raw_result = raw
  )

  result = strategy$get_child_objectives(
    Z = NULL, Y = Y, split_info = split_info,
    idx_left = idx_left, idx_right = idx_right,
    grid_left = grid_left, grid_right = grid_right
  )
  y_left = strategy$node_transform(Y = Y, idx = idx_left, grid = grid_left)
  y_right = strategy$node_transform(Y = Y, idx = idx_right, grid = grid_right)
  expected_left = strategy$heterogeneity(y_left)
  expected_right = strategy$heterogeneity(y_right)

  expect_equal(result$left_objective_value_j, expected_left)
  expect_equal(result$right_objective_value_j, expected_right)
  expect_equal(result$left_objective_value, sum(expected_left))
  expect_equal(result$right_objective_value, sum(expected_right))
})

test_that("AleStrategy can be created", {
  strategy = AleStrategy$new()
  expect_true(inherits(strategy, "AleStrategy"))
  expect_equal(strategy$name, "ale")
  expect_equal(strategy$categorical_split, "ordered_prefix")
  expect_equal(strategy$max_exhaustive_levels, 12L)
  exhaustive_strategy = AleStrategy$new(categorical_split = "exhaustive", max_exhaustive_levels = 8L)
  expect_equal(exhaustive_strategy$categorical_split, "exhaustive")
  expect_equal(exhaustive_strategy$max_exhaustive_levels, 8L)
})

test_that("AleStrategy heterogeneity returns numeric for ALE-like list", {
  tryCatch({
    dt = data.table::data.table(row_id = 1:5, interval_index = rep(1L, 5), d_l = 0, int_n = 5L, int_s1 = 0, int_s2 = 0)
    xplaineff:::calculate_ale_heterogeneity_list_cpp(list(x = dt))
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("ALE C++ symbols not loaded (install package with compile)")
    }
  })
  n = 20
  dt = data.table::data.table(
    row_id = seq_len(n),
    interval_index = rep(1:4, length.out = n),
    d_l = rnorm(n),
    int_n = 5L, int_s1 = 0, int_s2 = 1
  )
  Y = list(f1 = dt)
  strategy = AleStrategy$new()
  h = strategy$heterogeneity(Y)
  expect_true(is.numeric(h))
  expect_length(h, 1)
  expect_true(!is.na(h))
  expect_true(h >= 0)
})

test_that("AleStrategy child objectives use the selected split when best splits tie", {
  strategy = AleStrategy$new()
  raw = data.frame(
    split_feature = rep(c("x", "z"), each = 2L),
    split_point = rep(c(0.5, 1.5), each = 2L),
    feature = rep(c("f1", "f2"), times = 2L),
    left_objective_value_j = c(1, 2, 100, 200),
    right_objective_value_j = c(3, 4, 300, 400),
    best_split = TRUE
  )
  split_info = list(
    split_feature = "x",
    split_value = 0.5,
    raw_result = raw
  )

  obj = strategy$get_child_objectives(
    Z = NULL, Y = NULL, split_info = split_info,
    idx_left = NULL, idx_right = NULL, grid_left = NULL, grid_right = NULL
  )

  expect_equal(obj$left_objective_value_j, c(f1 = 1, f2 = 2))
  expect_equal(obj$right_objective_value_j, c(f1 = 3, f2 = 4))
  expect_equal(obj$left_objective_value, 3)
  expect_equal(obj$right_objective_value, 7)
})
