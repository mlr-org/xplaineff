skip_cpp_if_unavailable = function() {
  tryCatch({
    search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
}

test_that("search_best_split_cpp works with numeric data", {
  skip_cpp_if_unavailable()
  set.seed(1)
  n = 25
  Z = data.frame(x = runif(n), y = runif(n))
  Y = list(matrix(rnorm(n * 2), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 5)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 1)
  expect_true(all(c("split_feature", "is_categorical", "split_point",
        "split_objective", "best_split") %in% names(result)))
  expect_true(all(result$is_categorical %in% c(TRUE, FALSE)))
})

test_that("search_best_split_cpp works with single feature", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:10)
  Y = list(matrix(rnorm(20), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$split_feature[1], "x")
})

test_that("search_best_split_cpp respects min_node_size", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:10)
  Y = list(matrix(rnorm(20), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 10)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_false(result$best_split[1])
})

test_that("search_best_split_cpp with multiple Y matrices", {
  skip_cpp_if_unavailable()
  n = 20
  Z = data.frame(x = runif(n))
  Y = list(
    matrix(rnorm(n * 2), ncol = 2),
    matrix(rnorm(n * 2), ncol = 2)
  )
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("search_best_split_cpp returns child objective list columns", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:6)
  Y = list(
    f1 = matrix(c(rep(2, 3L), rep(-2, 3L)), ncol = 1L),
    f2 = matrix(0, nrow = 6L, ncol = 1L)
  )

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2L)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1L)
  expect_true(is.list(result$left_objective_value_j))
  expect_true(is.list(result$right_objective_value_j))
  expect_named(result$left_objective_value_j[[1L]], c("f1", "f2"))
  expect_named(result$right_objective_value_j[[1L]], c("f1", "f2"))
  expect_equal(result$left_objective_value_j[[1L]], c(f1 = 0, f2 = 0))
  expect_equal(result$right_objective_value_j[[1L]], c(f1 = 0, f2 = 0))
})

test_that("search_best_split_cpp prunes numerically inactive effect matrices", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:10)
  Y = list(
    signal = matrix(c(rep(2, 5L), rep(-2, 5L)), ncol = 1L),
    inactive = matrix(seq_len(10L) * 1e-14, ncol = 1L)
  )

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2L)

  best = which(result$best_split)[1L]
  expect_named(result$left_objective_value_j[[best]], c("signal", "inactive"))
  expect_named(result$right_objective_value_j[[best]], c("signal", "inactive"))
  expect_equal(unname(result$left_objective_value_j[[best]]["inactive"]), 0)
  expect_equal(unname(result$right_objective_value_j[[best]]["inactive"]), 0)
})

test_that("search_best_split_cpp handles categorical split and one-level factor", {
  skip_cpp_if_unavailable()
  Z = data.frame(
    group = factor(c("a", "a", "a", "b", "b", "b")),
    one = factor(rep("z", 6L))
  )
  Y = list(matrix(c(5, 5, 5, -5, -5, -5), ncol = 1L))

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2L)

  expect_true(result$is_categorical[1L])
  expect_equal(result$split_point[1L], "a")
  expect_true(result$best_split[1L])
  expect_equal(result$split_objective[2L], Inf)
  expect_false(result$best_split[2L])
})

test_that("search_best_split_cpp handles categorical self effect grids", {
  skip_cpp_if_unavailable()
  Z = data.frame(cat = factor(c("a", "a", "b", "b", "c", "c")))
  Y = list(cat = matrix(
    c(
      1, 2, 3,
      2, 3, 4,
      4, 3, 2,
      5, 4, 3,
      3, 5, 7,
      4, 6, 8
    ),
    nrow = 6L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  ))

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 1L)

  expect_true(is.data.frame(result))
  expect_true(result$is_categorical[1L])
  expect_true(is.finite(result$split_objective[1L]))
})

test_that("search_best_split_cpp halves numeric self-effect grids", {
  skip_cpp_if_unavailable()
  z = 1:6
  Y = matrix(
    c(
      1, 2, 3, 9, 10, 11,
      1, 2, 3, 9, 10, 11,
      1, 2, 3, 9, 10, 11,
      10, 11, 12, 1, 2, 3,
      10, 11, 12, 1, 2, 3,
      10, 11, 12, 1, 2, 3
    ),
    nrow = 6L,
    byrow = TRUE,
    dimnames = list(NULL, as.character(z))
  )
  objective = function(split_value) {
    grid = as.numeric(colnames(Y))
    idx_left = z <= split_value
    grid_left = grid <= split_value
    left = Y[idx_left, grid_left, drop = FALSE]
    right = Y[!idx_left, !grid_left, drop = FALSE]
    parent_ss = sum(colSums(Y^2))
    sum(colSums(left^2) - colSums(left)^2 / sum(idx_left)) +
      sum(colSums(right^2) - colSums(right)^2 / sum(!idx_left)) -
      parent_ss
  }
  expected = vapply(z[-length(z)], objective, numeric(1L))

  result = search_best_split_cpp(Z = data.frame(x = z), Y = list(x = Y), min_node_size = 2L)

  expect_equal(as.numeric(result$split_point[1L]), 3.5)
  expect_equal(result$split_objective[1L], min(expected), tolerance = 1e-10)
})

test_that("search_best_split_cpp handles character numeric features and quantile candidates", {
  skip_cpp_if_unavailable()
  Z = data.frame(
    x = as.character(1:8),
    stringsAsFactors = FALSE
  )
  Y = list(matrix(c(rep(2, 4L), rep(-2, 4L)), ncol = 1L))

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2L, n_quantiles = 3L)

  expect_false(result$is_categorical[1L])
  expect_equal(as.numeric(result$split_point[1L]), 4.5)
  expect_true(result$best_split[1L])
})

test_that("search_best_split_cpp treats NaN effect values as zero", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:6)
  Y = list(matrix(c(1, NaN, 1, -1, NaN, -1), ncol = 1L))

  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2L)

  expect_true(is.finite(result$split_objective[1L]))
  expect_true(result$best_split[1L])
})
