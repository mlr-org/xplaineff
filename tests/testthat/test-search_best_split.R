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
