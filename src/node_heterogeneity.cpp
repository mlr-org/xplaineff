/**
 * @file node_heterogeneity.cpp
 * @brief Heterogeneity (sum of variances) per feature from ICE matrices - C++ implementation.
 */

#include <Rcpp.h>

// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

// -----------------------------------------------------------------------------
// node_heterogeneity_cpp
// Purpose:
//   Compute heterogeneity (sum of variances) per feature from effect matrices.
//   For each matrix in Y: sum over columns of (sum of squares - sum^2 / n),
//   i.e., sum_j [ sum_i Y[i,j]^2 - (sum_i Y[i,j])^2 / n ].
//   Used by PdStrategy$heterogeneity().
// Inputs:
//   - Y: List of numeric matrices (e.g., ICE effect matrices per feature).
// Notes:
//   - NA/NaN values are excluded (na.rm = TRUE) per column.
//   - Empty or all-NA columns contribute 0.
// Output:
//   NumericVector of length length(Y).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector node_heterogeneity_cpp(List Y) {
  int n_mat = Y.size();
  NumericVector result(n_mat);

  for (int m = 0; m < n_mat; ++m) {
    NumericMatrix M = as<NumericMatrix>(Y[m]);
    int nrow = M.nrow();
    int ncol = M.ncol();
    double total = 0.0;

    for (int j = 0; j < ncol; ++j) {
      /* Per column: accumulate sum, sum of squares, and valid count (exclude NA/NaN). */
      double sum_val = 0.0;
      double sum_sq = 0.0;
      int n_valid = 0;

      for (int i = 0; i < nrow; ++i) {
        double v = M(i, j);
        if (!NumericMatrix::is_na(v)) {
          sum_val += v;
          sum_sq += v * v;
          n_valid++;
        }
      }

      /* Add n*Var = sum(x^2) - sum(x)^2/n for this column. */
      if (n_valid > 0) {
        total += sum_sq - (sum_val * sum_val) / static_cast<double>(n_valid);
      }
    }

    result[m] = total;
  }

  return result;
}
