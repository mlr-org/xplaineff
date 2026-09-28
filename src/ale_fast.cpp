#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

// Binary search helper to find interval index for x in left-open, right-closed
// intervals (q_k, q_{k+1}]. Intervals are defined by cutpoints in q (assumed to
// be sorted ascending) and the return value is in {0, ..., qlen-2}.
//
// Boundary handling:
// - x <= q[0] is mapped to interval 0.
// - x >= q[qlen-1] is mapped to the last interval (qlen-2).
static int find_interval(double x, const double *q, R_xlen_t qlen) {
  if (x <= q[0]) {
    return 0;
  }
  if (x >= q[qlen - 1]) {
    return (int)(qlen - 2);
  }

  R_xlen_t lo = 0;
  R_xlen_t hi = qlen - 1;
  // Invariant: q[lo] < x <= q[hi]
  while (hi - lo > 1) {
    const R_xlen_t mid = lo + (hi - lo) / 2;
    if (q[mid] < x) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (int)lo;
}

static double quantile_type7_sorted(const NumericVector& x_sorted, double p) {
  const int n = x_sorted.size();
  if (n <= 0) return NA_REAL;
  if (n == 1) return x_sorted[0];
  if (p <= 0.0) return x_sorted[0];
  if (p >= 1.0) return x_sorted[n - 1];
  const double h = (n - 1) * p;
  const int lo = static_cast<int>(h);
  const double f = h - lo;
  return x_sorted[lo] + f * (x_sorted[lo + 1] - x_sorted[lo]);
}

// -----------------------------------------------------------------------------
// cpp_ale_numeric_breaks
// Purpose:
//   Compute ALE interval breaks for a numeric feature as type-7 quantiles at
//   probabilities 0, 1/n_intervals, ..., 1 of the finite values of x, with
//   duplicate breaks collapsed.
// Inputs:
//   x: numeric feature values (non-finite values are ignored)
//   n_intervals: requested number of intervals (>= 1)
// Output:
//   Numeric vector of unique, ascending breaks with at least two entries
//   (c(NA, NA) if x has no finite value; c(q, q) if all finite values equal).
// -----------------------------------------------------------------------------
static NumericVector cpp_ale_numeric_breaks(NumericVector x, int n_intervals) {
  if (n_intervals < 1) {
    stop("n_intervals must be >= 1.");
  }
  std::vector<double> x_finite;
  x_finite.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      x_finite.push_back(x[i]);
    }
  }
  if (x_finite.size() == 0) {
    return NumericVector::create(NA_REAL, NA_REAL);
  }
  std::sort(x_finite.begin(), x_finite.end());
  NumericVector x_sorted(x_finite.begin(), x_finite.end());
  NumericVector q(n_intervals + 1);
  for (int i = 0; i <= n_intervals; ++i) {
    const double p = static_cast<double>(i) / static_cast<double>(n_intervals);
    q[i] = quantile_type7_sorted(x_sorted, p);
  }
  q.erase(std::unique(q.begin(), q.end()), q.end());
  if (q.size() < 2) {
    q = NumericVector::create(q[0], q[0]);
  }
  return q;
}

// Largest positive (non-NA) value of a 1-based interval index vector; 0 if none.
static int max_valid_index(const IntegerVector& interval_index) {
  int k = 0;
  for (int i = 0; i < interval_index.size(); ++i) {
    if (interval_index[i] != NA_INTEGER && interval_index[i] > k) {
      k = interval_index[i];
    }
  }
  return k;
}

// Accumulate per-interval sufficient statistics (count, sum, sum of squares)
// of d over observations with a valid 1-based index in 1..k and finite d.
static void accumulate_interval_stats(
  const IntegerVector& interval_index,
  const NumericVector& d,
  int k,
  NumericVector& interval_n,
  NumericVector& interval_s1,
  NumericVector& interval_s2
) {
  const int n = d.size();
  for (int i = 0; i < n; ++i) {
    const int idx = interval_index[i];
    const double di = d[i];
    if (idx == NA_INTEGER || idx < 1 || idx > k || !R_finite(di)) {
      continue;
    }
    const int pos = idx - 1;
    interval_n[pos] += 1.0;
    interval_s1[pos] += di;
    interval_s2[pos] += di * di;
  }
}

// -----------------------------------------------------------------------------
// cpp_ale_interval_aggregate
// Purpose:
//   Aggregate local effects into per-interval sufficient statistics.
// Inputs:
//   d_l: local effects per observation (non-finite values are skipped)
//   interval_index: 1-based interval id per observation (NA skipped)
// Output:
//   List of length-k vectors (k = max interval id): interval_n (counts),
//   interval_s1 (sum of d_l), interval_s2 (sum of d_l^2). Empty vectors if no
//   valid index.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_ale_interval_aggregate(NumericVector d_l, IntegerVector interval_index) {
  int n = d_l.size();
  if (interval_index.size() != n) {
    stop("d_l and interval_index must have the same length.");
  }

  const int k = max_valid_index(interval_index);
  if (k <= 0) {
    return List::create(
      _["interval_n"] = NumericVector(0),
      _["interval_s1"] = NumericVector(0),
      _["interval_s2"] = NumericVector(0)
    );
  }

  NumericVector interval_n(k, 0.0);
  NumericVector interval_s1(k, 0.0);
  NumericVector interval_s2(k, 0.0);
  accumulate_interval_stats(interval_index, d_l, k, interval_n, interval_s1, interval_s2);

  return List::create(
    _["interval_n"] = interval_n,
    _["interval_s1"] = interval_s1,
    _["interval_s2"] = interval_s2
  );
}

// -----------------------------------------------------------------------------
// ale_finalize_local_effects
// Purpose:
//   Shared body of the numeric and categorical ALE finalization: compute local
//   effects d_l = upper - lower, accumulate per-interval sufficient statistics,
//   and broadcast the statistics of each observation's interval back to rows.
// Inputs:
//   interval_index: 1-based interval (or level) id per observation
//   upper, lower: predictions at the upper/lower interval boundary per observation
//   index_name: name of the index argument, used in the error message
// Output:
//   List of length-n vectors: d_l (NA where the index is invalid or a prediction
//   is non-finite), int_n, int_s1, int_s2 (statistics of the row's interval).
// -----------------------------------------------------------------------------
static List ale_finalize_local_effects(
  const IntegerVector& interval_index,
  const NumericVector& upper,
  const NumericVector& lower,
  const char* index_name
) {
  const int n = interval_index.size();
  const int k = max_valid_index(interval_index);
  if (k <= 0) {
    stop("%s must contain at least one valid positive index.", index_name);
  }

  NumericVector d_l(n, NA_REAL);
  for (int i = 0; i < n; ++i) {
    const int idx = interval_index[i];
    if (idx == NA_INTEGER || idx < 1 || idx > k) {
      continue;
    }
    const double lo = lower[i];
    const double up = upper[i];
    if (!R_finite(lo) || !R_finite(up)) {
      continue;
    }
    d_l[i] = up - lo;
  }

  NumericVector interval_n(k, 0.0);
  NumericVector interval_s1(k, 0.0);
  NumericVector interval_s2(k, 0.0);
  accumulate_interval_stats(interval_index, d_l, k, interval_n, interval_s1, interval_s2);

  NumericVector int_n_row(n, NA_REAL);
  NumericVector int_s1_row(n, NA_REAL);
  NumericVector int_s2_row(n, NA_REAL);
  for (int i = 0; i < n; ++i) {
    const int idx = interval_index[i];
    if (idx == NA_INTEGER || idx < 1 || idx > k) {
      continue;
    }
    const int pos = idx - 1;
    int_n_row[i] = interval_n[pos];
    int_s1_row[i] = interval_s1[pos];
    int_s2_row[i] = interval_s2[pos];
  }

  return List::create(
    _["d_l"] = d_l,
    _["int_n"] = int_n_row,
    _["int_s1"] = int_s1_row,
    _["int_s2"] = int_s2_row
  );
}

// -----------------------------------------------------------------------------
// cpp_ale_numeric_effect_table
// Purpose:
//   Build the per-observation ALE effect table for a numeric feature.
// Inputs:
//   feat_val: feature value per observation
//   x_left, x_right: lower/upper break of the observation's interval
//   interval_index: 1-based interval id per observation
//   preds_lower, preds_upper: predictions with the feature set to x_left/x_right
// Output:
//   List of length-n columns: row_id, feat_val, x_left, x_right,
//   d_l (= preds_upper - preds_lower), interval_index, int_n, int_s1, int_s2.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_ale_numeric_effect_table(
  NumericVector feat_val,
  NumericVector x_left,
  NumericVector x_right,
  IntegerVector interval_index,
  NumericVector preds_lower,
  NumericVector preds_upper
) {
  const int n = feat_val.size();
  if (x_left.size() != n || x_right.size() != n || interval_index.size() != n ||
      preds_lower.size() != n || preds_upper.size() != n) {
    stop("All inputs must have the same length.");
  }

  List finalized = ale_finalize_local_effects(interval_index, preds_upper, preds_lower, "interval_index");
  IntegerVector row_id(n);
  for (int i = 0; i < n; ++i) {
    row_id[i] = i + 1;
  }

  return List::create(
    _["row_id"] = row_id,
    _["feat_val"] = feat_val,
    _["x_left"] = x_left,
    _["x_right"] = x_right,
    _["d_l"] = finalized["d_l"],
    _["interval_index"] = interval_index,
    _["int_n"] = finalized["int_n"],
    _["int_s1"] = finalized["int_s1"],
    _["int_s2"] = finalized["int_s2"]
  );
}

// -----------------------------------------------------------------------------
// cpp_ale_categorical_effect_table
// Purpose:
//   Build the per-observation ALE effect table for a categorical feature whose
//   levels have been ordered.
// Inputs:
//   feat_val: observed level id per observation
//   x_left, x_right: neighbouring level ids used for the finite difference
//   interval_index: 1-based level id per observation
//   y_hat_plus, y_hat_neg: predictions with the feature set to x_right/x_left
// Output:
//   List of length-n columns: row_id, feat_val, x_left, x_right,
//   d_l (= y_hat_plus - y_hat_neg), interval_index, int_n, int_s1, int_s2.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_ale_categorical_effect_table(
  IntegerVector feat_val,
  IntegerVector x_left,
  IntegerVector x_right,
  IntegerVector interval_index,
  NumericVector y_hat_plus,
  NumericVector y_hat_neg
) {
  const int n = feat_val.size();
  if (x_left.size() != n || x_right.size() != n || interval_index.size() != n ||
      y_hat_plus.size() != n || y_hat_neg.size() != n) {
    stop("All inputs must have the same length.");
  }

  List finalized = ale_finalize_local_effects(interval_index, y_hat_plus, y_hat_neg, "levels_id");
  IntegerVector row_id(n);
  for (int i = 0; i < n; ++i) {
    row_id[i] = i + 1;
  }

  return List::create(
    _["row_id"] = row_id,
    _["feat_val"] = feat_val,
    _["x_left"] = x_left,
    _["x_right"] = x_right,
    _["d_l"] = finalized["d_l"],
    _["interval_index"] = interval_index,
    _["int_n"] = finalized["int_n"],
    _["int_s1"] = finalized["int_s1"],
    _["int_s2"] = finalized["int_s2"]
  );
}

// -----------------------------------------------------------------------------
// cpp_ale_numeric_prepare
// Purpose:
//   Prepare the ALE interval layout for a numeric feature: compute quantile
//   breaks and assign each observation to its interval.
// Inputs:
//   x: numeric feature values (non-finite values get NA index/bounds)
//   n_intervals: requested number of intervals (>= 1)
// Output:
//   List: zero_effect (TRUE if x is empty or all breaks coincide, in which case
//   interval_index/x_left/x_right are all NA), breaks, interval_index (1-based),
//   x_left, x_right (lower/upper break of each observation's interval).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_ale_numeric_prepare(NumericVector x, int n_intervals) {
  if (n_intervals < 1) {
    stop("n_intervals must be >= 1.");
  }
  const int n = x.size();
  NumericVector x_left(n, NA_REAL);
  NumericVector x_right(n, NA_REAL);
  IntegerVector interval_index(n, NA_INTEGER);

  if (n <= 0) {
    return List::create(
      _["zero_effect"] = true,
      _["breaks"] = NumericVector::create(NA_REAL, NA_REAL),
      _["interval_index"] = interval_index,
      _["x_left"] = x_left,
      _["x_right"] = x_right
    );
  }

  NumericVector breaks = cpp_ale_numeric_breaks(x, n_intervals);
  const int nb = breaks.size();
  if (nb < 2 || !R_finite(breaks[0]) || !R_finite(breaks[nb - 1]) || breaks[0] == breaks[nb - 1]) {
    return List::create(
      _["zero_effect"] = true,
      _["breaks"] = breaks,
      _["interval_index"] = interval_index,
      _["x_left"] = x_left,
      _["x_right"] = x_right
    );
  }

  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    if (!R_finite(xi)) {
      continue;
    }
    const int idx = find_interval(xi, REAL(breaks), nb) + 1; // 1-based
    interval_index[i] = idx;
    x_left[i] = breaks[idx - 1];
    x_right[i] = breaks[idx];
  }

  return List::create(
    _["zero_effect"] = false,
    _["breaks"] = breaks,
    _["interval_index"] = interval_index,
    _["x_left"] = x_left,
    _["x_right"] = x_right
  );
}

// -----------------------------------------------------------------------------
// cpp_ale_categorical_prepare
// Purpose:
//   For an ordered categorical feature, find the neighbouring level ids used for
//   the ALE finite difference: left = lv - 1 (or lv at the first level) and
//   right = lv + 1 (or lv at the last level).
// Inputs:
//   levels_id: 1-based ordered level id per observation
//   n_levels: number of levels (>= 1); ids outside 1..n_levels get NA
// Output:
//   List: left_id, right_id (integer vectors of length n).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_ale_categorical_prepare(IntegerVector levels_id, int n_levels) {
  if (n_levels < 1) {
    stop("n_levels must be >= 1.");
  }
  const int n = levels_id.size();
  IntegerVector left_id(n, NA_INTEGER);
  IntegerVector right_id(n, NA_INTEGER);

  for (int i = 0; i < n; ++i) {
    const int lv = levels_id[i];
    if (lv == NA_INTEGER || lv < 1 || lv > n_levels) {
      continue;
    }
    left_id[i] = lv;
    right_id[i] = lv;
    if (lv < n_levels) {
      right_id[i] = lv + 1;
    }
    if (lv > 1) {
      left_id[i] = lv - 1;
    }
  }

  return List::create(
    _["left_id"] = left_id,
    _["right_id"] = right_id
  );
}
