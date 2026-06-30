/**
 * @file ale_sweep.cpp
 * @brief Fast ALE main sweep for tree splitting - C++ implementation
 *
 * Implements the core loop of ALE-based split search: rows are ordered by the
 * split feature z; we sweep through split positions t=1..n_obs-1, moving one
 * row at a time from the right child to the left. At each position, we compute
 * the total heterogeneity (sum of SSE over intervals and features) and track
 * the best split. Uses sufficient statistics (n, s1, s2 per interval) for O(1)
 * incremental updates.
 */

#include <Rcpp.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

/* Expected raw self gain inside a crossed interval is one within-interval
 * variance unit under the null. The split score subtracts one estimated unit. */
static const double kSelfBiasLambda = 1.0;

// -----------------------------------------------------------------------------
// risk_from_stats
// Purpose:
//   Compute SSE from sufficient statistics: Var = E[X^2] - E[X]^2, so
//   SSE = sum(x^2) - sum(x)^2/n = s2 - s1^2/n. Used for per-interval risk.
// Notes:
//   Returns 0 when n <= 1 (no variance for single point).
// -----------------------------------------------------------------------------
inline double risk_from_stats(double n, double s1, double s2) {
  if (n <= 1.0) return 0.0;
  return s2 - (s1 * s1) / n;
}

inline double sample_var_from_stats(double n, double s1, double s2) {
  if (n <= 1.0) return 0.0;
  return risk_from_stats(n, s1, s2) / (n - 1.0);
}

// -----------------------------------------------------------------------------
// ale_sweep_cpp
// Purpose:
//   Find best split position t by sweeping rows (ordered by split feature) from
//   right to left. At each t, left child = rows 1..t, right = t+1..n_obs. Uses
//   sufficient statistics (n, s1, s2 per interval) for O(1) risk updates when
//   moving a row. For self-ALE splits, applies the appendix bias-corrected
//   objective to the split feature's own risk.
// Inputs:
//   ord_idx: 1-based row indices in sorted order (by z); length n_obs
//   d_l_mat: p x N, d_l_mat(j,i) = local effect for feature j, sample i
//   interval_idx_mat: p x N, interval index per feature per sample (1-based)
//   offsets: length p, start offset per feature in flattened tot_n/s1/s2 arrays
//   tot_n, tot_s1, tot_s2: full-node totals; r_n = right-node copy, l = tot - r
//   r_risks: initial right risks per feature (SSE per feature)
//   is_cand: length n_obs-1, TRUE where split position is a candidate (e.g. z changes)
//   min_node_size: minimum observations per child
//   split_feat_j: 1-based ALE feature index for split var, or 0 if no self-ALE
//   z_sorted: z values in ord order; for categorical, all 0
//   n_obs: number of observations in node
// Output:
//   List: best_t (1-based), best_risks_sum, best_left_risks, best_right_risks.
//   If no valid split: best_t = NA_INTEGER, best_risks_sum = Inf.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List ale_sweep_cpp(
    IntegerVector ord_idx,
    NumericMatrix d_l_mat,
    IntegerMatrix interval_idx_mat,
    IntegerVector offsets,
    NumericVector tot_n,
    NumericVector tot_s1,
    NumericVector tot_s2,
    NumericVector r_risks,
    LogicalVector is_cand,
    int min_node_size,
    int split_feat_j,
    NumericVector z_sorted,
    int n_obs
) {
  const int p = d_l_mat.nrow();
  const int M = tot_n.size();
  const int j0 = split_feat_j - 1;  /* 0-based index of split feature. */

  /* Working copies: r_n, r_s1, r_s2 = right-node sufficient stats per interval. */
  std::vector<double> r_n(tot_n.begin(), tot_n.end());
  std::vector<double> r_s1(tot_s1.begin(), tot_s1.end());
  std::vector<double> r_s2(tot_s2.begin(), tot_s2.end());
  std::vector<double> left_risks(p, 0.0);
  std::vector<double> right_risks(r_risks.begin(), r_risks.end());

  double risks_sum = 0.0;
  for (int j = 0; j < p; ++j) risks_sum += right_risks[j];

  double best_risks_sum = R_PosInf;
  int best_t = -1;
  std::vector<double> best_left_risks(p), best_right_risks(p);

  bool has_self_ale = (split_feat_j >= 1 && split_feat_j <= p);

  const double self_root_risk = has_self_ale ? r_risks[j0] : 0.0;

  /* Precompute interval indices in sweep order for self-bias correction. */
  IntegerVector interval_idx_sorted(n_obs);
  const int N = d_l_mat.ncol();
  if (has_self_ale && z_sorted.size() >= (R_xlen_t)n_obs) {
    for (int i = 0; i < n_obs; ++i) {
      int row = ord_idx[i] - 1;
      if (row < 0 || row >= N) stop("ord_idx contains invalid row index");
      interval_idx_sorted[i] = interval_idx_mat(j0, row);
    }
  }
  /* Sweep: at each t, move row ord_idx[t-1] from right to left. */
  for (int t = 1; t <= n_obs - 1; ++t) {
    int row = ord_idx[t - 1] - 1;  /* 1-based to 0-based. */
    if (row < 0 || row >= N) stop("ord_idx contains invalid row index");

    /* Update sufficient stats and risks for each feature. */
    for (int j = 0; j < p; ++j) {
      double d = d_l_mat(j, row);
      int interval_idx_val = interval_idx_mat(j, row);
      int m = offsets[j] + interval_idx_val - 1;  /* Flattened interval index. */
      if (m < 0 || m >= M) continue;

      /* Remove row from right: r_n -= 1, r_s1 -= d, r_s2 -= d^2. */
      double r_n_old = r_n[m];
      double r_s1_old = r_s1[m];
      double r_s2_old = r_s2[m];
      double r_risk_old = risk_from_stats(r_n_old, r_s1_old, r_s2_old);

      double r_n_new = r_n_old - 1.0;
      double r_s1_new = r_s1_old - d;
      double r_s2_new = r_s2_old - d * d;
      double r_risk_new = risk_from_stats(r_n_new, r_s1_new, r_s2_new);

      r_n[m] = r_n_new;
      r_s1[m] = r_s1_new;
      r_s2[m] = r_s2_new;
      right_risks[j] -= r_risk_old;
      right_risks[j] += r_risk_new;

      /* Left stats = tot - right (invariant: tot = left + right). */
      double l_n_old = tot_n[m] - r_n_old;
      double l_s1_old = tot_s1[m] - r_s1_old;
      double l_s2_old = tot_s2[m] - r_s2_old;
      double l_risk_old = risk_from_stats(l_n_old, l_s1_old, l_s2_old);

      double l_n_new = tot_n[m] - r_n_new;
      double l_s1_new = tot_s1[m] - r_s1_new;
      double l_s2_new = tot_s2[m] - r_s2_new;
      double l_risk_new = risk_from_stats(l_n_new, l_s1_new, l_s2_new);

      left_risks[j] -= l_risk_old;
      left_risks[j] += l_risk_new;
      risks_sum += (-l_risk_old - r_risk_old + l_risk_new + r_risk_new);
    }

    /* Skip if not a candidate position or violates min_node_size. */
    if (!is_cand[t - 1] || t < min_node_size || (n_obs - t) < min_node_size)
      continue;

    /* Left/right constant? (all same z value) -> drop self risk. */
    bool l_const = (z_sorted.size() > 0 && std::abs(z_sorted[0] - z_sorted[t - 1]) < 1e-15);
    bool r_const = (z_sorted.size() >= (R_xlen_t)n_obs && std::abs(z_sorted[t] - z_sorted[n_obs - 1]) < 1e-15);

    /* Self-ALE splits are ranked by other-feature child risk minus the
     * bias-corrected self gain. */
    double total = R_PosInf;

    if (!has_self_ale) {
      total = risks_sum;
    } else {
      double self_risk = left_risks[j0] + right_risks[j0];
      double drop = 0.0;
      if (l_const) drop += left_risks[j0];
      if (r_const) drop += right_risks[j0];
      double other_risks = risks_sum - self_risk;
      double self_risk_effective = self_risk - drop;
      double delta_raw = self_root_risk - self_risk_effective;
      double self_bias = 0.0;
      bool cut_self_interval =
        !l_const && !r_const &&
        z_sorted.size() >= (R_xlen_t)n_obs &&
        interval_idx_sorted[t - 1] == interval_idx_sorted[t];
      if (cut_self_interval) {
        int m_self = offsets[j0] + interval_idx_sorted[t - 1] - 1;
        if (m_self >= 0 && m_self < M) {
          self_bias = kSelfBiasLambda * sample_var_from_stats(
            tot_n[m_self], tot_s1[m_self], tot_s2[m_self]
          );
        }
      }
      double delta_corr = std::max(0.0, delta_raw - self_bias);
      total = other_risks - delta_corr;
    }

    if (!R_FINITE(total) || total >= best_risks_sum) continue;

    /* New best: store split index and risk vectors. */
    best_risks_sum = total;
    best_t = t;
    best_left_risks = left_risks;
    best_right_risks = right_risks;

    /* For self-ALE feature: store adjusted risks (or 0 if dropped). */
    if (has_self_ale) {
      if (l_const) best_left_risks[j0] = 0.0;
      if (r_const) best_right_risks[j0] = 0.0;
    }
  }

  /* No valid split found. */
  if (best_t < 0) {
    return List::create(
      _["best_t"] = NA_INTEGER,
      _["best_risks_sum"] = R_PosInf,
      _["best_left_risks"] = NumericVector(p, NA_REAL),
      _["best_right_risks"] = NumericVector(p, NA_REAL)
    );
  }

  return List::create(
    _["best_t"] = best_t,  /* 1-based: split after row best_t. */
    _["best_risks_sum"] = best_risks_sum,
    _["best_left_risks"] = NumericVector(best_left_risks.begin(), best_left_risks.end()),
    _["best_right_risks"] = NumericVector(best_right_risks.begin(), best_right_risks.end())
  );
}
