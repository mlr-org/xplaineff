/**
 * @file search_best_split.cpp
 * @brief Fast tree splitting for PD strategy - C++/Armadillo implementation.
 *
 * Categorical and numerical splitting with preprocessed effect matrices.
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Helper functions (internal)
// -----------------------------------------------------------------------------

/* Remove consecutive duplicates; assumes x is sorted. O(n). Copies input so caller's x is unchanged. */
inline NumericVector unique_cpp(const NumericVector& x) {
  NumericVector out = Rcpp::clone(x);
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

/* R type-7 quantile for sorted vector (R default). */
inline double quantile_type7(const NumericVector& x, double p) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];
  if (p <= 0.0) return x[0];
  if (p >= 1.0) return x[n - 1];

  const double h = (n - 1) * p;
  const int lo   = static_cast<int>(h);
  const double f = h - lo;

  return x[lo] + f * (x[lo + 1] - x[lo]);
}

/* Convert R matrix to Armadillo (zero-copy when possible). */
inline arma::mat arma_view(SEXP obj) {
  if (!Rf_isMatrix(obj) || TYPEOF(obj) != REALSXP) {
    Rcpp::Function as_matrix("as.matrix");
    obj = as_matrix(obj);
  }
  Rcpp::NumericMatrix M(obj);
  return arma::mat(M.begin(), M.nrow(), M.ncol(), false);
}

/* Count of grid points <= x via binary search (O(log n)). Precondition: grid sorted ascending
   (checked once per feature by the caller, not here, since this runs per split candidate).
   Mirrors R's findInterval default. */
inline int find_grid_interval(double x, const std::vector<double>& grid) {
  return static_cast<int>(std::upper_bound(grid.begin(), grid.end(), x) - grid.begin());
}

inline bool parse_double_strict(const std::string& s, double& out) {
  if (s.empty()) return false;
  const char* start = s.c_str();
  char* end = NULL;
  errno = 0;
  out = std::strtod(start, &end);
  if (start == end || errno == ERANGE || !R_finite(out)) return false;
  while (*end != '\0') {
    if (!std::isspace(static_cast<unsigned char>(*end))) return false;
    ++end;
  }
  return true;
}

/* Extract the numeric grid (column names) of an R matrix. Empty if names are absent. */
inline std::vector<double> grid_from_matrix_colnames(SEXP m) {
  std::vector<double> g;
  SEXP dn = Rf_getAttrib(m, R_DimNamesSymbol);
  if (Rf_isNull(dn)) return g;
  SEXP cn = VECTOR_ELT(dn, 1);
  if (Rf_isNull(cn)) return g;
  CharacterVector colnames(cn);
  g.reserve(colnames.size());
  for (R_xlen_t k = 0; k < colnames.size(); ++k) {
    if (colnames[k] == NA_STRING) return std::vector<double>();
    double value = NA_REAL;
    if (!parse_double_strict(as<std::string>(colnames[k]), value)) return std::vector<double>();
    g.push_back(value);
  }
  return g;
}

List child_objectives_from_flat_sums(
    const arma::vec& SL,
    const arma::vec& QL,
    const arma::vec& S_tot,
    const arma::vec& Q_tot,
    const std::vector<int>& offsets,
    int NL,
    int NR,
    int split_feat_effect = -1,      // effect index of the split feature's own effect (-1 = none)
    int split_feat_pos = 0           // # of that feature's grid points <= the split value
) {
  const int Ly = offsets.size() - 1;
  NumericVector left_obj(Ly, NA_REAL);
  NumericVector right_obj(Ly, NA_REAL);
  if (NL <= 0 || NR <= 0) {
    return List::create(
      _["left_objective_value_j"] = left_obj,
      _["right_objective_value_j"] = right_obj
    );
  }
  for (int l = 0; l < Ly; ++l) {
    if (offsets[l] == offsets[l + 1]) {
      left_obj[l] = 0.0;
      right_obj[l] = 0.0;
      continue;
    }
    const int start = offsets[l];
    const int end = offsets[l + 1] - 1;
    // Split feature: its ICE grid is divided across the children, so the left child's
    // objective covers only grid columns <= split (positions [start, mid)) and the right
    // child's only columns > split ([mid, end]). See search_best_split_point_cpp_internal.
    if (l == split_feat_effect) {
      const int mid = start + split_feat_pos;   // first right-surviving column
      if (mid > start) {
        const arma::vec SL_l = SL.subvec(start, mid - 1);
        const arma::vec QL_l = QL.subvec(start, mid - 1);
        left_obj[l] = arma::accu(QL_l - SL_l%SL_l / NL);
      } else {
        left_obj[l] = 0.0;
      }
      if (mid <= end) {
        const arma::vec SR_l = S_tot.subvec(mid, end) - SL.subvec(mid, end);
        const arma::vec QR_l = Q_tot.subvec(mid, end) - QL.subvec(mid, end);
        right_obj[l] = arma::accu(QR_l - SR_l%SR_l / NR);
      } else {
        right_obj[l] = 0.0;
      }
      continue;
    }
    const arma::vec SL_l = SL.subvec(start, end);
    const arma::vec QL_l = QL.subvec(start, end);
    const arma::vec SR_l = S_tot.subvec(start, end) - SL_l;
    const arma::vec QR_l = Q_tot.subvec(start, end) - QL_l;
    left_obj[l] = arma::accu(QL_l - SL_l%SL_l / NL);
    right_obj[l] = arma::accu(QR_l - SR_l%SR_l / NR);
  }
  return List::create(
    _["left_objective_value_j"] = left_obj,
    _["right_objective_value_j"] = right_obj
  );
}

List child_objectives_numeric_flat(
    NumericVector z_num,
    double split_value,
    const arma::mat& Y_by_obs,
    const arma::vec& S_tot,
    const arma::vec& Q_tot,
    const std::vector<int>& offsets,
    int split_feat_effect = -1,
    int split_feat_pos = 0
) {
  const int N = Y_by_obs.n_cols;
  arma::vec SL(Y_by_obs.n_rows, arma::fill::zeros);
  arma::vec QL(Y_by_obs.n_rows, arma::fill::zeros);
  int NL = 0;
  for (int i = 0; i < N; ++i) {
    const double v = z_num[i];
    if (R_IsNA(v) || v > split_value) continue;
    ++NL;
    const arma::subview_col<double> yi = Y_by_obs.col(i);
    SL += yi;
    QL += yi % yi;
  }
  return child_objectives_from_flat_sums(SL, QL, S_tot, Q_tot, offsets, NL, N - NL,
    split_feat_effect, split_feat_pos);
}

List child_objectives_categorical_flat(
    IntegerVector z_fac,
    int level_index,
    const arma::mat& Y_by_obs,
    const arma::vec& S_tot,
    const arma::vec& Q_tot,
    const std::vector<int>& offsets
) {
  const int N = Y_by_obs.n_cols;
  arma::vec SL(Y_by_obs.n_rows, arma::fill::zeros);
  arma::vec QL(Y_by_obs.n_rows, arma::fill::zeros);
  int NL = 0;
  for (int i = 0; i < N; ++i) {
    if (z_fac[i] == NA_INTEGER || z_fac[i] != level_index + 1) continue;
    ++NL;
    const arma::subview_col<double> yi = Y_by_obs.col(i);
    SL += yi;
    QL += yi % yi;
  }
  // TODO: as in the categorical objective above, no split-feature grid halving is applied here
  // (default -1). The categorical splitting feature's own effect is treated like the others.
  return child_objectives_from_flat_sums(SL, QL, S_tot, Q_tot, offsets, NL, N - NL);
}

// -----------------------------------------------------------------------------
// search_best_split_point_cpp_internal
// Purpose:
//   Find best split for a single feature (categorical or numerical).
//   Accepts preprocessed Ym and S_tot for efficiency.
// Inputs:
//   z: Feature vector (numeric or categorical)
//   Ym: Preprocessed effect matrices (NaN replaced with 0)
//   S_tot: Column sums per matrix
//   n_quantiles, is_categorical, min_node_size
// Output:
//   List: split_point, split_objective
// -----------------------------------------------------------------------------
List search_best_split_point_cpp_internal(
    SEXP              z,
    const arma::mat&  Y_by_obs,              // effect rows as contiguous observation columns
    const arma::vec&  S_tot,                 // & = reference, const = read-only, avoids copying large vector
    const arma::vec&  Q_tot,                 // & = reference, const = read-only, avoids copying large vector
    const std::vector<int>& offsets,         // Column offsets for per-effect child objectives
    Nullable<int>     n_quantiles   = R_NilValue,
    bool              is_categorical = false,
    int               min_node_size  = 1,
    bool              compute_child_objectives = true,
    int               split_feat_effect = -1,      // effect index of the split feature's own
    const std::vector<double>& split_feat_grid = std::vector<double>())  // effect (-1 = none)
{
  const int Ly = offsets.size() - 1; // p
  const int N = Y_by_obs.n_cols; // n
  const int split_feat_n_grid = split_feat_effect >= 0 ?
    offsets[split_feat_effect + 1] - offsets[split_feat_effect] : 0;
  // The split feature's own effect (if any) has its ICE grid divided across the children;
  // its flattened column range comes straight from offsets.
  const bool is_own_effect_halved = split_feat_effect >= 0 &&
    split_feat_n_grid > 0 &&
    split_feat_grid.size() == static_cast<size_t>(split_feat_n_grid);
  const arma::uword sf_start = is_own_effect_halved ? offsets[split_feat_effect] : 0;      // a
  const arma::uword sf_end   = is_own_effect_halved ? offsets[split_feat_effect + 1] : 0;  // b
  // Split feature's parent total sum-of-squares, computed once here (before evaluating any
  // split candidate). WHY it is subtracted: for every OTHER feature the SS term cancels
  // (SS_parent = SS_left + SS_right, so it drops as a constant and only S terms remain). For
  // the split feature the grid is halved, so its SS does NOT cancel and must be kept -- but
  // then its objective carries an extra +SS_parent,j that differs per feature. Subtracting
  // sf_const removes that per-feature offset (leaving only the global constant Sum_k SS_k,
  // identical for every split choice) so split_objective stays comparable ACROSS features.
  const double sf_const = is_own_effect_halved ?
    arma::accu(Q_tot.subvec(sf_start, sf_end - 1)) : 0.0;

  double best_obj = R_PosInf, best_split = NA_REAL;
  std::string best_level;
  int best_level_idx = -1;
  NumericVector best_left_obj(Ly, NA_REAL);
  NumericVector best_right_obj(Ly, NA_REAL);

  /* ================= Categorical Variable Splitting ================= */
  if (is_categorical) {
    // Extract factor levels and convert to integer vector
    IntegerVector z_fac(z);
    CharacterVector lev = z_fac.attr("levels");
    const int K = lev.size();

    // If only one level, no valid split possible
    if (K <= 1)
      return List::create(_["split_point"] = R_NaString,
        _["split_objective"] = R_PosInf,
        _["left_objective_value_j"] = best_left_obj,
        _["right_objective_value_j"] = best_right_obj);

    // Calculate per-level left sums for the flattened effect matrix.
    std::vector<arma::vec> SumL(K);
    std::vector<int> countL(K, 0);
    for (int k = 0; k < K; ++k) {
      SumL[k] = arma::vec(Y_by_obs.n_rows, arma::fill::zeros);
    }

    // Accumulate sums for each level (skip NA)
    for (int i = 0; i < N; ++i) {
      if (z_fac[i] == NA_INTEGER) continue;
      int k = z_fac[i] - 1;
      ++countL[k];
      SumL[k] += Y_by_obs.col(i);  // Direct accumulation, NaN already processed in preprocessing
    }

    // Evaluate each level as potential split (one vs rest)
    for (int k = 0; k < K; ++k) {
      int NL = countL[k], NR = N - NL;

      // Skip if either child node would be too small
      if (NL < min_node_size || NR < min_node_size) continue;

      // Calculate objective function for this split
      // TODO: split-feature handling for a CATEGORICAL splitting feature is not done here.
      // For the numeric case, the split feature's own ICE grid is divided across the children,
      // so its risk uses the true child SSE over each surviving half minus its parent SS
      // (is_own_effect_halved branch below). The analogous treatment for a categorical splitting
      // feature is still open (same TODO as in the R reference search_best_split_point_pd.R);
      // currently its own effect is treated like any other feature (S-only, full grid).
      const arma::vec SL = SumL[k];
      const arma::vec SR = S_tot - SL;
      double obj = arma::accu( - SL%SL / NL - SR%SR / NR );

      // Update best split if this one is better
      if (obj < best_obj) {
        best_obj = obj;
        best_level = as<std::string>(lev[k]);
        best_level_idx = k;
      }
    }

    if (best_obj == R_PosInf)
      return List::create(_["split_point"]     = R_NaString,
        _["split_objective"] = R_PosInf,
        _["left_objective_value_j"] = best_left_obj,
        _["right_objective_value_j"] = best_right_obj);

    if (compute_child_objectives) {
      List child_obj = child_objectives_categorical_flat(z_fac, best_level_idx, Y_by_obs, S_tot, Q_tot, offsets);
      best_left_obj = child_obj["left_objective_value_j"];
      best_right_obj = child_obj["right_objective_value_j"];
    }
    return List::create(_["split_point"]     = best_level,
      _["split_objective"] = best_obj,
      _["left_objective_value_j"] = best_left_obj,
      _["right_objective_value_j"] = best_right_obj);
  }

  /* ================= Numerical Variable Splitting =================== */

  // Convert to numeric and sort the feature vector once
  NumericVector z_num;
  if (TYPEOF(z) == STRSXP) {
    // Convert character to numeric if possible
    CharacterVector z_char(z);
    z_num = NumericVector(z_char.size());
    for (int i = 0; i < z_char.size(); ++i) {
      if (z_char[i] == NA_STRING) {
        z_num[i] = NA_REAL;
      } else {
        double value = NA_REAL;
        if (!parse_double_strict(as<std::string>(z_char[i]), value)) {
          Rcpp::stop("search_best_split: character split feature contains non-numeric values.");
        }
        z_num[i] = value;
      }
    }
  } else {
    z_num = NumericVector(z);
  }

  // Create sorted index and sorted values
  IntegerVector ord = Rcpp::seq(0, N - 1);
  std::sort(ord.begin(), ord.end(), [&](int i, int j){ return z_num[i] < z_num[j]; });

  NumericVector z_sorted(N);
  for (int i = 0; i < N; ++i) z_sorted[i] = z_num[ord[i]];

  // Generate candidate split points
  NumericVector splits;
  if (n_quantiles.isNotNull()) {
    // Use quantile-based grid
    int nq = Rcpp::as<int>(n_quantiles);
    NumericVector uniq = unique_cpp(z_sorted);
    if (uniq.size() < nq) {
      splits = uniq;
    } else {
      // Calculate quantile positions
      NumericVector q(nq);
      for (int i = 0; i < nq; ++i) q[i] = (i + 1.0)/(nq + 1.0);
      splits = NumericVector(nq);
      for (int i = 0; i < nq; ++i)
        splits[i] = quantile_type7(z_sorted, q[i]);
      splits = unique_cpp(splits);
    }
  } else {
    // Use all unique values as candidates
    splits = unique_cpp(z_sorted);
  }

  // Check if we have any valid split candidates
  if (splits.size() == 0)
    return List::create(_["split_point"]     = NA_REAL,
      _["split_objective"] = R_PosInf,
      _["left_objective_value_j"] = best_left_obj,
      _["right_objective_value_j"] = best_right_obj);

  // Grid must be ascending for find_grid_interval; checked once here, not per candidate.
  if (is_own_effect_halved && !std::is_sorted(split_feat_grid.begin(), split_feat_grid.end())) {
    Rcpp::stop("search_best_split: split-feature ICE grid is not sorted ascending.");
  }

  // Stream through split candidates and accumulate left sums (incremental; splits are sorted)
  arma::vec SL(Y_by_obs.n_rows, arma::fill::zeros);
  // For the split feature, also track the left-child sum-of-squares over its own grid columns
  // (sf_M = its grid size). This mirrors the SL accumulation but is restricted to that one
  // feature's columns, so it adds O(sf_M) per candidate -- the same order as the existing
  // O(M) objective sweep, never a new asymptotic term even when the grid scales with n.
  const arma::uword sf_M = is_own_effect_halved ? (sf_end - sf_start) : 0;
  arma::vec QL_split(sf_M, arma::fill::zeros);

  int idx = 0;
  for (double sp : splits) {
    while (idx < N && z_sorted[idx] <= sp) {
      int r = ord[idx++];
      SL += Y_by_obs.col(r);
      if (is_own_effect_halved) {
        const double* yr = Y_by_obs.colptr(r) + sf_start;  // split feature's own rows for obs r
        for (arma::uword c = 0; c < sf_M; ++c) QL_split[c] += yr[c] * yr[c];
      }
    }
    int NL = idx, NR = N - NL;

    // Skip if either child node would be too small
    if (NL < min_node_size || NR < min_node_size) continue;

    // Calculate objective function for this split
    const arma::vec SR = S_tot - SL;
    double obj;
    if (!is_own_effect_halved) {
      obj = arma::accu( - SL%SL / NL - SR%SR / NR );
    } else {
      // The split feature's own ICE grid is divided between the children: the left child
      // keeps grid columns <= sp, the right child keeps columns > sp. Every other feature
      // keeps its full grid in both children.
      const arma::uword a = sf_start;   // split feature's first column
      const arma::uword b = sf_end;     // one past its last column
      const arma::uword M = SL.n_elem;
      const arma::uword mid = a + find_grid_interval(sp, split_feat_grid);  // first right-surviving col
      obj = 0.0;
      // Other features (columns outside the split feature's block): their sum-of-squares
      // cancels across split points, so only the S term is needed (both children).
      if (a > 0) {
        obj -= arma::dot(SL.head(a), SL.head(a)) / NL;
        obj -= arma::dot(SR.head(a), SR.head(a)) / NR;
      }
      if (b < M) {
        obj -= arma::dot(SL.tail(M - b), SL.tail(M - b)) / NL;
        obj -= arma::dot(SR.tail(M - b), SR.tail(M - b)) / NR;
      }
      // Split feature: the SS term does NOT cancel here, so use the true child SSE
      // (Q - S^2/n) over each surviving half. Subtracting sf_const (its parent SS) leaves the
      // score comparable across features. QL_split is indexed locally (column c -> c - a).
      if (mid > a) {  // left half, grid <= sp
        const arma::vec SL_l = SL.subvec(a, mid - 1);
        obj += arma::accu(QL_split.head(mid - a)) - arma::dot(SL_l, SL_l) / NL;
      }
      if (mid < b) {  // right half, grid > sp
        const arma::vec QR_r = Q_tot.subvec(mid, b - 1) - QL_split.subvec(mid - a, b - 1 - a);
        const arma::vec SR_r = SR.subvec(mid, b - 1);
        obj += arma::accu(QR_r) - arma::dot(SR_r, SR_r) / NR;
      }
      obj -= sf_const;
    }
    if (obj < best_obj) {
      best_obj = obj;
      best_split = sp;
    }
  }

  if (best_obj == R_PosInf || R_IsNA(best_split))
    return List::create(_["split_point"]     = NA_REAL,
      _["split_objective"] = R_PosInf,
      _["left_objective_value_j"] = best_left_obj,
      _["right_objective_value_j"] = best_right_obj);

  // Refine split point to midpoint between adjacent values
  double Lft = -std::numeric_limits<double>::infinity(),
    Rgt =  std::numeric_limits<double>::infinity();
  for (int i = 0; i < N; ++i) {
    double v = z_num[i];
    if (!R_IsNA(v) && v <= best_split && v > Lft) Lft = v;
    if (!R_IsNA(v) && v >  best_split && v < Rgt) Rgt = v;
  }
  double mid = std::isinf(Rgt) ? Lft : (Lft + Rgt) / 2.0;

  if (compute_child_objectives) {
    const int child_pos = is_own_effect_halved ? find_grid_interval(best_split, split_feat_grid) : 0;
    List child_obj = child_objectives_numeric_flat(z_num, best_split, Y_by_obs, S_tot, Q_tot, offsets,
      is_own_effect_halved ? split_feat_effect : -1, child_pos);
    best_left_obj = child_obj["left_objective_value_j"];
    best_right_obj = child_obj["right_objective_value_j"];
  }
  return List::create(_["split_point"]     = mid,
    _["split_objective"] = best_obj,
    _["left_objective_value_j"] = best_left_obj,
    _["right_objective_value_j"] = best_right_obj);
}

// -----------------------------------------------------------------------------
// search_best_split_cpp
// Purpose:
//   Evaluate all features in Z, find best split per feature, return full results.
// Inputs:
//   Z: DataFrame of features
//   Y: List of effect matrices
//   min_node_size, n_quantiles
// Output:
//   DataFrame: split_feature, is_categorical, split_point, split_objective, etc.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame search_best_split_cpp(
    DataFrame       Z,
    List            Y,
    int             min_node_size,
    Nullable<int>   n_quantiles = R_NilValue,
    double          active_effect_rel_tol = 1e-14)
{
  // Initialize output vectors
  const int p = Z.size();
  CharacterVector feat_names = Z.names();

  CharacterVector split_feature(p);
  LogicalVector   is_cat_vec(p);
  CharacterVector split_point_out(p);
  NumericVector   split_obj(p);

  // Preprocess and flatten all Y matrices' NaN values to avoid repeated processing.
  const int Ly = Y.size();
  std::vector<arma::mat> Ym(Ly);
  std::vector<int> offsets(Ly + 1, 0);
  std::vector<bool> active_effect(Ly, false);
  std::vector<double> effect_objective(Ly, 0.0);
  std::vector<arma::rowvec> S_effect(Ly);
  std::vector<arma::rowvec> Q_effect(Ly);
  int N = 0;
  for (int l = 0; l < Ly; ++l) {
    Ym[l] = arma_view(Y[l]);
    Ym[l].replace(arma::datum::nan, 0.0);  // Process all NaN values once
    if (l == 0) {
      N = Ym[l].n_rows;
    }
    S_effect[l] = arma::sum(Ym[l], 0);
    Q_effect[l] = arma::sum(Ym[l] % Ym[l], 0);
    effect_objective[l] = arma::accu(Q_effect[l] - S_effect[l]%S_effect[l] / N);
  }

  double total_effect_objective = 0.0;
  for (int l = 0; l < Ly; ++l) {
    if (effect_objective[l] > 0.0) {
      total_effect_objective += effect_objective[l];
    }
  }
  const double rel_tol = std::max(0.0, active_effect_rel_tol);
  const double active_effect_tol = total_effect_objective * rel_tol;

  int M = 0;
  int n_active_effects = 0;
  int single_active_effect = -1;
  for (int l = 0; l < Ly; ++l) {
    active_effect[l] = effect_objective[l] > active_effect_tol;
    if (active_effect[l]) {
      ++n_active_effects;
      single_active_effect = l;
    }
  }
  if (n_active_effects == 0 && Ly > 0) {
    int best_effect = 0;
    double best_effect_objective = effect_objective[0];
    for (int l = 1; l < Ly; ++l) {
      if (effect_objective[l] > best_effect_objective) {
        best_effect = l;
        best_effect_objective = effect_objective[l];
      }
    }
    active_effect[best_effect] = true;
    n_active_effects = 1;
    single_active_effect = best_effect;
  }
  for (int l = 0; l < Ly; ++l) {
    offsets[l] = M;
    if (active_effect[l]) {
      M += Ym[l].n_cols;
    }
  }
  offsets[Ly] = M;
  arma::mat Y_flat;
  const arma::mat* Y_search = nullptr;
  arma::rowvec S_tot(M, arma::fill::zeros);
  arma::rowvec Q_tot(M, arma::fill::zeros);
  if (n_active_effects == 1) {
    Y_search = &Ym[single_active_effect];
    S_tot = S_effect[single_active_effect];
    Q_tot = Q_effect[single_active_effect];
  } else {
    Y_flat.set_size(N, M);
    for (int l = 0; l < Ly; ++l) {
      if (active_effect[l]) {
        const int start = offsets[l];
        const int end = offsets[l + 1] - 1;
        Y_flat.cols(start, end) = Ym[l];
        S_tot.subvec(start, end) = S_effect[l];
        Q_tot.subvec(start, end) = Q_effect[l];
      }
    }
    Y_search = &Y_flat;
  }
  // Armadillo stores matrices column-major. Split scanning repeatedly accumulates
  // all effect columns for one observation, so scanning a transposed view gives
  // contiguous memory access for each observation vector.
  arma::mat Y_by_obs = Y_search->t();
  arma::vec S_tot_col = S_tot.t();
  arma::vec Q_tot_col = Q_tot.t();
  SEXP effect_names_obj = Y.attr("names");
  const bool has_effect_names = !Rf_isNull(effect_names_obj);
  CharacterVector effect_names = has_effect_names ? CharacterVector(effect_names_obj) : CharacterVector(0);
  List left_obj_list(p);
  List right_obj_list(p);
  NumericVector empty_left_obj(Ly, NA_REAL);
  NumericVector empty_right_obj(Ly, NA_REAL);
  if (has_effect_names) {
    empty_left_obj.attr("names") = effect_names;
    empty_right_obj.attr("names") = effect_names;
  }
  // Map each split feature (Z column) to its own (active) effect index and grid. When a
  // feature is split on, its ICE grid is divided across the children; the effect index locates
  // its columns in offsets, and the grid tells the point search where to divide.
  std::vector<int> feat_effect(p, -1);
  std::vector<std::vector<double>> feat_grid(p);
  if (has_effect_names) {
    for (int j = 0; j < p; ++j) {
      if (Rf_isFactor(Z[j])) continue;
      const std::string fname = as<std::string>(feat_names[j]);
      for (int l = 0; l < Ly; ++l) {
        if (!active_effect[l]) continue;
        if (as<std::string>(effect_names[l]) == fname) {
          feat_effect[j] = l;
          feat_grid[j] = grid_from_matrix_colnames(Y[l]);
          break;
        }
      }
    }
  }

  // Evaluate each feature
  for (int j = 0; j < p; ++j) {
    SEXP z_j  = Z[j];
    bool is_c = Rf_isFactor(z_j);

    // First pass: find the best split objective for each split feature.
    // Child objective vectors are only needed for the globally selected split and are computed below.
    List res = search_best_split_point_cpp_internal(
      z_j, Y_by_obs, S_tot_col, Q_tot_col, offsets, n_quantiles, is_c, min_node_size, false,
      feat_effect[j], feat_grid[j]);

    // Store results
    split_feature[j]   = feat_names[j];
    is_cat_vec[j]      = is_c;
    split_obj[j]       = res["split_objective"];
    split_point_out[j] = as<CharacterVector>(wrap(res["split_point"]))[0];
    left_obj_list[j] = Rcpp::clone(empty_left_obj);
    right_obj_list[j] = Rcpp::clone(empty_right_obj);
  }
  // Find best valid split.
  // The internal function always returns R_PosInf when no split is found, so
  // split_obj[j] < best_val (with best_val = R_PosInf) is sufficient: R_PosInf < R_PosInf
  // is false, correctly excluding invalid features without inspecting the split point string.
  double best_val = R_PosInf;
  int best_idx = -1;

  for (int j = 0; j < p; ++j) {
    if (split_obj[j] < best_val) {
      best_val = split_obj[j];
      best_idx = j;
    }
  }
  
  LogicalVector best_split(p, false);
  if (best_idx >= 0) {
    best_split[best_idx] = true;

    SEXP z_best = Z[best_idx];
    List best_res = search_best_split_point_cpp_internal(
      z_best, Y_by_obs, S_tot_col, Q_tot_col, offsets, n_quantiles, is_cat_vec[best_idx], min_node_size, true,
      feat_effect[best_idx], feat_grid[best_idx]);
    split_obj[best_idx] = best_res["split_objective"];
    split_point_out[best_idx] = as<CharacterVector>(wrap(best_res["split_point"]))[0];

    NumericVector left_obj = best_res["left_objective_value_j"];
    NumericVector right_obj = best_res["right_objective_value_j"];
    if (has_effect_names) {
      left_obj.attr("names") = effect_names;
      right_obj.attr("names") = effect_names;
    }
    left_obj_list[best_idx] = left_obj;
    right_obj_list[best_idx] = right_obj;
  }

  // Return results as DataFrame. Attach child objectives afterwards so they stay list-columns.
  DataFrame out = DataFrame::create(
    _["split_feature"]   = split_feature,
    _["is_categorical"]  = is_cat_vec,
    _["split_point"]     = split_point_out,
    _["split_objective"] = split_obj,
    _["best_split"]      = best_split
  );
  out["left_objective_value_j"] = left_obj_list;
  out["right_objective_value_j"] = right_obj_list;
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -p);
  return out;
}
