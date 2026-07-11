#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

static inline double ranger_predict_one_tree(
    const double* X,
    int n,
    int row,
    int override_feature,
    double override_value,
    const std::vector<int>& left_child,
    const std::vector<int>& right_child,
    const std::vector<int>& split_var_ids,
    const std::vector<double>& split_values
) {
  int node = 0;
  while (true) {
    const int left = left_child[node];
    const int right = right_child[node];
    if (left == 0 && right == 0) {
      return split_values[node];
    }

    const int split_var = split_var_ids[node];
    const double x = split_var == override_feature ? override_value : X[row + n * split_var];
    node = x <= split_values[node] ? left : right;
  }
}

static std::vector<int> as_int_vector(const NumericVector& x) {
  const int n = x.size();
  std::vector<int> out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = static_cast<int>(x[i]);
  }
  return out;
}

static std::vector<double> as_double_vector(const NumericVector& x) {
  const int n = x.size();
  std::vector<double> out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = x[i];
  }
  return out;
}

// [[Rcpp::export]]
List ranger_pd_numeric_cpp(List forest, NumericMatrix X, IntegerVector feature_indices, List grids) {
  const int n = X.nrow();
  const int num_trees = as<int>(forest["num.trees"]);
  List child_node_ids = forest["child.nodeIDs"];
  List split_var_ids = forest["split.varIDs"];
  List split_values = forest["split.values"];
  std::vector<std::vector<int> > left_children(num_trees);
  std::vector<std::vector<int> > right_children(num_trees);
  std::vector<std::vector<int> > tree_split_var_ids(num_trees);
  std::vector<std::vector<double> > tree_split_values(num_trees);
  for (int t = 0; t < num_trees; ++t) {
    List tree_children = child_node_ids[t];
    left_children[t] = as_int_vector(tree_children[0]);
    right_children[t] = as_int_vector(tree_children[1]);
    tree_split_var_ids[t] = as_int_vector(split_var_ids[t]);
    tree_split_values[t] = as_double_vector(split_values[t]);
  }
  const double* x_ptr = REAL(X);
  const int n_features = feature_indices.size();
  List out(n_features);

  for (int f = 0; f < n_features; ++f) {
    const int feature_index = feature_indices[f];
    NumericVector grid = grids[f];
    const int grid_len = grid.size();
    NumericMatrix ice(n, grid_len);

    for (int g = 0; g < grid_len; ++g) {
      const double grid_value = grid[g];
      for (int i = 0; i < n; ++i) {
        double pred = 0.0;
        for (int t = 0; t < num_trees; ++t) {
          pred += ranger_predict_one_tree(
            x_ptr,
            n,
            i,
            feature_index,
            grid_value,
            left_children[t],
            right_children[t],
            tree_split_var_ids[t],
            tree_split_values[t]
          );
        }
        ice(i, g) = pred / num_trees;
      }
    }
    out[f] = ice;
  }

  return out;
}

// [[Rcpp::export]]
List ranger_ale_numeric_cpp(List forest, NumericMatrix X, IntegerVector feature_indices, List x_left, List x_right) {
  const int n = X.nrow();
  const int num_trees = as<int>(forest["num.trees"]);
  List child_node_ids = forest["child.nodeIDs"];
  List split_var_ids = forest["split.varIDs"];
  List split_values = forest["split.values"];
  std::vector<std::vector<int> > left_children(num_trees);
  std::vector<std::vector<int> > right_children(num_trees);
  std::vector<std::vector<int> > tree_split_var_ids(num_trees);
  std::vector<std::vector<double> > tree_split_values(num_trees);
  for (int t = 0; t < num_trees; ++t) {
    List tree_children = child_node_ids[t];
    left_children[t] = as_int_vector(tree_children[0]);
    right_children[t] = as_int_vector(tree_children[1]);
    tree_split_var_ids[t] = as_int_vector(split_var_ids[t]);
    tree_split_values[t] = as_double_vector(split_values[t]);
  }
  const double* x_ptr = REAL(X);
  const int n_features = feature_indices.size();
  List out(n_features);

  for (int f = 0; f < n_features; ++f) {
    const int feature_index = feature_indices[f];
    NumericVector lower_values = x_left[f];
    NumericVector upper_values = x_right[f];
    NumericVector lower(n);
    NumericVector upper(n);

    for (int i = 0; i < n; ++i) {
      double pred_lower = 0.0;
      double pred_upper = 0.0;
      const double lower_value = lower_values[i];
      const double upper_value = upper_values[i];
      for (int t = 0; t < num_trees; ++t) {
        pred_lower += ranger_predict_one_tree(
          x_ptr,
          n,
          i,
          feature_index,
          lower_value,
          left_children[t],
          right_children[t],
          tree_split_var_ids[t],
          tree_split_values[t]
        );
        pred_upper += ranger_predict_one_tree(
          x_ptr,
          n,
          i,
          feature_index,
          upper_value,
          left_children[t],
          right_children[t],
          tree_split_var_ids[t],
          tree_split_values[t]
        );
      }
      lower[i] = pred_lower / num_trees;
      upper[i] = pred_upper / num_trees;
    }
    out[f] = List::create(
      _["lower"] = lower,
      _["upper"] = upper
    );
  }

  return out;
}
