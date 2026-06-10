#!/usr/bin/env Rscript
# Quick smoke test: compare gadget PDP with and without ice packing
# Run from package root

Sys.setenv(
  OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", DATATABLE_NUM_THREADS = "1"
)

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: gadget") {
  stop("Run from GADGET package root")
}

pkgload::load_all(".", quiet = TRUE)
library(data.table)
setDTthreads(1L)
library(ranger)

cat("=== Loading data ===\n")
dat = fread("simulation/data/global_r_runtime/benchmark_N1000_D10_seed21.csv")
dat = as.data.frame(dat)

cat("=== Fitting RF ===\n")
model = ranger(y ~ ., data = dat, num.trees = 100L, mtry = 10L,
               min.node.size = 1L, replace = TRUE, sample.fraction = 1.0,
               splitrule = "variance", num.threads = 1L, seed = 21L)

rf_pred_fun = function(model, newdata) {
  as.numeric(predict(model, data = newdata, num.threads = 1L)$predictions)
}

features = setdiff(colnames(dat), "y")
n_grid = 20L
reps = 5L

cat("\n=== Version A: full calculate_pd (with ice packing) ===\n")
times_full = replicate(reps, {
  gc(FALSE)
  system.time(
    gadget:::calculate_pd(
      model = model, data = dat, target_feature_name = "y",
      feature_set = NULL, predict_fun = rf_pred_fun,
      n_grid = n_grid, pd_engine = "r"
    )
  )[["elapsed"]]
})
cat("Times:", paste(round(times_full, 3), collapse=", "), "s\n")
cat("Median:", round(median(times_full), 3), "s\n")

cat("\n=== Version B: ice_only (no packing) ===\n")
x_features_dt = data.table::as.data.table(dat[, features, drop = FALSE])
x_cols_list = as.list(x_features_dt)
n_obs = nrow(x_features_dt)

grids = stats::setNames(
  nm = features,
  lapply(features, function(f) gadget:::pd_feature_grid(x_features_dt[[f]], n_grid = n_grid))
)

max_g_r = max(lengths(grids))
stacked_pd_cache = list(
  stacked = data.table::as.data.table(lapply(x_features_dt, rep, times = max_g_r)),
  max_g = max_g_r,
  n_obs = n_obs
)

times_ice = replicate(reps, {
  gc(FALSE)
  system.time({
    for (feat in features) {
      grid = grids[[feat]]
      feat_index = match(feat, names(x_features_dt))
      ice_matrix = gadget:::compute_ice(
        model = model, data = x_features_dt, feature = feat, grid = grid,
        predict_fun = rf_pred_fun, pd_engine = "r",
        base_data_dt = x_features_dt, cols_list = x_cols_list,
        feature_index = feat_index, stacked_pd_cache = stacked_pd_cache
      )
    }
  })[["elapsed"]]
})
cat("Times:", paste(round(times_ice, 3), collapse=", "), "s\n")
cat("Median:", round(median(times_ice), 3), "s\n")

cat("\n=== Speedup ===\n")
speedup = median(times_full) / median(times_ice)
cat("Version B is", round(speedup, 3), "x faster than A\n")
cat("Packing overhead:", round((speedup - 1) * 100, 1), "%\n")
