options(warn = -1)

required_packages = c("ISLR2", "mlr3", "mlr3learners", "ggplot2", "patchwork")
missing_packages = required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

devtools::load_all(".", quiet = TRUE)
library(ISLR2)
library(mlr3)
library(mlr3learners)

out_dir = file.path(
  "simulation", "results", "runtime_runs", "bikeshare_exhaustive_categorical_20260722"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

format_value = function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("")
  }
  if (is.numeric(x)) {
    return(format(round(x, 4L), scientific = FALSE, trim = TRUE))
  }
  as.character(x)
}

write_markdown_table = function(x, path) {
  cols = names(x)
  rows = vapply(seq_len(nrow(x)), function(i) {
    values = vapply(cols, function(col) format_value(x[[col]][i]), character(1L))
    paste0("| ", paste(values, collapse = " | "), " |")
  }, character(1L))
  lines = c(
    paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|"),
    rows
  )
  writeLines(lines, path)
}

collapse_levels = function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  paste(as.character(x), collapse = ", ")
}

condition_or_empty = function(child) {
  if (is.null(child) || is.null(child$parent) || is.null(child$parent$split_condition)) {
    return("")
  }
  child$parent$split_condition
}

split_table = function(tree, label) {
  nodes = unlist(tree$get_tree_list(), recursive = FALSE)
  rows = lapply(nodes, function(node) {
    if (is.null(node) || is.null(node$split)) {
      return(NULL)
    }
    data.frame(
      run = label,
      id = as.integer(node$id),
      depth = as.integer(node$depth),
      n_obs = length(node$subset_idx),
      split_feature = as.character(node$split$feature),
      split_value = as.character(node$split$value),
      split_levels = collapse_levels(node$split$levels),
      left_condition = condition_or_empty(node$children$left_child),
      right_condition = condition_or_empty(node$children$right_child),
      int_imp = if (is.null(node$importance)) NA_real_ else node$importance$imp,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  rownames(out) = NULL
  out
}

save_tree_plot = function(tree, filename) {
  png(file.path(out_dir, filename), width = 1100, height = 650)
  print(tree$plot_tree_structure(label_wrap_width = 38L, node_spread_x = 1.7, node_spread_y = 1.2))
  dev.off()
}

plot_file_node_id = function(node_name) {
  as.integer(sub("^Node_", "", node_name))
}

save_effect_plots = function(plot_list, prefix) {
  if (!length(plot_list)) {
    return(invisible(NULL))
  }
  for (depth_name in names(plot_list)) {
    nodes = plot_list[[depth_name]]
    if (!is.list(nodes) || !length(nodes)) next
    depth_id = as.integer(sub("^Depth_", "", depth_name))
    for (node_name in names(nodes)) {
      node_id = plot_file_node_id(node_name)
      filename = file.path(out_dir, sprintf("%s_depth%d_id%d.png", prefix, depth_id, node_id))
      png(filename, width = 1100, height = 650)
      print(nodes[[node_name]])
      dev.off()
    }
  }
  invisible(NULL)
}

set.seed(20260722)
data("Bikeshare")
bike = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000L), ]
factor_features = c("season", "mnth", "holiday", "weekday", "workingday", "weathersit")
bike[factor_features] = lapply(bike[factor_features], as.factor)
bike_data = bike[, c(
  "hr", "temp", "workingday", "season", "mnth", "holiday", "weekday",
  "weathersit", "atemp", "hum", "windspeed", "bikers"
)]
names(bike_data)[names(bike_data) == "bikers"] = "target"

effect_features = c("hr", "temp", "season", "mnth", "workingday", "weathersit")
split_features = c("workingday", "season", "mnth", "holiday", "weekday", "weathersit")
effect_plot_features = c("hr", "temp", "season", "mnth")

task = TaskRegr$new(id = "bike_exhaustive_probe", backend = bike_data, target = "target")
learner = lrn("regr.ranger", num.trees = 80L, num.threads = 1L)
learner$train(task)

fit_pd = function(strategy, label) {
  tree = GadgetTree$new(strategy = strategy, n_split = 2L, min_node_size = 50L, impr_par = 0)
  tree$fit(
    data = bike_data,
    target_feature_name = "target",
    model = learner,
    feature_set = effect_features,
    split_feature = split_features,
    n_grid = 12L
  )
  save_tree_plot(tree, paste0(label, "_tree_structure.png"))
  plots = tree$plot(
    data = bike_data,
    target_feature_name = "target",
    features = effect_plot_features,
    show_plot = FALSE
  )
  save_effect_plots(plots, paste0(label, "_effects"))
  tree
}

fit_ale = function(strategy, label) {
  tree = GadgetTree$new(strategy = strategy, n_split = 2L, min_node_size = 50L, impr_par = 0)
  tree$fit(
    data = bike_data,
    target_feature_name = "target",
    model = learner,
    feature_set = effect_features,
    split_feature = split_features,
    n_intervals = 8L
  )
  save_tree_plot(tree, paste0(label, "_tree_structure.png"))
  plots = tree$plot(
    data = bike_data,
    target_feature_name = "target",
    features = effect_plot_features,
    show_plot = FALSE,
    mean_center = TRUE
  )
  save_effect_plots(plots, paste0(label, "_effects"))
  tree
}

cat("Fitting PD default categorical split...\n")
pd_default = fit_pd(PdStrategy$new(), "pd_default")
cat("Fitting PD exhaustive categorical split...\n")
pd_exhaustive = fit_pd(
  PdStrategy$new(categorical_split = "exhaustive", max_exhaustive_levels = 12L),
  "pd_exhaustive"
)
cat("Fitting ALE default categorical split...\n")
ale_default = fit_ale(AleStrategy$new(), "ale_default")
cat("Fitting ALE exhaustive categorical split...\n")
ale_exhaustive = fit_ale(
  AleStrategy$new(categorical_split = "exhaustive", max_exhaustive_levels = 12L),
  "ale_exhaustive"
)

summary_table = rbind(
  split_table(pd_default, "pd_default_one_vs_rest"),
  split_table(pd_exhaustive, "pd_exhaustive"),
  split_table(ale_default, "ale_default_ordered_prefix"),
  split_table(ale_exhaustive, "ale_exhaustive")
)
write.csv(summary_table, file.path(out_dir, "split_summary.csv"), row.names = FALSE)
write_markdown_table(summary_table, file.path(out_dir, "split_summary.md"))

metadata = data.frame(
  key = c(
    "n_rows", "effect_features", "split_features", "effect_plot_features",
    "pd_n_grid", "ale_n_intervals", "max_exhaustive_levels"
  ),
  value = c(
    nrow(bike_data),
    paste(effect_features, collapse = ", "),
    paste(split_features, collapse = ", "),
    paste(effect_plot_features, collapse = ", "),
    "12",
    "8",
    "12"
  )
)
write.csv(metadata, file.path(out_dir, "metadata.csv"), row.names = FALSE)
write_markdown_table(metadata, file.path(out_dir, "metadata.md"))

cat("Done. Outputs written to ", out_dir, "\n", sep = "")
