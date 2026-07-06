# Generate README example outputs
# Run: Rscript scripts/generate_readme_examples.R

options(warn = -1)
dir.create("figures", showWarnings = FALSE)
dir.create(file.path("paper", "figures"), recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob("figures/ale_bike_depth*_node*.png"))
unlink(Sys.glob("figures/ale_bike_depth*_id*.png"))
unlink(Sys.glob("figures/pd_bike_depth*_node*.png"))
unlink(Sys.glob("figures/pd_bike_depth*_id*.png"))
unlink(Sys.glob("figures/split_info_*_bike.*"))

devtools::load_all(".", quiet = TRUE)
library(mlr3)
library(mlr3learners)
library(ISLR2)

format_readme_value = function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  if (is.numeric(x)) {
    if (isTRUE(abs(x) >= 1000)) {
      return(format(round(x), big.mark = "", scientific = FALSE, trim = TRUE))
    }
    return(format(round(x, 3L), nsmall = 0L, scientific = FALSE, trim = TRUE))
  }
  as.character(x)
}

write_markdown_table = function(x, path) {
  cols = names(x)
  rows = vapply(seq_len(nrow(x)), function(i) {
    values = vapply(cols, function(col) format_readme_value(x[[col]][i]), character(1L))
    paste0("| ", paste(values, collapse = " | "), " |")
  }, character(1L))
  lines = c(
    paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|"),
    rows
  )
  writeLines(lines, path)
}

paper_split_info_cols = c(
  "id", "depth", "n_obs", "node_type", "split_feature", "split_value", "int_imp", "is_final"
)

select_paper_split_info = function(x) {
  x[, paper_split_info_cols, drop = FALSE]
}

plot_file_node_id = function(node_name) {
  as.integer(sub("^Node_", "", node_name))
}

save_bike_plots = function(plot_list, prefix, output_dir = "figures", feature = NULL) {
  if (!length(plot_list)) {
    return(invisible(NULL))
  }
  feature_suffix = if (is.null(feature)) "" else paste0("_", feature)
  for (depth_name in names(plot_list)) {
    nodes = plot_list[[depth_name]]
    if (!is.list(nodes) || !length(nodes)) next
    depth_id = as.integer(sub("^Depth_", "", depth_name))
    for (node_name in names(nodes)) {
      node_id = plot_file_node_id(node_name)
      fname = file.path(
        output_dir,
        sprintf("%s_bike%s_depth%d_id%d.png", prefix, feature_suffix, depth_id, node_id)
      )
      png(fname, width = 800, height = 500)
      print(nodes[[node_name]])
      dev.off()
    }
  }
  invisible(NULL)
}

save_readme_plots = function(plot_list, prefix) {
  save_bike_plots(plot_list = plot_list, prefix = prefix, output_dir = "figures")
}

# ---- ALE + Bike (README example) ----
cat("\n=== ALE + Bike ===\n")
set.seed(123)
bike = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000), ]
bike$workingday = as.factor(bike$workingday)
bike_data = bike[, c("hr", "temp", "workingday", "bikers")]
names(bike_data)[names(bike_data) == "bikers"] = "target"

task = TaskRegr$new(id = "bike", backend = bike_data, target = "target")
learner = lrn("regr.ranger")
learner$train(task)

tree_ale_bike = GadgetTree$new(
  strategy = AleStrategy$new(),
  n_split = 2,
  impr_par = 0.01,
  min_node_size = 50
)
tree_ale_bike$fit(
  data = bike_data,
  target_feature_name = "target",
  model = learner,
  n_intervals = 10
)

split_ale_bike = tree_ale_bike$extract_split_info()
cat("Split info (ALE Bike):\n")
print(split_ale_bike)
sink("figures/split_info_ale_bike.txt")
print(split_ale_bike)
sink()
write_markdown_table(split_ale_bike, "figures/split_info_ale_bike.md")

# save tree structure plot
png("figures/ale_bike_tree_structure.png", width = 800, height = 500)
tree_ale_bike$plot_tree_structure()
dev.off()

# collect all regional ALE plots (depth x node)
pl_ale_bike_all = tree_ale_bike$plot(
  data = bike_data,
  target_feature_name = "target",
  features = c("hr"),
  mean_center = TRUE,
  show_plot = FALSE
)

# save per-node plots
save_readme_plots(pl_ale_bike_all, "ale")

# ---- PD + Bike (README example) ----
cat("\n=== PD + Bike ===\n")

tree_pd_bike = GadgetTree$new(strategy = PdStrategy$new(), n_split = 2, min_node_size = 50)
tree_pd_bike$fit(data = bike_data, target_feature_name = "target", model = learner, n_grid = 20L)

split_pd_bike = tree_pd_bike$extract_split_info()
cat("Split info (PD Bike):\n")
print(split_pd_bike)
sink("figures/split_info_pd_bike.txt")
print(split_pd_bike)
sink()
write_markdown_table(split_pd_bike, "figures/split_info_pd_bike.md")

# save tree structure plot
png("figures/pd_bike_tree_structure.png", width = 800, height = 500)
tree_pd_bike$plot_tree_structure()
dev.off()

# collect all regional PD/ICE plots (depth x node)
pl_pd_bike_all = tree_pd_bike$plot(
  data = bike_data,
  target_feature_name = "target",
  features = c("hr", "temp"),
  show_plot = FALSE
)

# save per-node plots
save_readme_plots(pl_pd_bike_all, "pd")

# ---- PD + Bike (paper example) ----
cat("\n=== PD + Bike (paper) ===\n")

set.seed(123)
bike_paper = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000), ]
factor_features = c("season", "mnth", "weekday", "workingday", "holiday", "weathersit")
bike_paper[factor_features] = lapply(bike_paper[factor_features], as.factor)
bike_paper_data = bike_paper[, c("hr", "temp", "workingday", "season",
  "mnth", "day", "holiday", "weekday",
  "weathersit", "atemp", "hum", "windspeed",
  "bikers")]
names(bike_paper_data)[names(bike_paper_data) == "bikers"] = "target"
effect_features_paper = c("hr", "temp", "workingday", "season")
split_features_paper = c("temp", "workingday", "season")

task_paper = TaskRegr$new(id = "bike_paper", backend = bike_paper_data, target = "target")
learner_paper = lrn("regr.ranger")
learner_paper$train(task_paper)

tree_pd_paper = GadgetTree$new(strategy = PdStrategy$new(), n_split = 2, min_node_size = 50)
tree_pd_paper$fit(
  data = bike_paper_data,
  target_feature_name = "target",
  model = learner_paper,
  feature_set = effect_features_paper,
  split_feature = split_features_paper,
  n_grid = 20L
)

split_pd_paper = select_paper_split_info(tree_pd_paper$extract_split_info())
cat("Split info (PD Bike paper):\n")
print(split_pd_paper, row.names = FALSE)
sink(file.path("paper", "figures", "split_info_pd_bike_paper.txt"))
print(split_pd_paper, row.names = FALSE)
sink()
write_markdown_table(split_pd_paper, file.path("paper", "figures", "split_info_pd_bike_paper.md"))

png(file.path("paper", "figures", "pd_bike_tree_structure.png"), width = 800, height = 500)
tree_pd_paper$plot_tree_structure()
dev.off()

pl_pd_paper_all = tree_pd_paper$plot(
  data = bike_paper_data,
  target_feature_name = "target",
  features = c("hr"),
  show_plot = FALSE
)

save_bike_plots(
  plot_list = pl_pd_paper_all,
  prefix = "pd",
  output_dir = file.path("paper", "figures"),
  feature = "hr"
)

tree_ale_paper = GadgetTree$new(
  strategy = AleStrategy$new(),
  n_split = 2,
  impr_par = 0.01,
  min_node_size = 50
)
tree_ale_paper$fit(
  data = bike_paper_data,
  target_feature_name = "target",
  model = learner_paper,
  feature_set = effect_features_paper,
  split_feature = split_features_paper,
  n_intervals = 10
)

split_ale_paper = select_paper_split_info(tree_ale_paper$extract_split_info())
cat("Split info (ALE Bike paper):\n")
print(split_ale_paper, row.names = FALSE)
sink(file.path("paper", "figures", "split_info_ale_bike_paper.txt"))
print(split_ale_paper, row.names = FALSE)
sink()
write_markdown_table(split_ale_paper, file.path("paper", "figures", "split_info_ale_bike_paper.md"))

png(file.path("paper", "figures", "ale_bike_tree_structure.png"), width = 800, height = 500)
tree_ale_paper$plot_tree_structure()
dev.off()

pl_ale_paper_all = tree_ale_paper$plot(
  data = bike_paper_data,
  target_feature_name = "target",
  features = c("hr"),
  mean_center = TRUE,
  show_plot = FALSE
)

save_bike_plots(
  plot_list = pl_ale_paper_all,
  prefix = "ale",
  output_dir = file.path("paper", "figures")
)

cat("\nDone. Outputs in figures/ and paper/figures/\n")
