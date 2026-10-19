# Runs the two bike-sharing examples of the R Journal article with the CRAN release of
# xplaineff on the current platform and prints the split tables, first with the code as
# submitted (no explicit ranger seed) and then with the revised code (ranger seed 123,
# one thread). Used by .github/workflows/platform-check.yaml to compare platforms.
suppressPackageStartupMessages({
  library(xplaineff)
  library(mlr3)
  library(mlr3learners)
})
cat(R.version.string, "|", Sys.info()[["sysname"]], Sys.info()[["machine"]], "\n")
for (p in c("xplaineff", "mlr3", "mlr3learners", "ranger", "ISLR2")) {
  cat(p, as.character(packageVersion(p)), " ")
}
cat("\n")

# Minimal ranger-only demonstration from https://github.com/imbs-hl/ranger/issues/533:
# with seed = NULL ranger takes its seed from the R stream; the issue reports
# OOB error 10.77% on macOS and 7.69% on Windows for this call.
cat("== ranger issue #533 snippet (set.seed(1), seed = NULL, 1 tree)\n")
set.seed(1)
fit = ranger::ranger(Species ~ ., data = iris, seed = NULL, num.trees = 1,
  importance = "impurity", num.threads = 1)
cat("OOB error:", fit$prediction.error, "\n")
print(round(fit$variable.importance, 2))
cat("== ranger with set.seed(123) and seed = 123, 500 trees\n")
set.seed(123)
fit = ranger::ranger(Species ~ ., data = iris, seed = 123, num.trees = 500,
  importance = "impurity", num.threads = 1)
cat("OOB error:", fit$prediction.error, "\n")
print(round(fit$variable.importance, 2))

data("Bikeshare", package = "ISLR2")
set.seed(123)
bike = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000), ]
factor_features = c("season", "mnth", "weekday", "workingday", "holiday", "weathersit")
bike[factor_features] = lapply(bike[factor_features], as.factor)
bike_data = bike[, c("hr", "temp", "workingday", "season", "mnth", "day", "holiday", "weekday",
                     "weathersit", "atemp", "hum", "windspeed", "bikers")]
names(bike_data)[names(bike_data) == "bikers"] = "target"
effect_features = c("hr", "temp", "workingday", "season")
split_features = c("temp", "workingday", "season")
cols = c("id", "depth", "n_obs", "node_type", "split_feature", "split_value", "int_imp", "is_final")

run_examples = function(learner, label) {
  task = TaskRegr$new(id = "bike", backend = bike_data, target = "target")
  learner$train(task)
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 2, min_node_size = 50)
  tree$fit(data = bike_data, target_feature_name = "target", model = learner,
    feature_set = effect_features, split_feature = split_features, n_grid = 20)
  cat("\n==", label, "| Example 2 (PD tree)\n")
  print(tree$extract_split_info()[, cols])
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 2, impr_par = 0.01, min_node_size = 50)
  tree$fit(data = bike_data, target_feature_name = "target", model = learner,
    feature_set = effect_features, split_feature = split_features, n_intervals = 10)
  cat("==", label, "| Example 3 (ALE tree)\n")
  print(tree$extract_split_info()[, cols])
}

# As submitted: ranger draws its seed from the R random number stream after set.seed(123).
set.seed(123)
bike = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000), ]
run_examples(lrn("regr.ranger"), "submitted code, no explicit seed")

# Revised: explicit ranger seed and a single thread.
run_examples(lrn("regr.ranger", seed = 123, num.threads = 1), "revised code, seed = 123")
