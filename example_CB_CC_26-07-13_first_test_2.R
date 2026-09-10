### First tests / small benchmark of selective early stopping for GADGET 2.0, 26-07-13.
### Structure mirrors example_CB_25-10-12_toy_tests.R from before.

source("example_CB_CC_26-07-13_input.R")

cat("\nTook ~ 2 min on 23rd Aug 26\n")



## Test 2: two independent interactions x1:x3 and x2:x4 (x5 pure noise) --------

set.seed(123)
n = 2000
dat2 = data.frame(
  x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1),
  x4 = runif(n, -1, 1), x5 = runif(n, -1, 1)
)
dat2$y = ifelse(dat2$x3 > 0, 3 * dat2$x1, -3 * dat2$x1) +
  ifelse(dat2$x4 > 0, 2 * dat2$x2, -2 * dat2$x2) + dat2$x3 -
  dat2$x2 + 0.5 * dat2$x1 + 2 * dat2$x4 + rnorm(n, sd = 0.3)
effect2 = make_effect(dat2)

cat("\n===== Test 2: disabled (baseline) =====\n")
tree2_off = fit_tree(dat2, effect2, verbose = 3)
show_tree(tree2_off)
cat("remaining_features per node:\n")
walk_remaining(tree2_off$root)

for (improvement_method in early_stopping_methods) {
  tau_1 = 0.05
  cat(sprintf("\n===== Test 2: %s, tau = %s =====\n", improvement_method, tau_1))
  tree = fit_tree(dat2, effect2, improvement_method, tau = tau_1, verbose = 3)
  show_tree(tree)
  cat("remaining_features per node:\n")
  walk_remaining(tree$root)
}

# for (improvement_method in early_stopping_methods) {
#   for (tau in c(0.005, 0.05, 0.5)) {
#     cat(sprintf("\n===== Test 2: %s, tau = %s =====\n", improvement_method, tau))
#     tree = fit_tree(dat2, effect2, improvement_method, tau = tau)
#     show_tree(tree)
#     cat("remaining_features per node:\n")
#     walk_remaining(tree$root)
#   }
# }
