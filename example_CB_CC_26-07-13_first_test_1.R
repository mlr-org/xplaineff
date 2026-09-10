### First tests of selective early stopping for GADGET 2.0, 26-07-13.
### Structure mirrors example_CB_25-10-12_toy_tests.R from before.

source("example_CB_CC_26-07-13_input.R")

cat("\nTook ~ 2 min on 23rd Aug 26\n")



## Test 1: single interaction x1:x3 (x2, x4 are pure noise features) -----------
## This is the first main example from the GADGET paper (Section 4.3 and Figure 3).
# In this test: Highest verbosity level for debugging purposes

set.seed(12345)
n = 1000
dat1 = data.frame(
  x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
  x3 = runif(n, -1, 1), x4 = runif(n, -1, 1)
)
dat1$y = ifelse(dat1$x3 > 0, 3 * dat1$x1, -3 * dat1$x1) + dat1$x3 + dat1$x1 - 2 * dat1$x4 + rnorm(n, sd = 0.3)
effect1 = make_effect(dat1)

cat("\n===== Test 1: Early stopping disabled (baseline) =====\n")
tree1_off = fit_tree(dat1, effect1, verbose = 3)
# show_tree(tree1_off)
print(tree1_off$extract_split_info())
cat("remaining_features per node:\n")
walk_remaining(tree1_off$root)

for (improvement_method in early_stopping_methods) {
  cat(sprintf("\n===== Test 1: Early stopping Method: %s, tau = 0.05 =====\n", improvement_method))
  tree1_on = fit_tree(dat1, effect1, improvement_method, tau = 0.05, verbose = 3)
  # show_tree(tree1_on)
  print(tree1_on$extract_split_info())
  cat("remaining_features per node:\n")
  walk_remaining(tree1_on$root)
}



## Test 2: two independent interactions x1:x3 and x2:x4 (x5 pure noise) --------

set.seed(123)
n = 1000
dat2 = data.frame(
  x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1),
  x4 = runif(n, -1, 1), x5 = runif(n, -1, 1)
)
dat2$y = ifelse(dat2$x3 > 0, 3 * dat2$x1, -3 * dat2$x1) +
  ifelse(dat2$x4 > 0, 2 * dat2$x2, -2 * dat2$x2) + dat2$x3 -
  dat2$x2 + 0.5 * dat2$x1 + 2 * dat2$x4 + rnorm(n, sd = 0.3)
effect2 = make_effect(dat2)

# cat("\n===== Test 2: disabled (baseline) =====\n")
# tree2_off = fit_tree(dat2, effect2)
# show_tree(tree2_off)
# cat("remaining_features per node:\n")
# walk_remaining(tree2_off$root)

# for (improvement_method in early_stopping_methods) {
#   for (tau in c(0.005, 0.05, 0.5)) {
#     cat(sprintf("\n===== Test 2: %s, tau = %s =====\n", improvement_method, tau))
#     tree = fit_tree(dat2, effect2, improvement_method, tau = tau)
#     show_tree(tree)
#     cat("remaining_features per node:\n")
#     walk_remaining(tree$root)
#   }
# }







# -------- Sanity checks and corner case tests ------------------------



## Test 3: effect-pruning threshold warning -----------------------------------
### prune_effects_for_split_search() drops every feature whose root risk falls below
### rel_tol * (total root risk). With p features the average share is only about 1/p,
### so a fixed rel_tol gets more aggressive as p grows. A warning is therefore issued
### once p * rel_tol >= 0.1 (e.g. p = 1000 at the default rel_tol = 1e-4).
### Here we force it with 5 features and rel_tol = 0.05 -> p * rel_tol = 0.25 >= 0.1.

cat("\n===== Tests 3 and 4: Smoke test / validation for single methods of selective early-stopping, testing effect-pruning at root, and reporting columns in extract_split_info() =====\n")

# Save before the loop in order to restore later
old_rel_tol = getOption("xplaineff.active_effect_rel_tol")
options(xplaineff.active_effect_rel_tol = 0.05)

# Columns to be tested for test 4
report_cols = c("node_objective", "node_objective_remaining", "int_imp", "int_imp_remaining")

for (improvement_method in early_stopping_methods) {
  cat(sprintf("\n===== Test 3: effect-pruning threshold warning, Early stopping Method: %s =====\n", improvement_method))

  tree_warn = withCallingHandlers(
    fit_tree(dat2, effect2, improvement_method, tau = 0.05),
    warning = function(w) {
      cat("  caught warning: ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    }
  )

  # The fit must still succeed: the early-stopping bookkeeping is aligned with the pruned
  # feature set, so vecb_remaining_features matches the (shorter) pruned Y.
  cat("fit succeeded; nodes =", nrow(tree_warn$extract_split_info()), "\n")
  cat("features kept after pruning:", length(tree_warn$root$vecb_remaining_features), "\n")
  cat("still interacting at root:",
    paste(names(tree_warn$root$vecb_remaining_features)[tree_warn$root$vecb_remaining_features],
      collapse = ","), "\n")



## Test 4: all early-stopping methods run and report the expected columns --------
### Smoke test over the implemented selective early stopping methods:
###   - disabled                    : no early stopping (baseline)
###   - "plain_risk"                : Method 1, absolute normalized risk
###   - "risk_reduction"            : Method 2, relative risk reduction (drops earliest after first split)
###   - "interaction_fraction"      : Method 3, per-feature interaction share
###   - "interaction_fraction_total": Method 4, interaction share against the total variance
### Every method must fit, and extract_split_info() must carry both the total and the
### remaining-only objective / relative improvement. An unknown method must be rejected.

  cat(sprintf("\n===== Test 4: Validation and checking reporting columns for improvement method: %s =====\n", improvement_method))
  split_info = tree$extract_split_info()
  method_label = if (is.null(method)) "disabled" else method
  cols_present = all(report_cols %in% colnames(split_info))
  cat(sprintf("  %-28s nodes=%d  all report columns present: %s\n", method_label, nrow(split_info), cols_present))
}

options(xplaineff.active_effect_rel_tol = old_rel_tol)

# For test 4: An unknown improvement method must be rejected rather than silently ignored.
unknown_rejected = tryCatch({
  fit_tree(dat1, effect1, "nonsense", tau = 0.05)
  FALSE
}, error = function(e) TRUE)
cat("  unknown method rejected:", unknown_rejected, "\n")
