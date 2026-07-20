### First tests of selective early stopping (GADGET 2.0, Method 1 = "plain_risk"), 26-07-13.
### Structure mirrors example_CB_25-10-12_toy_tests.R, but uses the current API
### (GadgetTree / PdStrategy / target_feature_name, and the gadget_improvements argument).

source("example_CB_CC_26-07-13_input.R")



## Test 1: single interaction x1:x3 (x2, x4 are pure noise features) -----------

set.seed(123)
n = 1000
dat1 = data.frame(
  x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
  x3 = runif(n, -1, 1), x4 = runif(n, -1, 1)
)
dat1$y = ifelse(dat1$x3 > 0, 3 * dat1$x1, -3 * dat1$x1) + dat1$x3 + rnorm(n, sd = 0.3)
eff1 = make_effect(dat1)

cat("\n===== Test 1: disabled (baseline) =====\n")
t1_off = fit_tree(dat1, eff1)
show_tree(t1_off)

cat("\n===== Test 1: plain_risk, tau = 0.05 =====\n")
t1_on = fit_tree(dat1, eff1, "plain_risk", tau = 0.05)
show_tree(t1_on)
cat("remaining_features per node:\n")
walk_remaining(t1_on$root)



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
eff2 = make_effect(dat2)

cat("\n===== Test 2: disabled (baseline) =====\n")
t2_off = fit_tree(dat2, eff2)
show_tree(t2_off)

for (tau in c(0.005, 0.05, 0.5)) {
  cat(sprintf("\n===== Test 2: plain_risk, tau = %s =====\n", tau))
  tr = fit_tree(dat2, eff2, "plain_risk", tau = tau)
  show_tree(tr)
  cat("remaining_features per node:\n")
  walk_remaining(tr$root)
}



## Test 3: effect-pruning threshold warning -----------------------------------
### prune_effects_for_split_search() drops every feature whose root risk falls below
### rel_tol * (total root risk). With p features the average share is only about 1/p,
### so a fixed rel_tol gets more aggressive as p grows. A warning is therefore issued
### once p * rel_tol >= 0.1 (e.g. p = 1000 at the default rel_tol = 1e-4).
### Here we force it with 5 features x rel_tol = 0.05 -> 0.25 >= 0.1.

cat("\n===== Test 3: effect-pruning threshold warning =====\n")
old_rel_tol = getOption("xplaineff.active_effect_rel_tol")
options(xplaineff.active_effect_rel_tol = 0.05)

tr_warn = withCallingHandlers(
  fit_tree(dat2, eff2, "plain_risk", tau = 0.05),
  warning = function(w) {
    cat("  caught warning: ", conditionMessage(w), "\n", sep = "")
    invokeRestart("muffleWarning")
  }
)

# The fit must still succeed: the early-stopping bookkeeping is aligned with the pruned
# feature set, so vecb_remaining_features matches the (shorter) pruned Y.
cat("fit succeeded; nodes =", nrow(tr_warn$extract_split_info()), "\n")
cat("features kept after pruning:", length(tr_warn$root$vecb_remaining_features), "\n")
cat("still interacting at root:",
  paste(names(tr_warn$root$vecb_remaining_features)[tr_warn$root$vecb_remaining_features],
    collapse = ","), "\n")

options(xplaineff.active_effect_rel_tol = old_rel_tol)



## Test 4: all early-stopping modes run and report the expected columns --------
### Smoke test over the implemented selective early stopping methods (Section 5.1):
###   - disabled          : no early stopping (baseline)
###   - "plain_risk"      : Method 1, absolute normalized risk (drops already at the root)
###   - "risk_reduction"  : Method 2, relative risk reduction (drops only after a split)
### Every mode must fit, and extract_split_info() must carry both the total and the
### remaining-only objective / relative improvement. An unknown method must be rejected.

cat("\n===== Test 4: modes, reporting columns, and validation =====\n")

report_cols = c("node_objective", "node_objective_remaining", "int_imp", "int_imp_remaining")
for (method in list(NULL, "plain_risk", "risk_reduction")) {
  tr = fit_tree(dat1, eff1, method, tau = 0.05)
  si = tr$extract_split_info()
  method_label = if (is.null(method)) "disabled" else method
  cols_present = all(report_cols %in% colnames(si))
  cat(sprintf("  %-16s nodes=%d  all report columns present: %s\n", method_label, nrow(si), cols_present))
}

# An unknown improvement method must be rejected rather than silently ignored.
unknown_rejected = tryCatch({
  fit_tree(dat1, eff1, "nonsense", tau = 0.05)
  FALSE
}, error = function(e) TRUE)
cat("  unknown method rejected:", unknown_rejected, "\n")
