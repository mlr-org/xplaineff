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
