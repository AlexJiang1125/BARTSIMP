#!/usr/bin/env Rscript
# test_depthfloor_unit.R — unit tests for the semi-dense depth-floor edit to
# dev/rbart.R. Verifies:
#   (T1) d0 = 0 reduces to legacy EXACTLY (byte-identical RNG stream + trees)
#   (T2) tree_complete(d0) builds a full binary tree of depth d0 (2^d0 leaves,
#        all coarse layers present)
#   (T3) ensemble_sweep with d0 > 0 NEVER prunes a node at depth <= d0 over many
#        sweeps (the coarse block is invariant)
#   (T4) the sweep runs without error and leaf mu's adapt under the floor
#
# Run from bartsimp_cluster/:
#   env -u MallocNanoZone /Library/Frameworks/R.framework/Resources/bin/Rscript \
#       <thispath>/test_depthfloor_unit.R
suppressPackageStartupMessages({})
CLUSTER <- "/Users/alexziyujiang/Documents/Claude/Projects/BARTSIMP_glm/bartsimp_cluster"
source(file.path(CLUSTER, "dev/rbart.R"))

pass <- 0L; fail <- 0L
ok <- function(cond, msg) {
  if (isTRUE(cond)) { cat(sprintf("  [PASS] %s\n", msg)); pass <<- pass + 1L }
  else              { cat(sprintf("  [FAIL] %s\n", msg)); fail <<- fail + 1L }
}

# ---- synthetic leaf data ----
set.seed(1)
n <- 400L; p <- 5L
X <- matrix(runif(n * p), n, p)
r <- 0.8 * (X[,1] - 0.5) + 0.5 * sin(2*pi*X[,3]) + rnorm(n, 0, 0.3)
Omega <- rep(4, n)
tau2  <- (3 / (2 * sqrt(20)))^2

cat("== T1: d0 = 0 is byte-identical to legacy ==\n")
# legacy = ensemble_init without d0 arg, sweep without d0 arg
set.seed(123)
ens_a <- ensemble_init(20L, init_val = 0)            # default d0 = 0
set.seed(123); .Random.seed -> s_a
sweep_a <- ensemble_sweep(ens_a, X, r, Omega, tau2, alpha = 0.95, beta = 2)
set.seed(123)
ens_b <- ensemble_init(20L, init_val = 0, d0 = 0L)   # explicit d0 = 0
sweep_b <- ensemble_sweep(ens_b, X, r, Omega, tau2, alpha = 0.95, beta = 2, d0 = 0L)
ok(identical(ens_a, ens_b), "ensemble_init d0=0 identical to default")
ok(isTRUE(all.equal(sweep_a$fit, sweep_b$fit)), "ensemble_sweep d0=0 fit identical")
ok(identical(sweep_a$accepts, sweep_b$accepts), "ensemble_sweep d0=0 accept count identical")

cat("== T2: tree_complete builds full binary tree to depth d0 ==\n")
for (d0 in 1:3) {
  tr <- tree_complete(d0, init_mu = 0)
  n_leaves <- sum(tr$is_leaf)
  max_depth <- max(tr$depth)
  all_leaf_depth <- all(tr$depth[tr$is_leaf] == d0)
  ok(n_leaves == 2^d0, sprintf("d0=%d -> 2^%d=%d leaves (got %d)", d0, d0, 2^d0, n_leaves))
  ok(max_depth == d0, sprintf("d0=%d -> max depth %d (got %d)", d0, d0, max_depth))
  ok(all_leaf_depth, sprintf("d0=%d -> every leaf at depth d0", d0))
}
ok(identical(tree_complete(0L), new_tree(init_mu = 0)), "tree_complete(0) == stump")

cat("== T3: floored sweeps never prune the coarse block (depth <= d0) ==\n")
for (d0 in 2:3) {
  set.seed(7)
  ens <- ensemble_init(20L, init_val = 0, d0 = d0)
  # confirm init: every tree has the complete coarse block
  init_ok <- all(vapply(ens, function(t) sum(t$is_leaf) == 2^d0, logical(1)))
  ok(init_ok, sprintf("d0=%d init: all trees complete to depth d0", d0))
  viol <- 0L
  for (it in 1:40) {
    ens <- ensemble_sweep(ens, X, r, Omega, tau2, alpha = 0.95, beta = 2, d0 = d0)$ensemble
    # invariant: every internal node at depth < d0 must still be internal, and
    # the set of nodes at depth <= d0 must form the complete block (>= 2^k nodes
    # per level k for k <= d0). Simplest robust check: no leaf exists at depth < d0.
    for (t in ens) {
      if (any(t$is_leaf & t$depth < d0)) viol <- viol + 1L
    }
  }
  ok(viol == 0L, sprintf("d0=%d: no leaf ever appeared at depth < d0 over 40 sweeps (viol=%d)", d0, viol))
}

cat("== T4: floored sweep adapts (fit changes, finite, leaf mu's move) ==\n")
set.seed(11)
ens <- ensemble_init(20L, init_val = 0, d0 = 2L)
fit0 <- ensemble_predict(ens, X)
for (it in 1:30) ens <- ensemble_sweep(ens, X, r, Omega, tau2, d0 = 2L)$ensemble
fit1 <- ensemble_predict(ens, X)
ok(all(is.finite(fit1)), "d0=2 fit is finite after 30 sweeps")
ok(mean(abs(fit1 - fit0)) > 1e-6, "d0=2 fit adapts away from init")
# directional sanity: fit should correlate positively with the signal r
ok(cor(fit1, r) > 0.2, sprintf("d0=2 fit tracks signal (cor=%.2f)", cor(fit1, r)))

cat(sprintf("\n==== UNIT TESTS: %d passed, %d failed ====\n", pass, fail))
if (fail > 0L) quit(status = 1L)
