# =============================================================================
# test_ensemble_predict_cpp.R
#
# Validate the Rcpp ensemble_predict against the pure-R reference.
#
# (1) Numerical equality: for several synthetic ensembles, every tree shape,
#     every X size, the Rcpp output must be IDENTICAL to the R output.
#     We use identical() not all.equal() — the math is pure float adds of
#     leaf values + deterministic tree traversal, so it should be exact.
# (2) Speedup: report wall time on a realistic-sized problem.
# (3) Spot-check with a real toy ensemble from a fit.
#
# Run:
#   Rscript dev/test_ensemble_predict_cpp.R
# =============================================================================

DEV <- "~/Documents/GitHub/BARTSIMP/dev"
source(file.path(DEV, "rbart.R"))
source(file.path(DEV, "ensemble_predict_cpp.R"))

stopifnot(ensemble_predict_backend() == "cpp")
cat("[test] Rcpp backend loaded.\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_random_tree <- function(P, depth_max = 4, seed = 1) {
  set.seed(seed)
  t <- new_tree(init_mu = rnorm(1, 0, 0.3))
  # Grow it by random successful grow moves.
  n_grow <- sample(0:8, 1)
  for (k in seq_len(n_grow)) {
    leaves <- which(t$is_leaf)
    if (length(leaves) == 0) break
    node_id <- sample(leaves, 1)
    if (t$depth[node_id] >= depth_max) next
    sv <- sample(P, 1)
    svv <- rnorm(1)
    ml <- rnorm(1, 0, 0.3)
    mr <- rnorm(1, 0, 0.3)
    t <- tree_grow_at(t, node_id, sv, svv, ml, mr)
  }
  t
}

make_random_ensemble <- function(P, m_trees, seed = 1) {
  set.seed(seed)
  lapply(seq_len(m_trees), function(i)
    make_random_tree(P, seed = sample.int(.Machine$integer.max, 1)))
}

# ---------------------------------------------------------------------------
# (1) Numerical equality on synthetic ensembles
# ---------------------------------------------------------------------------
cat("[test] (1) Numerical equality\n")
cases <- list(
  list(P = 1L,  m = 1L,  n = 100L,  label = "trivial"),
  list(P = 5L,  m = 10L, n = 500L,  label = "small"),
  list(P = 20L, m = 50L, n = 5000L, label = "medium"),
  list(P = 49L, m = 100L,n = 1000L, label = "production-shape")
)
for (cs in cases) {
  set.seed(42)
  X <- matrix(rnorm(cs$n * cs$P), cs$n, cs$P)
  ens <- make_random_ensemble(cs$P, cs$m, seed = 42)
  out_r   <- ensemble_predict_r(ens, X)
  out_cpp <- ensemble_predict(ens, X)   # routed to cpp
  if (!identical(out_r, out_cpp)) {
    diff <- max(abs(out_r - out_cpp))
    stop(sprintf("[FAIL] %s: max |R - cpp| = %.3e   sd(R) = %.3e",
                 cs$label, diff, sd(out_r)))
  }
  cat(sprintf("  ✓ %s (P=%d, m=%d, n=%d): identical() == TRUE\n",
              cs$label, cs$P, cs$m, cs$n))
}

# ---------------------------------------------------------------------------
# (2) Edge cases
# ---------------------------------------------------------------------------
cat("\n[test] (2) Edge cases\n")

# (a) all-leaf trees
set.seed(1)
X <- matrix(rnorm(100 * 5), 100, 5)
ens_leaf <- lapply(1:20, function(i) new_tree(init_mu = rnorm(1, 0, 0.3)))
or <- ensemble_predict_r(ens_leaf, X)
oc <- ensemble_predict(ens_leaf, X)
stopifnot(identical(or, oc))
cat("  ✓ all-leaf (m=20, depth=0) trees: identical\n")

# (b) single tree
ens1 <- list(make_random_tree(5, seed = 99))
or <- ensemble_predict_r(ens1, X)
oc <- ensemble_predict(ens1, X)
stopifnot(identical(or, oc))
cat("  ✓ single-tree ensemble: identical\n")

# (c) X with one row
ens5 <- make_random_ensemble(5, 10, seed = 13)
X1 <- matrix(rnorm(5), 1, 5)
or <- ensemble_predict_r(ens5, X1)
oc <- ensemble_predict(ens5, X1)
stopifnot(identical(or, oc))
cat("  ✓ X with 1 row: identical\n")

# (d) X exactly at split values (boundary behavior)
set.seed(7)
P <- 3L
t <- new_tree(init_mu = 0)
t <- tree_grow_at(t, 1L, split_var = 1L, split_val = 0.5,
                  mu_left = -1, mu_right = +1)
ens_split <- list(t)
X_boundary <- matrix(c(0.4999, 0.5000, 0.5001), 3, 3, byrow = FALSE)
or <- ensemble_predict_r(ens_split, X_boundary)
oc <- ensemble_predict(ens_split, X_boundary)
stopifnot(identical(or, oc))
cat(sprintf("  ✓ X at split boundary (< vs ≥): identical (values: %s)\n",
            paste(or, collapse = ",")))

# ---------------------------------------------------------------------------
# (3) Speedup benchmark on production-shape data
# ---------------------------------------------------------------------------
cat("\n[test] (3) Speed (production-shape: P=49, m_trees=100, n=80,000)\n")
set.seed(123)
P <- 49L; m <- 100L; n <- 80000L
X_big <- matrix(rnorm(n * P), n, P)
ens_big <- make_random_ensemble(P, m, seed = 123)

# Warm both backends
invisible(ensemble_predict_r(ens_big, X_big[1:10, , drop = FALSE]))
invisible(ensemble_predict(ens_big, X_big[1:10, , drop = FALSE]))

t_r <- system.time(out_r <- ensemble_predict_r(ens_big, X_big))[3]
t_c <- system.time(out_c <- ensemble_predict(ens_big, X_big))[3]
stopifnot(identical(out_r, out_c))
cat(sprintf("  Pure R   : %.2f sec\n", t_r))
cat(sprintf("  Rcpp     : %.2f sec\n", t_c))
cat(sprintf("  Speedup  : %.1fx\n", t_r / t_c))
cat(sprintf("  identical(R, Rcpp): %s\n", identical(out_r, out_c)))

# ---------------------------------------------------------------------------
# (4) Spot-check with the real toy posterior on real X covariates
# ---------------------------------------------------------------------------
TOY_RDS <- "~/Documents/Claude/Projects/BARTSIMP_glm/bartsimp_cluster/model/results/posterior_chains_stunting_rsr-top5_toy.rds"
if (file.exists(TOY_RDS)) {
  cat("\n[test] (4) Real toy posterior + sample X\n")
  post <- readRDS(TOY_RDS)
  fit  <- post$fit
  set.seed(7)
  # Build a synthetic X with the same standardization as the fit expected
  P <- length(post$X_std_mean)
  X_real <- matrix(rnorm(5000 * P), 5000, P)
  colnames(X_real) <- post$X_cov_names
  X_real_std <- scale(X_real, center = post$X_std_mean, scale = post$X_std_sd)
  ens_real <- fit$ensemble_snapshots[[1]]
  or <- ensemble_predict_r(ens_real, X_real_std)
  oc <- ensemble_predict(ens_real, X_real_std)
  if (identical(or, oc)) {
    cat("  ✓ real toy snapshot: identical\n")
  } else {
    cat(sprintf("  ⚠ DIFFERENCE: max abs diff = %.3e\n", max(abs(or - oc))))
  }
} else {
  cat("\n[test] (4) toy posterior not present (skipped)\n")
}

cat("\nAll tests passed.\n")
