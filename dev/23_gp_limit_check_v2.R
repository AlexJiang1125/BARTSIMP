# =============================================================================
# Empirical verification of Theorem 4.9 (BARTSIMP-PG -> sum-of-kernels GP),
# revision 2.
#
# Why v2?
# -------
# v1 (dev/22) tested Gaussianity at a SINGLE test point. That's vacuous:
# the marginal of f_BART(x*) at any fixed x* is exactly N(0, tau^2) under
# the BARTSIMP-PG scaling, for every m -- because the leaf value mu is
# Gaussian regardless of which tree it comes from. The CLT (Lemma 4.8) is
# about the JOINT distribution at multiple points, not marginals. Likewise,
# Cov(eta(x_1, s_1), eta(x_2, s_2)) is exactly k_BART(x_1, x_2) + k_Mat(s_1, s_2)
# for every m (sum-of-independent-components rule); what changes with m is
# the SHAPE of the joint distribution, not the second moments.
#
# v2 tests the right thing:
#   (A) Bivariate scatter at (x_1, x_2) for a CLOSE pair vs a FAR pair, at
#       m in {1, 5, 20, 100, 500}. The close pair under small m falls in the
#       same leaf with non-negligible probability -> mixture of (Y = X) on
#       the diagonal and independent off-diagonal blob. As m grows the mixture
#       averages to a single ellipse.
#   (B) Multivariate normality test (energy::mvnorm.etest) on a 4-point joint
#       distribution, vs m. The test statistic should decay with m.
#   (C) QQ of a SUM statistic T = sum_i eta_i: by the joint CLT this converges
#       to Gaussian as m grows; non-Gaussian at small m.
#   (D) Sanity additivity check (covariance bookkeeping): Cov(eta) = Cov(BART)
#       + Cov(Matern) exactly, for all m. Should be 45-deg line.
# =============================================================================

source("dev/bartsimp_pg.R")
suppressPackageStartupMessages(library(energy))

set.seed(2026)

# ---- BART prior sampling (self-contained) ---------------------------------

sample_tree_prior <- function(X, alpha = 0.95, beta = 2, tau = 1.0) {
  P <- ncol(X)
  ranges <- apply(X, 2, range)
  build <- function(depth) {
    if (runif(1) >= alpha * (1 + depth)^(-beta)) {
      return(list(is_leaf = TRUE, mu = rnorm(1, 0, tau)))
    }
    v <- sample.int(P, 1)
    cut <- runif(1, ranges[1, v], ranges[2, v])
    list(is_leaf = FALSE, var = v, val = cut,
         left  = build(depth + 1L),
         right = build(depth + 1L))
  }
  build(0L)
}

eval_tree <- function(tree, X) {
  out <- numeric(nrow(X))
  walk <- function(node, idx) {
    if (node$is_leaf) { out[idx] <<- node$mu; return(invisible()) }
    left_mask <- X[idx, node$var] <= node$val
    walk(node$left,  idx[ left_mask])
    walk(node$right, idx[!left_mask])
  }
  walk(tree, seq_len(nrow(X)))
  out
}

sample_bart_prior <- function(X, m, B, alpha = 0.95, beta = 2, tau_total = 1.0,
                              verbose = FALSE) {
  tau <- tau_total / sqrt(m)
  out <- matrix(0, B, nrow(X))
  for (b in seq_len(B)) {
    if (verbose && b %% 500 == 0) cat(sprintf("    bart m=%d sample %d/%d\n", m, b, B))
    fb <- numeric(nrow(X))
    for (l in seq_len(m)) fb <- fb + eval_tree(sample_tree_prior(X, alpha, beta, tau), X)
    out[b, ] <- fb
  }
  out
}

matern_cov_at <- function(locs, sigma_m = 1.0, rho = 0.3, nu = 1) {
  D <- as.matrix(dist(locs))
  kappa <- sqrt(8 * nu) / rho
  sigma_m^2 * matern_corr(D, kappa, nu)
}

sample_matern_prior <- function(locs, B, sigma_m = 1.0, rho = 0.3, nu = 1) {
  S <- matern_cov_at(locs, sigma_m, rho, nu) + 1e-8 * diag(nrow(locs))
  L <- t(chol(S))
  t(L %*% matrix(rnorm(B * nrow(locs)), nrow(locs), B))
}

# ---- Test points: a CLOSE pair, a FAR pair, plus 6 random points ----------

# Build with controlled geometry. First 4 points are special:
#   (1, 2): a CLOSE pair in covariate space (x_2 very near x_1, same loc)
#   (1, 3): a FAR pair (very different x_3, same loc as x_1)
# The remaining 4 points are random; we use these for joint Mardia test.

X_test <- rbind(
  c(0.50, 0.50),   # x_1
  c(0.51, 0.49),   # x_2 -- VERY close to x_1
  c(0.10, 0.95),   # x_3 -- far from x_1
  c(0.85, 0.20),   # x_4 -- different again
  matrix(runif(6 * 2), 6, 2)   # 6 random
)
N_TEST <- nrow(X_test)

locs_test <- matrix(runif(N_TEST * 2), N_TEST, 2)

ALPHA     <- 0.95
BETA      <- 2
TAU_TOTAL <- 1.0
SIGMA_M   <- 1.0
RHO       <- 0.3

M_GRID <- c(1L, 5L, 20L, 100L)
B      <- 3000L

# ---- Sweep over m ---------------------------------------------------------

results <- list()
for (m in M_GRID) {
  cat(sprintf("\n=== Sampling m = %d ===\n", m))
  t0 <- Sys.time()

  bart_draws <- sample_bart_prior(X_test, m, B, ALPHA, BETA, TAU_TOTAL,
                                   verbose = TRUE)
  mat_draws  <- sample_matern_prior(locs_test, B, SIGMA_M, RHO, 1)
  eta_draws  <- bart_draws + mat_draws

  # (B) Multivariate normality on first 4 test points
  mvn_eta  <- energy::mvnorm.etest(eta_draws [, 1:4], R = 199)
  mvn_bart <- energy::mvnorm.etest(bart_draws[, 1:4], R = 199)

  # (D) Additivity scalar (should be ~ 0 always)
  cov_diff <- max(abs(cov(eta_draws) - cov(bart_draws) - cov(mat_draws)))

  # (C) Sum statistic for QQ
  T_sum <- rowSums(bart_draws[, 1:4])

  results[[as.character(m)]] <- list(
    m              = m,
    elapsed_sec    = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    bart_draws     = bart_draws,   # keep for plots
    eta_draws      = eta_draws,
    T_sum          = T_sum,
    mvn_eta_stat   = mvn_eta$statistic,
    mvn_eta_p      = mvn_eta$p.value,
    mvn_bart_stat  = mvn_bart$statistic,
    mvn_bart_p     = mvn_bart$p.value,
    cov_diff       = cov_diff
  )

  cat(sprintf("  m=%3d  elapsed=%.1fs  mvn(eta) stat=%.3f p=%.3f  mvn(bart) stat=%.3f p=%.3f  cov_diff=%.3e\n",
              m, results[[as.character(m)]]$elapsed_sec,
              mvn_eta$statistic, mvn_eta$p.value,
              mvn_bart$statistic, mvn_bart$p.value,
              cov_diff))
}

# ---- Summary --------------------------------------------------------------

tab <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    m              = r$m,
    elapsed_sec    = r$elapsed_sec,
    mvn_eta_stat   = r$mvn_eta_stat,
    mvn_eta_p      = r$mvn_eta_p,
    mvn_bart_stat  = r$mvn_bart_stat,
    mvn_bart_p     = r$mvn_bart_p,
    cov_diff       = r$cov_diff,
    stringsAsFactors = FALSE
  )
}))
cat("\n=== Summary across m ===\n")
print(tab, row.names = FALSE, digits = 3)

# ---- Figure ---------------------------------------------------------------

pdf("dev/23_gp_limit_check_v2.pdf", width = 12, height = 12)
layout(matrix(c(
  1, 2, 3, 4,    # row 1: scatter of BART (eta_1, eta_2) for CLOSE pair, m varies
  5, 6, 7, 8,    # row 2: scatter of BART (eta_1, eta_3) for FAR pair,   m varies
  9, 10, 11, 12, # row 3: QQ of T = sum of first 4 BART components
  13,13, 14,14   # row 4: MVN-test statistic vs m;  cov-additivity scatter (m=100)
), nrow = 4, byrow = TRUE), heights = c(1, 1, 1, 1))

# Row 1: close pair  (cols 1 vs 2 of BART)
xr_close <- c(-4, 4)
for (m in M_GRID) {
  r <- results[[as.character(m)]]
  plot(r$bart_draws[, 1], r$bart_draws[, 2],
       pch = 16, cex = 0.25, col = rgb(0.1, 0.5, 0.8, 0.3),
       xlim = xr_close, ylim = xr_close,
       xlab = "BART(x_1)", ylab = "BART(x_2)  (x_2 close to x_1)",
       main = sprintf("BART CLOSE-pair scatter, m=%d", m))
  abline(0, 1, lty = 2, col = "red", lwd = 1)
  abline(h = 0, v = 0, lty = 3, col = "gray60")
}

# Row 2: far pair (cols 1 vs 3 of BART)
for (m in M_GRID) {
  r <- results[[as.character(m)]]
  plot(r$bart_draws[, 1], r$bart_draws[, 3],
       pch = 16, cex = 0.25, col = rgb(0.6, 0.2, 0.4, 0.3),
       xlim = xr_close, ylim = xr_close,
       xlab = "BART(x_1)", ylab = "BART(x_3)  (x_3 far from x_1)",
       main = sprintf("BART FAR-pair scatter, m=%d", m))
  abline(0, 1, lty = 2, col = "red", lwd = 1)
  abline(h = 0, v = 0, lty = 3, col = "gray60")
}

# Row 3: QQ of sum statistic on BART
for (m in M_GRID) {
  r <- results[[as.character(m)]]
  qqnorm(r$T_sum, pch = 16, cex = 0.3, col = "#2b7bba",
         main = sprintf("QQ of T = sum_{i=1}^4 BART(x_i), m=%d", m))
  qqline(r$T_sum, col = "red", lwd = 2)
}

# Row 4 cell A: MVN test statistic vs m (log-log)
plot(M_GRID, sapply(M_GRID, function(m) results[[as.character(m)]]$mvn_bart_stat),
     log = "xy", type = "b", pch = 16, lwd = 2, col = "#d95f02",
     xlab = "m (number of trees)", ylab = "Energy MVN test statistic (on BART, 4 pts)",
     main = "Joint Gaussianity test: BART -> GP as m grows (Lemma 4.8)")
lines(M_GRID, sapply(M_GRID, function(m) results[[as.character(m)]]$mvn_eta_stat),
      type = "b", pch = 16, lwd = 2, col = "#1b9e77")
ref0 <- results[[as.character(M_GRID[1])]]$mvn_bart_stat
lines(M_GRID, ref0 / sqrt(M_GRID / M_GRID[1]),
      lty = 3, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c("BART only (joint at 4 pts)",
                  "eta = BART + Matern   (joint at 4 pts)",
                  expression(1/sqrt(m)~reference)),
       col = c("#d95f02", "#1b9e77", "gray40"),
       lwd = 2, lty = c(1, 1, 3), pch = c(16, 16, NA), bty = "n", cex = 0.85)

# Row 4 cell B: covariance additivity (m=100)
r_top <- results[[as.character(M_GRID[length(M_GRID)])]]
cov_eta  <- cov(r_top$eta_draws)
cov_bart <- cov(r_top$bart_draws)
cov_mat  <- matern_cov_at(locs_test, SIGMA_M, RHO, 1)   # analytical
plot(as.vector(cov_bart + cov_mat),
     as.vector(cov_eta),
     pch = 16, cex = 0.4, col = "#7570b3",
     xlab = "Cov(BART)_emp + Cov(Matern)_analytical",
     ylab = "Cov(eta)_emp",
     main = sprintf("Sum-of-kernels additivity sanity check (m=%d)",
                    M_GRID[length(M_GRID)]))
abline(0, 1, col = "red", lty = 2, lwd = 2)
mtext("(holds for ALL m by independence; serves as a bookkeeping check)",
      side = 3, line = -1.5, cex = 0.8)

dev.off()

saveRDS(list(
  config = list(M_GRID = M_GRID, B = B, ALPHA = ALPHA, BETA = BETA,
                TAU_TOTAL = TAU_TOTAL, SIGMA_M = SIGMA_M, RHO = RHO),
  test_points = list(X = X_test, locs = locs_test),
  results = results,
  summary = tab
), "dev/23_gp_limit_check_v2.rds")

cat("\nSaved dev/23_gp_limit_check_v2.{pdf,rds}\n")
