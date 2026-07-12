# =============================================================================
# Empirical verification of Theorem 4.9 (BARTSIMP-PG -> sum-of-kernels GP).
#
# What we test
# ------------
# Draw many realizations of eta(x, s) = f_BART(x) + z(s) from the BARTSIMP-PG
# PRIOR (no data) at a fixed set of test points, for several values of m. We
# verify three predictions:
#
#   (A) CLT visualization (Lemma 4.8). The marginal distribution of eta at a
#       single test point converges to N(0, k_BART(x,x) + k_Mat(s,s)) as
#       m -> infty. Diagnostic: histograms, QQ vs Normal, Shapiro-Wilk W.
#
#   (B) Covariance recovery (Theorem 4.9). Empirical pairwise covariance of
#       eta over test points should converge to k_BART + k_Mat. Diagnostic:
#       Frobenius-norm distance ||Sigma_emp(m) - Sigma_target||_F vs m, where
#       Sigma_target is built from k_BART (estimated from m=1000 BART-only
#       prior) plus k_Mat (analytical Matern at the test locations).
#
#   (C) Additivity (sum-of-kernels). Same covariance should equal the sum of
#       BART-only and Matern-only covariances. Diagnostic: scatter of the
#       three covariance matrices' entries against each other.
#
# Output: dev/22_gp_limit_check.pdf and dev/22_gp_limit_check.rds.
# Expected runtime: a few minutes.
# =============================================================================

# We need matern_corr from bartsimp_pg.R but NOT the sampler itself; the file
# is safe to source and only adds INLA imports.
source("dev/bartsimp_pg.R")

set.seed(2026)

# ---- BART prior sampling (self-contained; does not need rbart.R for this) --

sample_tree_prior <- function(X, alpha = 0.95, beta = 2, tau = 1.0,
                              cuts_per_var = 100L) {
  # Returns a (recursive) tree drawn from the Chipman-George-McCulloch tree
  # prior with depth-decay alpha*(1+depth)^(-beta) and N(0, tau^2) leaves.
  P <- ncol(X)
  ranges <- apply(X, 2, range)
  build <- function(depth) {
    p_split <- alpha * (1 + depth)^(-beta)
    if (runif(1) >= p_split) {
      return(list(is_leaf = TRUE, mu = rnorm(1, 0, tau)))
    }
    v <- sample.int(P, 1)
    cut <- runif(1, ranges[1, v], ranges[2, v])
    list(
      is_leaf = FALSE, var = v, val = cut,
      left  = build(depth + 1L),
      right = build(depth + 1L)
    )
  }
  build(0L)
}

eval_tree <- function(tree, X) {
  # Vectorized: walk each row down the tree, return numeric vec length nrow(X).
  out <- numeric(nrow(X))
  walk <- function(node, idx) {
    if (node$is_leaf) {
      out[idx] <<- node$mu
      return(invisible())
    }
    left_mask  <- X[idx, node$var] <= node$val
    walk(node$left,  idx[ left_mask])
    walk(node$right, idx[!left_mask])
  }
  walk(tree, seq_len(nrow(X)))
  out
}

sample_bart_prior_one <- function(X, m, alpha, beta, tau_total) {
  # tau_total = sd(f_BART(x)) target; per-leaf sd is tau_total/sqrt(m).
  tau <- tau_total / sqrt(m)
  out <- numeric(nrow(X))
  for (l in seq_len(m)) {
    tree <- sample_tree_prior(X, alpha, beta, tau)
    out  <- out + eval_tree(tree, X)
  }
  out
}

sample_bart_prior_many <- function(X, m, B, alpha = 0.95, beta = 2,
                                   tau_total = 1.0, verbose = FALSE) {
  out <- matrix(0, B, nrow(X))
  for (b in seq_len(B)) {
    if (verbose && (b %% 200 == 0))
      cat(sprintf("    [bart m=%d] sample %d/%d\n", m, b, B))
    out[b, ] <- sample_bart_prior_one(X, m, alpha, beta, tau_total)
  }
  out
}

# ---- Matern (analytical) at the test locations -----------------------------

matern_cov_at <- function(locs, sigma_m = 1.0, rho = 0.3, nu = 1) {
  D     <- as.matrix(dist(locs))
  kappa <- sqrt(8 * nu) / rho
  sigma_m^2 * matern_corr(D, kappa, nu)
}

sample_matern_prior_many <- function(locs, B, sigma_m = 1.0, rho = 0.3, nu = 1) {
  S <- matern_cov_at(locs, sigma_m, rho, nu) + 1e-8 * diag(nrow(locs))
  L <- t(chol(S))
  zeta <- matrix(rnorm(B * nrow(locs)), nrow = nrow(locs), ncol = B)
  t(L %*% zeta)  # B x n
}

# ---- Test points ----------------------------------------------------------

N_TEST <- 30L                       # number of test points (joint covariate+location)
P      <- 2L                        # covariate dim
D      <- 2L                        # spatial dim

X_test    <- matrix(runif(N_TEST * P), N_TEST, P)
locs_test <- matrix(runif(N_TEST * D), N_TEST, D)

# Prior hyperparameters: keep the BART target sd and the Matern sd at 1 each,
# so the sum-of-kernels diagonal is around 2.
ALPHA      <- 0.95
BETA       <- 2
TAU_TOTAL  <- 1.0                   # target sd of f_BART
SIGMA_M    <- 1.0
RHO        <- 0.3

M_GRID <- c(1L, 5L, 20L, 100L, 500L)
B      <- 2000L                     # prior samples per m

# Estimate the "true" BART kernel from a large-m run (the GP-limit baseline).
M_REF <- 1000L
B_REF <- 4000L

cat("=== Sampling reference BART prior at m =", M_REF, "(this is k_BART target) ===\n")
t0 <- Sys.time()
bart_ref <- sample_bart_prior_many(
  X = X_test, m = M_REF, B = B_REF,
  alpha = ALPHA, beta = BETA, tau_total = TAU_TOTAL,
  verbose = TRUE
)
cat(sprintf("Reference BART done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

k_BART_target <- cov(bart_ref)                              # empirical k_BART
k_MAT_target  <- matern_cov_at(locs_test, SIGMA_M, RHO, 1)  # analytical
k_JOINT_target <- k_BART_target + k_MAT_target              # sum-of-kernels

# ---- Sweep over m ---------------------------------------------------------

results <- list()
for (m in M_GRID) {
  cat(sprintf("\n=== Sampling m = %d ===\n", m))
  t0 <- Sys.time()
  bart_draws <- sample_bart_prior_many(
    X = X_test, m = m, B = B,
    alpha = ALPHA, beta = BETA, tau_total = TAU_TOTAL,
    verbose = TRUE
  )
  mat_draws  <- sample_matern_prior_many(
    locs = locs_test, B = B,
    sigma_m = SIGMA_M, rho = RHO, nu = 1
  )
  eta_draws  <- bart_draws + mat_draws

  Sigma_emp <- cov(eta_draws)
  Sigma_bart_emp <- cov(bart_draws)
  Sigma_mat_emp  <- cov(mat_draws)

  # Frobenius distances to targets
  fro_joint  <- sqrt(sum((Sigma_emp     - k_JOINT_target)^2))
  fro_bart   <- sqrt(sum((Sigma_bart_emp - k_BART_target)^2))
  fro_mat    <- sqrt(sum((Sigma_mat_emp  - k_MAT_target )^2))

  # Single-point Gaussianity diagnostics at the FIRST test point
  focal_eta  <- eta_draws[, 1]
  shapiro_p  <- if (B <= 5000) shapiro.test(focal_eta)$p.value else NA_real_
  ks_p       <- ks.test(focal_eta,
                        "pnorm",
                        mean = mean(focal_eta), sd = sd(focal_eta))$p.value

  results[[as.character(m)]] <- list(
    m              = m,
    elapsed_sec    = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    Sigma_emp      = Sigma_emp,
    Sigma_bart_emp = Sigma_bart_emp,
    Sigma_mat_emp  = Sigma_mat_emp,
    fro_joint      = fro_joint,
    fro_bart       = fro_bart,
    fro_mat        = fro_mat,
    focal_eta      = focal_eta,
    shapiro_p      = shapiro_p,
    ks_p           = ks_p
  )

  cat(sprintf("  m=%4d  elapsed=%.1fs  fro(eta - target)=%.3f  fro(bart - k_BART)=%.3f\n",
              m, results[[as.character(m)]]$elapsed_sec, fro_joint, fro_bart))
  cat(sprintf("  Gaussianity at focal: shapiro p=%.3f  ks p=%.3f\n",
              shapiro_p, ks_p))
}

# ---- Summary table --------------------------------------------------------

tab <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    m         = r$m,
    fro_joint = r$fro_joint,
    fro_bart  = r$fro_bart,
    fro_mat   = r$fro_mat,
    shapiro_p = r$shapiro_p,
    ks_p      = r$ks_p,
    stringsAsFactors = FALSE
  )
}))
cat("\n=== Summary across m ===\n")
print(tab, row.names = FALSE, digits = 3)

# ---- Figure ---------------------------------------------------------------

pdf("dev/22_gp_limit_check.pdf", width = 11, height = 10)
layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 7, 7,
                8, 8, 8),
              nrow = 4, byrow = TRUE),
       heights = c(1, 1, 1, 1))

# Row 1: histograms of eta at focal point for 3 m values
m_show <- c(1L, 20L, 500L)
xr <- range(unlist(lapply(m_show, function(m) results[[as.character(m)]]$focal_eta)))
xs <- seq(xr[1], xr[2], length.out = 200)
for (m in m_show) {
  r <- results[[as.character(m)]]
  hist(r$focal_eta, breaks = 50, probability = TRUE,
       main = sprintf("m = %d:  hist(eta_focal)", m),
       xlab = expression(eta(x^"*", s^"*")),
       col = "lightblue", border = "white", xlim = xr,
       ylim = c(0, max(0.6, 1/sqrt(2*pi*var(r$focal_eta)))))
  lines(xs, dnorm(xs, mean = mean(r$focal_eta), sd = sd(r$focal_eta)),
        col = "red", lwd = 2)
  legend("topright", legend = c("empirical", "N(mu_hat, sd_hat)"),
         col = c("lightblue", "red"), lwd = c(NA, 2), pch = c(15, NA),
         bty = "n", cex = 0.9)
}

# Row 2: QQ plots vs Normal for the same m values
for (m in m_show) {
  r <- results[[as.character(m)]]
  qqnorm(r$focal_eta, pch = 16, cex = 0.4, col = "#2b7bba",
         main = sprintf("m = %d:  QQ vs Normal", m))
  qqline(r$focal_eta, col = "red", lwd = 2)
}

# Row 3: Frobenius distance vs m, log-log
plot(M_GRID, sapply(M_GRID, function(m) results[[as.character(m)]]$fro_joint),
     log = "xy", type = "b", pch = 16, lwd = 2, col = "#1b9e77",
     xlab = "m (number of trees)", ylab = "||Sigma_emp - Sigma_target||_F",
     main = "Covariance convergence to sum-of-kernels target (Theorem 4.9)",
     ylim = range(c(
       sapply(M_GRID, function(m) results[[as.character(m)]]$fro_joint),
       sapply(M_GRID, function(m) results[[as.character(m)]]$fro_bart)
     )))
lines(M_GRID, sapply(M_GRID, function(m) results[[as.character(m)]]$fro_bart),
      type = "b", pch = 16, lwd = 2, col = "#d95f02")
# 1/sqrt(m) reference line, scaled to pass through the first point
ref0 <- results[[as.character(M_GRID[1])]]$fro_joint
lines(M_GRID, ref0 / sqrt(M_GRID / M_GRID[1]),
      lty = 3, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c("eta = BART + Matern  vs  k_BART + k_Mat",
                  "BART only            vs  k_BART (limit)",
                  expression(1/sqrt(m)~reference)),
       col = c("#1b9e77", "#d95f02", "gray40"),
       lwd = 2, lty = c(1, 1, 3), pch = c(16, 16, NA), bty = "n", cex = 0.9)

# Row 4: additivity scatter -- joint cov vs BART cov + Matern cov, at m = 500
r_top <- results[[as.character(M_GRID[length(M_GRID)])]]
plot(as.vector(r_top$Sigma_bart_emp + r_top$Sigma_mat_emp),
     as.vector(r_top$Sigma_emp),
     pch = 16, cex = 0.4, col = "#7570b3",
     xlab = "k_BART_emp(x,x') + k_Mat_analytical(s,s')",
     ylab = "Sigma_emp(eta)(x,s; x',s')",
     main = sprintf("Sum-of-kernels additivity at m = %d (Theorem 4.9)",
                    M_GRID[length(M_GRID)]))
abline(0, 1, col = "red", lty = 2, lwd = 2)
mtext("(points should lie on the 45-degree line)",
      side = 3, line = -1.5, cex = 0.8)

dev.off()

# ---- Save raw -------------------------------------------------------------

saveRDS(list(
  config = list(N_TEST = N_TEST, P = P, D = D,
                ALPHA = ALPHA, BETA = BETA, TAU_TOTAL = TAU_TOTAL,
                SIGMA_M = SIGMA_M, RHO = RHO,
                M_GRID = M_GRID, B = B, M_REF = M_REF, B_REF = B_REF),
  test_points = list(X = X_test, locs = locs_test),
  targets = list(k_BART_emp_at_M_REF = k_BART_target,
                 k_MAT_analytical    = k_MAT_target,
                 k_JOINT             = k_JOINT_target),
  results = results,
  summary = tab
), "dev/22_gp_limit_check.rds")

cat("\nSaved dev/22_gp_limit_check.pdf and dev/22_gp_limit_check.rds\n")
