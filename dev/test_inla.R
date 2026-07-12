# =============================================================================
# INLA install smoke test
# =============================================================================
# Two fits on simulated spatial data with known truth:
#   (a) Gaussian SPDE -- the bread-and-butter use case
#   (b) Binomial SPDE -- the family we'll need for Stage 2
# We check that posterior summaries recover the true (range, sigma_m, beta).
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
})
INLA::inla.setOption(fmesher.evolution.warn = FALSE)

set.seed(42)
nu <- 1
sigma_m_true <- 0.8
range_true   <- 0.3
kappa_true   <- sqrt(8 * nu) / range_true
beta_true    <- c(intercept = -0.5, x1 = 1.2, x2 = -0.7)

n <- 400
locs <- cbind(runif(n), runif(n))
X    <- cbind(1, runif(n), runif(n))
colnames(X) <- c("intercept", "x1", "x2")

# True spatial field on observation points (via direct Matern, not SPDE)
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1; out
}
D <- as.matrix(dist(locs))
Sigma_z <- sigma_m_true^2 * matern_corr(D, kappa_true, nu) + 1e-6 * diag(n)
z_true  <- as.vector(t(chol(Sigma_z)) %*% rnorm(n))
eta_true <- as.vector(X %*% beta_true + z_true)

# Build mesh + SPDE once
mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(0.08, 0.25), cutoff = 0.03)
cat("Mesh nodes:", mesh$n, "\n")
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = c(0.1, 0.5),
                                  prior.sigma = c(1.0, 0.5))
A_obs <- INLA::inla.spde.make.A(mesh = mesh, loc = locs)

# Helper: pretty-print posterior summary row
print_summary <- function(label, fit) {
  cat("\n=== ", label, " ===\n", sep = "")
  cat("Fixed effects (truth, post.mean, 95% CI):\n")
  for (nm in colnames(X)) {
    if (nm %in% rownames(fit$summary.fixed)) {
      r <- fit$summary.fixed[nm, ]
      cat(sprintf("  %-9s truth=% .3f   mean=% .3f   CI=[% .3f, % .3f]\n",
                  nm, beta_true[nm], r$mean, r$`0.025quant`, r$`0.975quant`))
    }
  }
  cat("Spatial hyperparameters (truth, post.mean, 95% CI):\n")
  if (!is.null(fit$summary.hyperpar)) {
    hp <- fit$summary.hyperpar
    rng_row <- grep("^Range", rownames(hp), value = TRUE)
    sig_row <- grep("^Stdev", rownames(hp), value = TRUE)
    if (length(rng_row)) {
      r <- hp[rng_row, ]
      cat(sprintf("  range     truth=% .3f   mean=% .3f   CI=[% .3f, % .3f]\n",
                  range_true, r$mean, r$`0.025quant`, r$`0.975quant`))
    }
    if (length(sig_row)) {
      r <- hp[sig_row, ]
      cat(sprintf("  sigma_m   truth=% .3f   mean=% .3f   CI=[% .3f, % .3f]\n",
                  sigma_m_true, r$mean, r$`0.025quant`, r$`0.975quant`))
    }
  }
  cat("log mlik =", round(fit$mlik[1], 2), "\n")
}

# -----------------------------------------------------------------------------
# (a) Gaussian SPDE
# -----------------------------------------------------------------------------
sigma_e_true <- 0.3
y_gauss <- eta_true + rnorm(n, 0, sigma_e_true)

stk_g <- INLA::inla.stack(
  data    = list(y = y_gauss),
  A       = list(A_obs, 1),
  effects = list(i = 1:spde$n.spde,
                 data.frame(intercept = X[, "intercept"],
                            x1        = X[, "x1"],
                            x2        = X[, "x2"]))
)
t0 <- Sys.time()
fit_g <- INLA::inla(
  y ~ -1 + intercept + x1 + x2 + f(i, model = spde),
  data = INLA::inla.stack.data(stk_g),
  control.predictor = list(A = INLA::inla.stack.A(stk_g)),
  control.compute   = list(dic = FALSE, waic = FALSE),
  verbose = FALSE
)
t1 <- Sys.time()
cat("Gaussian fit:", format(t1 - t0), "\n")
print_summary("Gaussian SPDE", fit_g)

# -----------------------------------------------------------------------------
# (b) Binomial SPDE
# -----------------------------------------------------------------------------
n_trials <- as.integer(sample(5:15, n, replace = TRUE))
k_succ   <- rbinom(n, n_trials, plogis(eta_true))

stk_b <- INLA::inla.stack(
  data    = list(y = k_succ, n_trials = n_trials),
  A       = list(A_obs, 1),
  effects = list(i = 1:spde$n.spde,
                 data.frame(intercept = X[, "intercept"],
                            x1        = X[, "x1"],
                            x2        = X[, "x2"]))
)
t0 <- Sys.time()
fit_b <- INLA::inla(
  y ~ -1 + intercept + x1 + x2 + f(i, model = spde),
  family = "binomial",
  Ntrials = INLA::inla.stack.data(stk_b)$n_trials,
  data    = INLA::inla.stack.data(stk_b),
  control.predictor = list(A = INLA::inla.stack.A(stk_b), link = 1),
  control.compute   = list(dic = FALSE, waic = FALSE),
  verbose = FALSE
)
t1 <- Sys.time()
cat("\nBinomial fit:", format(t1 - t0), "\n")
print_summary("Binomial SPDE", fit_b)

cat("\n--- INLA smoke test complete ---\n")
