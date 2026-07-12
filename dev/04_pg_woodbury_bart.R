# =============================================================================
# Prototype 2b — PG + SPDE-Woodbury + BART
# =============================================================================
# Rewrite of Proto 2 using the same "extract Q, do Woodbury" trick the existing
# BARTSIMP C++ uses for its fast path. Goals:
#   - Validate the Woodbury identity gives sensible results in the PG context
#   - Verify the SPDE approximation agrees with the dense-Matern Proto 2
#   - Pin down all formulas so the C++ port is mechanical
#
# Notation:
#   A     : n_cl x n_mesh sparse projection matrix from mesh nodes to clusters
#   Q_psi : n_mesh x n_mesh sparse SPDE precision (depends on rho, sigma_m)
#   W     : diag(Omega) per-cluster precision after PG augmentation
#   Q_u   : A' W A + Q_psi    posterior precision over mesh nodes
#   yhat  : ybar* - f_BART_cluster
#
# Marginal log-likelihood of yhat after integrating out u:
#   yhat ~ N(0, W^{-1} + A Q_psi^{-1} A')                  (1)
# Woodbury identity:
#   Sigma_y^{-1} = W - W A Q_u^{-1} A' W                   (2)
# Sylvester determinant theorem:
#   log|Sigma_y| = log|Q_u| - log|W| - log|Q_psi|          (3)
# Quadratic form (with rhs = A' W yhat, m_u = Q_u^{-1} rhs):
#   yhat' Sigma_y^{-1} yhat = sum(Omega * yhat^2) - rhs' m_u   (4)
# So:
#   log p(yhat) = 0.5 [-N log(2 pi) + log|W| + log|Q_psi| - log|Q_u|
#                      - (sum(Omega yhat^2) - rhs' m_u)]
#
# Sampling z (= A u) given current state, using the same Cholesky of Q_u:
#   u | rest ~ N(m_u, Q_u^{-1})
#   draw u = m_u + L^{-T} eta,  eta ~ N(0, I),  where L L' = Q_u
#   then z = A u
#
# Note: the ψ-MH step now uses the MARGINAL likelihood (z integrated out),
# unlike Proto 2 which conditioned on the explicit z. This matches the
# INLA-within-MCMC structure of the published BARTSIMP and should mix better.
# =============================================================================

suppressPackageStartupMessages({
  library(BayesLogit)
  library(Matrix)
  library(INLA)
})
INLA::inla.setOption(fmesher.evolution.warn = FALSE)
source("dev/rbart.R")

set.seed(123)  # SAME SEED as Proto 2 so comparison is on identical data

# -----------------------------------------------------------------------------
# Simulate (identical to Proto 2)
# -----------------------------------------------------------------------------
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1; out
}

n_cl         <- 200
nu           <- 1
sigma_m_true <- 0.8
rho_true     <- 0.3
kappa_true   <- sqrt(8 * nu) / rho_true

x1 <- runif(n_cl); x2 <- runif(n_cl)
X  <- cbind(x1, x2)
colnames(X) <- c("x1", "x2")
f_true <- ifelse(x1 < 0.5,  1.5,
         ifelse(x2 < 0.5, -1.0,  0))

locs <- cbind(runif(n_cl), runif(n_cl))
n_i  <- as.numeric(sample(5:15, n_cl, replace = TRUE))

D <- as.matrix(dist(locs))
Sigma_z_true <- sigma_m_true^2 * matern_corr(D, kappa_true, nu)
z_true <- as.vector(t(chol(Sigma_z_true + 1e-6 * diag(n_cl))) %*% rnorm(n_cl))

eta_true <- f_true + z_true
p_true   <- plogis(eta_true)
y_list   <- mapply(function(ni, pi) rbinom(ni, 1, pi),
                   n_i, p_true, SIMPLIFY = FALSE)
k_i      <- vapply(y_list, sum, numeric(1))

cat("Simulated:", n_cl, "clusters, total N =", sum(n_i),
    ", overall y mean =", round(mean(unlist(y_list)), 3), "\n")

# -----------------------------------------------------------------------------
# Build mesh, SPDE object, projection A (all fixed across sweeps)
# -----------------------------------------------------------------------------
mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(0.08, 0.25), cutoff = 0.03)
cat("Mesh n =", mesh$n, "\n")
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = c(0.1, 0.5),
                                  prior.sigma = c(1.0, 0.5))
A <- INLA::inla.spde.make.A(mesh = mesh, loc = locs)
A <- as(A, "CsparseMatrix")

# -----------------------------------------------------------------------------
# Get sparse SPDE precision Q_psi for given (sigma_m, rho)
# Uses the SAME function the C++ uses (inla.spde.precision).
# theta convention: (log range, log sigma_m). range = sqrt(8 nu)/kappa.
# -----------------------------------------------------------------------------
make_Qpsi <- function(sigma_m, rho) {
  INLA::inla.spde.precision(spde, theta = c(log(rho), log(sigma_m)))
}

# -----------------------------------------------------------------------------
# Woodbury marginal log-likelihood + posterior factorization for z draw
#   yhat:  cluster-level residual = ybar* - f_BART
#   Omega: cluster-level PG sums
#   Returns mlik, m_u (posterior mean of u), Fchol (sparse Cholesky of Q_u)
# -----------------------------------------------------------------------------
margprob_woodbury <- function(yhat, Omega, A, Q_psi) {
  # A' diag(Omega) A   --- broadcasting Omega along rows of A
  AtWA  <- Matrix::crossprod(A, as.vector(Omega) * A)
  Qu    <- Matrix::forceSymmetric(AtWA + Q_psi)
  Fchol <- Matrix::Cholesky(Qu, LDL = FALSE, super = FALSE, perm = FALSE)
  rhs   <- as.numeric(Matrix::crossprod(A, Omega * yhat))
  m_u   <- as.numeric(Matrix::solve(Fchol, rhs, system = "A"))

  quad        <- sum(Omega * yhat^2) - sum(rhs * m_u)
  logdet_W    <- sum(log(Omega))
  logdet_Qpsi <- as.numeric(Matrix::determinant(Q_psi, logarithm = TRUE)$modulus)
  logdet_Qu   <- as.numeric(Matrix::determinant(Qu,    logarithm = TRUE)$modulus)
  N <- length(yhat)
  mlik <- 0.5 * (-N * log(2 * pi) + logdet_W + logdet_Qpsi - logdet_Qu - quad)

  list(mlik = mlik, m_u = m_u, Fchol = Fchol)
}

# Draw u ~ N(m_u, Q_u^{-1}) via the sparse Cholesky; project to clusters
draw_z_woodbury <- function(m_u, Fchol, A) {
  eta <- rnorm(length(m_u))
  u   <- m_u + as.numeric(Matrix::solve(Fchol, eta, system = "Lt"))
  as.numeric(A %*% u)
}

# -----------------------------------------------------------------------------
# Priors / hyperparams
# -----------------------------------------------------------------------------
m_trees    <- 20
k_bart     <- 2
tau        <- 3 / (k_bart * sqrt(m_trees))
tau2       <- tau^2
alpha_bart <- 0.95
beta_bart  <- 2

d_dim   <- 2
rho_0   <- 0.10
sigma_0 <- 2.0
alpha_1 <- 0.05
alpha_2 <- 0.05
lam1 <- -log(alpha_1) * rho_0^(d_dim/2)
lam2 <- -log(alpha_2) / sigma_0
log_pc_prior <- function(rho, sigma_m) {
  log(d_dim/2) + log(lam1) + log(lam2) -
    (d_dim/2 + 1) * log(rho) -
    lam1 * rho^(-d_dim/2) - lam2 * sigma_m
}

# -----------------------------------------------------------------------------
# MCMC
# -----------------------------------------------------------------------------
n_iter <- 1500
burn   <- 500

ensemble <- ensemble_init(m_trees, init_val = 0)
z        <- rep(0, n_cl)
sigma_m  <- 0.8
rho      <- 0.5

f_draws       <- matrix(NA_real_, n_iter, n_cl)
sigma_m_draws <- numeric(n_iter)
rho_draws     <- numeric(n_iter)
z_draws       <- matrix(NA_real_, n_iter, n_cl)
accept_psi    <- logical(n_iter)

sd_log_sigma <- 0.15
sd_log_rho   <- 0.15

f_current <- rep(0, n_cl)
t0 <- Sys.time()

for (it in seq_len(n_iter)) {
  # 1. PG augmentation
  eta   <- f_current + z
  Omega <- rpg(num = n_cl, h = n_i, z = eta)
  ybar  <- (k_i - n_i / 2) / Omega

  # 2. BART sweep (subtract z from pseudo-response so trees model f only)
  pseudo_y_for_bart <- ybar - z
  sweep <- ensemble_sweep(ensemble, X, pseudo_y_for_bart, Omega, tau2,
                          alpha = alpha_bart, beta = beta_bart)
  ensemble  <- sweep$ensemble
  f_current <- sweep$fit

  # 3. Hyperpar MH using Woodbury marginal likelihood (z integrated out)
  yhat   <- ybar - f_current
  Q_psi  <- make_Qpsi(sigma_m, rho)
  curr   <- margprob_woodbury(yhat, Omega, A, Q_psi)

  sigma_m_p <- exp(log(sigma_m) + rnorm(1, 0, sd_log_sigma))
  rho_p     <- exp(log(rho)     + rnorm(1, 0, sd_log_rho))
  Q_psi_p   <- make_Qpsi(sigma_m_p, rho_p)
  prop      <- margprob_woodbury(yhat, Omega, A, Q_psi_p)

  log_alpha <- (prop$mlik + log_pc_prior(rho_p, sigma_m_p)) -
               (curr$mlik + log_pc_prior(rho,   sigma_m))   +
               (log(sigma_m_p) + log(rho_p)) -
               (log(sigma_m)   + log(rho))
  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    sigma_m <- sigma_m_p
    rho     <- rho_p
    accepted <- prop
    accept_psi[it] <- TRUE
  } else {
    accepted <- curr
  }

  # 4. Refresh z reusing the accepted factorization
  z <- draw_z_woodbury(accepted$m_u, accepted$Fchol, A)

  f_draws[it, ]     <- f_current
  sigma_m_draws[it] <- sigma_m
  rho_draws[it]     <- rho
  z_draws[it, ]     <- z

  if (it %% 100 == 0) {
    rmse_f <- sqrt(mean((f_current - f_true)^2))
    cat(sprintf("iter %4d  sigma_m=%5.3f  rho=%5.3f  RMSE(f)=%.3f  psi_acc=%.2f\n",
                it, sigma_m, rho, rmse_f, mean(accept_psi[seq_len(it)])))
  }
}
t1 <- Sys.time()
cat("\nElapsed:", format(t1 - t0), "\n\n")

# -----------------------------------------------------------------------------
# Diagnostics + comparison to Proto 2 (dense)
# -----------------------------------------------------------------------------
keep <- (burn + 1):n_iter
f_post <- colMeans(f_draws[keep, ])
z_post <- colMeans(z_draws[keep, ])

cat("Recovery (Woodbury / SPDE form):\n")
cat(sprintf("  sigma_m: true=%.2f  post.mean=%.3f  CI=[%.3f, %.3f]\n",
            sigma_m_true, mean(sigma_m_draws[keep]),
            quantile(sigma_m_draws[keep], 0.025),
            quantile(sigma_m_draws[keep], 0.975)))
cat(sprintf("  rho:     true=%.2f  post.mean=%.3f  CI=[%.3f, %.3f]\n",
            rho_true, mean(rho_draws[keep]),
            quantile(rho_draws[keep], 0.025),
            quantile(rho_draws[keep], 0.975)))
cat(sprintf("  f(x):    cor(true, post)=%.3f  RMSE=%.3f\n",
            cor(f_true, f_post), sqrt(mean((f_true - f_post)^2))))
cat(sprintf("  z:       cor(true, post)=%.3f  RMSE=%.3f\n",
            cor(z_true, z_post), sqrt(mean((z_true - z_post)^2))))
cat(sprintf("  eta:     cor(true, post)=%.3f\n",
            cor(eta_true, f_post + z_post)))
cat(sprintf("  psi MH acceptance: %.3f\n", mean(accept_psi)))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
pdf("dev/04_pg_woodbury_bart_diagnostics.pdf", width = 9, height = 8)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

plot(sigma_m_draws, type = "l", ylab = "sigma_m", xlab = "iter",
     main = "sigma_m trace (Woodbury)"); abline(h = sigma_m_true, col = "red", lty = 2)
plot(rho_draws, type = "l", ylab = "rho", xlab = "iter",
     main = "rho trace (Woodbury)"); abline(h = rho_true, col = "red", lty = 2)
plot(f_true, f_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "f true", ylab = "f post mean", main = "BART surface recovery")
abline(0, 1, col = "red", lty = 2)
plot(z_true, z_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "z true", ylab = "z post mean", main = "spatial recovery")
abline(0, 1, col = "red", lty = 2)
plot(eta_true, f_post + z_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "eta true", ylab = "eta post mean", main = "linear predictor")
abline(0, 1, col = "red", lty = 2)
plot(x1, x2, col = ifelse(f_post > 0.5, "red",
                  ifelse(f_post < -0.5, "blue", "gray")),
     pch = 16, main = "BART surface (final state)")
abline(v = 0.5, lty = 3); abline(h = 0.5, lty = 3)

par(op); dev.off()
cat("\nDiagnostics saved to dev/04_pg_woodbury_bart_diagnostics.pdf\n")

# Save summary stats for downstream comparison vs Proto 2
saveRDS(list(
  sigma_m_draws = sigma_m_draws,
  rho_draws     = rho_draws,
  f_post        = f_post,
  z_post        = z_post,
  accept_psi    = accept_psi,
  elapsed       = as.numeric(difftime(t1, t0, units = "secs"))
), "dev/04_pg_woodbury_bart_results.rds")
