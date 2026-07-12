# =============================================================================
# Prototype 2 (Stage 1 target) — PG + Matern GP + BART covariates
# =============================================================================
# Full Stage-1 sampler: Polya-Gamma + spatial Matern GP + BART sum-of-trees.
# All in pure R (uses dev/rbart.R for the BART machinery). No INLA, no SPDE.
#
# Model:
#   y_{ik} ~ Bernoulli(sigmoid(eta_i))
#   eta_i  = f_BART(x_i) + z_i
#   f_BART = sum_l g(x; T_l, mu_l)             Chipman BART prior
#   z      ~ N(0, sigma_m^2 * Matern(rho))
#   (rho, sigma_m) ~ Fuglstad PC prior
#
# Sweep:
#   1. Omega_i ~ PG(n_i, |eta_i|),  ybar* = (k_i - n_i/2)/Omega_i
#   2. BART:  pseudo_y = ybar* - z,  weights Omega   -> ensemble_sweep
#   3. z:     pseudo_y = ybar* - f_BART,             Gaussian closed-form
#   4. (sigma_m, rho)  RW-MH on (log, log)
# =============================================================================

suppressPackageStartupMessages({
  library(BayesLogit)
  library(Matrix)
})
source("dev/rbart.R")

set.seed(123)

# -----------------------------------------------------------------------------
# Matern correlation, nu = 1
# -----------------------------------------------------------------------------
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1
  out
}

# -----------------------------------------------------------------------------
# Simulate spatial binary data with tree-structured f(x)
# -----------------------------------------------------------------------------
n_cl         <- 200
nu           <- 1
sigma_m_true <- 0.8
rho_true     <- 0.3
kappa_true   <- sqrt(8 * nu) / rho_true

# Covariate surface (tree-structured, à la the paper's sim study)
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
cat("sd(f_true) =", round(sd(f_true), 3),
    "  sd(z_true) =", round(sd(z_true), 3), "\n\n")

# -----------------------------------------------------------------------------
# Priors / hyperparams
# -----------------------------------------------------------------------------
# BART
m_trees <- 20
# Set tau so total prior SD of sum-of-trees is ~ half the plausible range of
# f_BART (which is roughly the same magnitude as plausible eta).  Plausible
# eta ~ [-3, 3], so half-range = 3.  Standard BART k = 2.
k_bart <- 2
tau    <- 3 / (k_bart * sqrt(m_trees))    # ~ 0.335 with m=20
tau2   <- tau^2
alpha_bart <- 0.95
beta_bart  <- 2
cat("BART tau =", round(tau, 3), "\n")

# PC prior on (rho, sigma_m)
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

make_Sigma_z <- function(sigma_m, rho) {
  kap <- sqrt(8 * nu) / rho
  sigma_m^2 * matern_corr(D, kap, nu) + 1e-6 * diag(n_cl)
}
log_p_z <- function(z, sigma_m, rho) {
  S <- make_Sigma_z(sigma_m, rho)
  U <- chol(S)
  q <- forwardsolve(t(U), z)
  -sum(log(diag(U))) - 0.5 * sum(q * q)
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
bart_accepts  <- integer(n_iter)
accept_psi    <- logical(n_iter)

sd_log_sigma <- 0.15
sd_log_rho   <- 0.15

f_current <- rep(0, n_cl)        # current sum-of-trees fit at each cluster
t0 <- Sys.time()
for (it in seq_len(n_iter)) {

  # ---- 1. PG augmentation ----
  eta   <- f_current + z
  Omega <- rpg(num = n_cl, h = n_i, z = eta)
  ybar  <- (k_i - n_i / 2) / Omega

  # ---- 2. BART sweep on (pseudo_y = ybar - z, weights = Omega) ----
  pseudo_y_for_bart <- ybar - z
  sweep <- ensemble_sweep(ensemble, X, pseudo_y_for_bart, Omega, tau2,
                          alpha = alpha_bart, beta = beta_bart)
  ensemble  <- sweep$ensemble
  f_current <- sweep$fit
  bart_accepts[it] <- sweep$accepts

  # ---- 3. z | rest  (Gaussian, dense; same as Proto 1) ----
  S      <- make_Sigma_z(sigma_m, rho)
  Q_psi  <- chol2inv(chol(S))
  Q_z    <- Q_psi + diag(Omega)
  b_vec  <- Omega * (ybar - f_current)
  U_z    <- chol(Q_z)
  mean_z <- backsolve(U_z, forwardsolve(t(U_z), b_vec))
  z      <- as.vector(mean_z + backsolve(U_z, rnorm(n_cl)))

  # ---- 4. (sigma_m, rho) RW-MH ----
  log_curr <- log_p_z(z, sigma_m, rho) + log_pc_prior(rho, sigma_m)
  sigma_m_p <- exp(log(sigma_m) + rnorm(1, 0, sd_log_sigma))
  rho_p     <- exp(log(rho)     + rnorm(1, 0, sd_log_rho))
  log_prop  <- log_p_z(z, sigma_m_p, rho_p) + log_pc_prior(rho_p, sigma_m_p)
  log_alpha <- log_prop - log_curr +
               (log(sigma_m_p) + log(rho_p)) -
               (log(sigma_m)   + log(rho))
  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    sigma_m <- sigma_m_p
    rho     <- rho_p
    accept_psi[it] <- TRUE
  }

  f_draws[it, ]     <- f_current
  sigma_m_draws[it] <- sigma_m
  rho_draws[it]     <- rho
  z_draws[it, ]     <- z

  if (it %% 100 == 0) {
    rmse_f <- sqrt(mean((f_current - f_true)^2))
    cat(sprintf("iter %4d  sigma_m=%5.3f  rho=%5.3f  RMSE(f)=%.3f  bart_acc=%d/%d  psi_acc=%.2f\n",
                it, sigma_m, rho, rmse_f, sweep$accepts, m_trees,
                mean(accept_psi[seq_len(it)])))
  }
}
t1 <- Sys.time()
cat("\nElapsed:", format(t1 - t0), "\n\n")

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
keep <- (burn + 1):n_iter
f_post <- colMeans(f_draws[keep, ])
z_post <- colMeans(z_draws[keep, ])

cat("Recovery:\n")
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

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
pdf("dev/03_pg_matern_bart_diagnostics.pdf", width = 9, height = 8)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

plot(sigma_m_draws, type = "l", ylab = "sigma_m", xlab = "iter",
     main = "sigma_m trace"); abline(h = sigma_m_true, col = "red", lty = 2)
plot(rho_draws, type = "l", ylab = "rho", xlab = "iter",
     main = "rho trace"); abline(h = rho_true, col = "red", lty = 2)
plot(f_true, f_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "f true", ylab = "f post mean", main = "BART surface recovery")
abline(0, 1, col = "red", lty = 2)
plot(z_true, z_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "z true", ylab = "z post mean", main = "spatial recovery")
abline(0, 1, col = "red", lty = 2)
plot(eta_true, f_post + z_post, pch = 20, col = rgb(0,0,0,0.5),
     xlab = "eta true", ylab = "eta post mean", main = "linear predictor recovery")
abline(0, 1, col = "red", lty = 2)

# Heatmap of f_post vs (x1, x2)
plot(x1, x2, col = ifelse(f_post > 0.5, "red",
                  ifelse(f_post < -0.5, "blue", "gray")),
     pch = 16, main = "BART surface (final state)")
abline(v = 0.5, lty = 3); abline(h = 0.5, lty = 3)

par(op); dev.off()
cat("\nDiagnostics saved to dev/03_pg_matern_bart_diagnostics.pdf\n")
