# =============================================================================
# Prototype 0 — PG + Matern GP, intercept-only Bernoulli model
# =============================================================================
# Goal: validate the Polya-Gamma + spatial Gaussian conditional posterior math
# WITHOUT any BART, WITHOUT SPDE approximation. Direct dense Matern covariance.
# If this doesn't recover (beta_0, sigma_m, rho) on simulated data, nothing
# downstream will.
#
# Model:
#   y_{ik} ~ Bernoulli(sigmoid(eta_i)),  k = 1..n_i,  i = 1..n
#   eta_i  = beta_0 + z_i
#   z      ~ N(0, sigma_m^2 * Corr_Matern(rho))    on cluster locations
#   beta_0 ~ N(0, 100)                              diffuse
#   (rho, sigma_m) ~ Fuglstad PC prior:
#       P(rho < rho_0) = alpha_1,  P(sigma_m > sigma_0) = alpha_2
#
# Sweep cycle:
#   1. Omega_i ~ PG(n_i, |eta_i|)                    [cluster-collapsed]
#      pseudo-response  ybar_i^* = (k_i - n_i/2) / Omega_i
#   2. beta_0 | rest    Gaussian closed-form
#   3. z | rest         Gaussian closed-form (dense Q = Sigma_z^{-1} + diag(Omega))
#   4. (sigma_m, rho)   RW-MH on (log sigma_m, log rho) with PC prior
# =============================================================================

suppressPackageStartupMessages({
  library(BayesLogit)
  library(Matrix)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Matern correlation (closed form, nu = 1 fixed throughout BARTSIMP)
# -----------------------------------------------------------------------------
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1
  out
}

# -----------------------------------------------------------------------------
# Simulate data
# -----------------------------------------------------------------------------
n_cl          <- 200                  # number of clusters
nu            <- 1                    # Matern smoothness, fixed
sigma_m_true  <- 1.0
rho_true      <- 0.3
kappa_true    <- sqrt(8 * nu) / rho_true
beta_0_true   <- -0.5

locs <- cbind(runif(n_cl), runif(n_cl))   # uniform on [0,1]^2
n_i  <- as.numeric(sample(5:15, n_cl, replace = TRUE))  # rpg needs double, not int

D <- as.matrix(dist(locs))
Sigma_z_true <- sigma_m_true^2 * matern_corr(D, kappa_true, nu)
z_true <- as.vector(t(chol(Sigma_z_true + 1e-6 * diag(n_cl))) %*% rnorm(n_cl))

eta_true <- beta_0_true + z_true
p_true   <- plogis(eta_true)
y_list   <- mapply(function(ni, pi) rbinom(ni, 1, pi),
                   n_i, p_true, SIMPLIFY = FALSE)
k_i      <- vapply(y_list, sum, integer(1))

cat("Simulated:", n_cl, "clusters, total N =", sum(n_i),
    ", overall y mean =", round(mean(unlist(y_list)), 3), "\n")

# -----------------------------------------------------------------------------
# PC prior on (rho, sigma_m): log density
#   pi(rho, sigma) = (d/2) * lam1 * lam2 * rho^{-d/2 - 1}
#                    * exp(-lam1 rho^{-d/2} - lam2 sigma)
#   lam1 = -log(alpha_1) * rho_0^{d/2}
#   lam2 = -log(alpha_2) / sigma_0
# -----------------------------------------------------------------------------
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
# log p(z | sigma_m, rho)  for the MH step on the hyperparameters
# -----------------------------------------------------------------------------
make_Sigma_z <- function(sigma_m, rho) {
  kap <- sqrt(8 * nu) / rho
  sigma_m^2 * matern_corr(D, kap, nu) + 1e-6 * diag(n_cl)
}
log_p_z <- function(z, sigma_m, rho) {
  S <- make_Sigma_z(sigma_m, rho)
  U <- chol(S)                                       # S = U' U
  q <- forwardsolve(t(U), z)                         # U' q = z  =>  z' S^{-1} z = q' q
  -sum(log(diag(U))) - 0.5 * sum(q * q)
}

# -----------------------------------------------------------------------------
# MCMC
# -----------------------------------------------------------------------------
n_iter <- 3000
burn   <- 1000

beta_0  <- 0
z       <- rep(0, n_cl)
sigma_m <- 0.8
rho     <- 0.5

beta_0_draws  <- numeric(n_iter)
sigma_m_draws <- numeric(n_iter)
rho_draws     <- numeric(n_iter)
z_draws       <- matrix(NA_real_, n_iter, n_cl)
accept_psi    <- logical(n_iter)

# RW step sizes on (log sigma_m, log rho) — tune by acceptance rate ~ 0.2-0.4
sd_log_sigma <- 0.15
sd_log_rho   <- 0.15

t0 <- Sys.time()
for (it in seq_len(n_iter)) {

  # ---- 1. PG augmentation, cluster-collapsed ----
  eta   <- beta_0 + z
  Omega <- rpg(num = n_cl, h = n_i, z = eta)               # PG(n_i, |eta_i|)
  ybar  <- (k_i - n_i / 2) / Omega                          # pseudo-response

  # ---- 2. beta_0 | rest  (Gaussian closed-form) ----
  prec_b <- 1/100 + sum(Omega)
  mean_b <- sum(Omega * (ybar - z)) / prec_b
  beta_0 <- rnorm(1, mean_b, 1 / sqrt(prec_b))

  # ---- 3. z | rest  (dense Gaussian, Q = Sigma_z^{-1} + diag(Omega)) ----
  S      <- make_Sigma_z(sigma_m, rho)
  Q_psi  <- chol2inv(chol(S))                               # Sigma_z^{-1}
  Q_z    <- Q_psi + diag(Omega)
  b_vec  <- Omega * (ybar - beta_0)
  U      <- chol(Q_z)                                       # Q_z = U' U
  mean_z <- backsolve(U, forwardsolve(t(U), b_vec))
  z      <- as.vector(mean_z + backsolve(U, rnorm(n_cl)))

  # ---- 4. (sigma_m, rho) RW-MH on log scale ----
  log_curr <- log_p_z(z, sigma_m, rho) + log_pc_prior(rho, sigma_m)
  sigma_m_p <- exp(log(sigma_m) + rnorm(1, 0, sd_log_sigma))
  rho_p     <- exp(log(rho)     + rnorm(1, 0, sd_log_rho))
  log_prop  <- log_p_z(z, sigma_m_p, rho_p) + log_pc_prior(rho_p, sigma_m_p)
  log_alpha <- log_prop - log_curr +
               (log(sigma_m_p) + log(rho_p)) -
               (log(sigma_m)   + log(rho))                  # log-scale Jacobian
  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    sigma_m <- sigma_m_p
    rho     <- rho_p
    accept_psi[it] <- TRUE
  }

  # ---- Store ----
  beta_0_draws[it]  <- beta_0
  sigma_m_draws[it] <- sigma_m
  rho_draws[it]     <- rho
  z_draws[it, ]     <- z

  if (it %% 500 == 0) {
    cat(sprintf("iter %5d   beta_0=%6.3f  sigma_m=%5.3f  rho=%5.3f  acc=%.2f\n",
                it, beta_0, sigma_m, rho, mean(accept_psi[seq_len(it)])))
  }
}
t1 <- Sys.time()
cat("\nElapsed:", format(t1 - t0), "\n\n")

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
keep <- (burn + 1):n_iter
summ <- function(x, truth) {
  c(truth = truth,
    mean  = mean(x),
    q025  = unname(quantile(x, 0.025)),
    q500  = unname(quantile(x, 0.500)),
    q975  = unname(quantile(x, 0.975)))
}

cat("Recovery summary (post-burn, ", length(keep), " draws):\n", sep = "")
print(round(rbind(
  beta_0  = summ(beta_0_draws[keep],  beta_0_true),
  sigma_m = summ(sigma_m_draws[keep], sigma_m_true),
  rho     = summ(rho_draws[keep],     rho_true)
), 3))

z_post <- colMeans(z_draws[keep, ])
cat("\nz: cor(true, post.mean) =", round(cor(z_true, z_post), 3), "\n")
cat("   RMSE(z) =", round(sqrt(mean((z_true - z_post)^2)), 3), "\n")
cat("psi MH acceptance rate:", round(mean(accept_psi), 3), "\n")

# -----------------------------------------------------------------------------
# Trace + recovery plot
# -----------------------------------------------------------------------------
pdf(file.path("dev", "01_pg_matern_intercept_diagnostics.pdf"), width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(beta_0_draws, type = "l", ylab = "beta_0", xlab = "iter",
     main = "beta_0 trace"); abline(h = beta_0_true, col = "red", lty = 2)
plot(sigma_m_draws, type = "l", ylab = "sigma_m", xlab = "iter",
     main = "sigma_m trace"); abline(h = sigma_m_true, col = "red", lty = 2)
plot(rho_draws, type = "l", ylab = "rho", xlab = "iter",
     main = "rho trace"); abline(h = rho_true, col = "red", lty = 2)
plot(z_true, z_post, xlab = "z true", ylab = "z post mean",
     main = "z recovery", pch = 20, col = rgb(0,0,0,0.5))
abline(0, 1, col = "red", lty = 2)

par(op); dev.off()
cat("Diagnostics saved to dev/01_pg_matern_intercept_diagnostics.pdf\n")
