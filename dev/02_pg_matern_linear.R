# =============================================================================
# Prototype 1 — PG + Matern GP + linear covariates, Bernoulli model
# =============================================================================
# Extends Prototype 0 with a linear covariate term. eta = X beta + z, with the
# intercept folded into X (first column of 1s). beta is a p-vector, updated
# jointly via Gaussian-Gaussian conjugacy given Omega.
#
# Model:
#   y_{ik}  ~ Bernoulli(sigmoid(eta_i))
#   eta_i   = x_i' beta + z_i
#   beta    ~ N(0, V_0)                    (V_0 = 100 * I)
#   z       ~ N(0, sigma_m^2 * Matern(rho))
#   (rho, sigma_m) ~ PC prior
#
# Sweep cycle:
#   1. Omega_i ~ PG(n_i, |eta_i|)
#      ybar_i^* = (k_i - n_i/2) / Omega_i
#   2. beta | rest  ~ N(m_beta, Q_beta^{-1})  with
#        Q_beta = V_0^{-1} + X' diag(Omega) X
#        m_beta = Q_beta^{-1} X' diag(Omega) (ybar - z)
#   3. z | rest      Gaussian closed-form, as in Prototype 0
#                    (data mean is now ybar - X beta instead of ybar - beta_0)
#   4. (sigma_m, rho) RW-MH (same as before)
# =============================================================================

suppressPackageStartupMessages({
  library(BayesLogit)
  library(Matrix)
})

set.seed(7)

# -----------------------------------------------------------------------------
# Matern correlation (nu = 1 fixed throughout BARTSIMP)
# -----------------------------------------------------------------------------
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1
  out
}

# -----------------------------------------------------------------------------
# Simulate data
# -----------------------------------------------------------------------------
n_cl          <- 200
nu            <- 1
sigma_m_true  <- 1.0
rho_true      <- 0.3
kappa_true    <- sqrt(8 * nu) / rho_true

# Covariates: intercept + 2 standardized predictors
p             <- 3
beta_true     <- c(-0.3, 0.8, -0.6)         # intercept, x1, x2
locs          <- cbind(runif(n_cl), runif(n_cl))
X             <- cbind(1, matrix(rnorm(n_cl * (p - 1)), n_cl, p - 1))
colnames(X)   <- c("(Intercept)", paste0("x", seq_len(p - 1)))

n_i  <- as.numeric(sample(5:15, n_cl, replace = TRUE))

D <- as.matrix(dist(locs))
Sigma_z_true <- sigma_m_true^2 * matern_corr(D, kappa_true, nu)
z_true <- as.vector(t(chol(Sigma_z_true + 1e-6 * diag(n_cl))) %*% rnorm(n_cl))

eta_true <- as.vector(X %*% beta_true + z_true)
p_true   <- plogis(eta_true)
y_list   <- mapply(function(ni, pi) rbinom(ni, 1, pi),
                   n_i, p_true, SIMPLIFY = FALSE)
k_i      <- vapply(y_list, sum, numeric(1))

cat("Simulated:", n_cl, "clusters, total N =", sum(n_i),
    ", overall y mean =", round(mean(unlist(y_list)), 3), "\n")
cat("Covariate signal (sd of X beta):", round(sd(X %*% beta_true), 3),
    "  spatial signal (sd of z):", round(sd(z_true), 3), "\n\n")

# -----------------------------------------------------------------------------
# PC prior on (rho, sigma_m)
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
# Spatial helpers
# -----------------------------------------------------------------------------
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
n_iter <- 3000
burn   <- 1000

# Initial values
beta    <- rep(0, p)
z       <- rep(0, n_cl)
sigma_m <- 0.8
rho     <- 0.5
V0_inv  <- diag(p) / 100          # diffuse prior on beta

beta_draws    <- matrix(NA_real_, n_iter, p, dimnames = list(NULL, colnames(X)))
sigma_m_draws <- numeric(n_iter)
rho_draws     <- numeric(n_iter)
z_draws       <- matrix(NA_real_, n_iter, n_cl)
accept_psi    <- logical(n_iter)

sd_log_sigma <- 0.15
sd_log_rho   <- 0.15

t0 <- Sys.time()
for (it in seq_len(n_iter)) {

  # ---- 1. PG augmentation ----
  eta   <- as.vector(X %*% beta + z)
  Omega <- rpg(num = n_cl, h = n_i, z = eta)
  ybar  <- (k_i - n_i / 2) / Omega

  # ---- 2. beta | rest  (multivariate Gaussian) ----
  # X' diag(Omega) computed via row-broadcast: (X * Omega) has row i scaled by Omega_i
  XtW   <- t(X * Omega)                        # p x n
  Q_b   <- V0_inv + XtW %*% X
  rhs_b <- XtW %*% (ybar - z)
  U_b   <- chol(Q_b)
  m_b   <- backsolve(U_b, forwardsolve(t(U_b), rhs_b))
  beta  <- as.vector(m_b + backsolve(U_b, rnorm(p)))

  # ---- 3. z | rest  (dense Gaussian) ----
  S      <- make_Sigma_z(sigma_m, rho)
  Q_psi  <- chol2inv(chol(S))
  Q_z    <- Q_psi + diag(Omega)
  b_vec  <- Omega * (ybar - as.vector(X %*% beta))
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

  beta_draws[it, ]  <- beta
  sigma_m_draws[it] <- sigma_m
  rho_draws[it]     <- rho
  z_draws[it, ]     <- z

  if (it %% 500 == 0) {
    cat(sprintf("iter %5d   beta=[%s]  sigma_m=%5.3f  rho=%5.3f  acc=%.2f\n",
                it,
                paste(sprintf("% .2f", beta), collapse = ","),
                sigma_m, rho, mean(accept_psi[seq_len(it)])))
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
recovery_tbl <- rbind(
  do.call(rbind, lapply(seq_len(p), function(j) summ(beta_draws[keep, j], beta_true[j]))),
  sigma_m = summ(sigma_m_draws[keep], sigma_m_true),
  rho     = summ(rho_draws[keep],     rho_true)
)
rownames(recovery_tbl)[seq_len(p)] <- colnames(X)
print(round(recovery_tbl, 3))

z_post <- colMeans(z_draws[keep, ])
cat("\nz: cor(true, post.mean) =", round(cor(z_true, z_post), 3), "\n")
cat("   RMSE(z) =", round(sqrt(mean((z_true - z_post)^2)), 3), "\n")
cat("psi MH acceptance rate:", round(mean(accept_psi), 3), "\n")

# Coverage check across all parameters
covers <- recovery_tbl[, "q025"] <= recovery_tbl[, "truth"] &
          recovery_tbl[, "truth"] <= recovery_tbl[, "q975"]
cat("\n95% CIs covering truth:", sum(covers), "/", length(covers), "\n")

# -----------------------------------------------------------------------------
# Trace + recovery plot
# -----------------------------------------------------------------------------
pdf(file.path("dev", "02_pg_matern_linear_diagnostics.pdf"), width = 9, height = 8)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (j in seq_len(p)) {
  plot(beta_draws[, j], type = "l", ylab = colnames(X)[j], xlab = "iter",
       main = paste(colnames(X)[j], "trace"))
  abline(h = beta_true[j], col = "red", lty = 2)
}
plot(sigma_m_draws, type = "l", ylab = "sigma_m", xlab = "iter",
     main = "sigma_m trace"); abline(h = sigma_m_true, col = "red", lty = 2)
plot(rho_draws, type = "l", ylab = "rho", xlab = "iter",
     main = "rho trace"); abline(h = rho_true, col = "red", lty = 2)
plot(z_true, z_post, xlab = "z true", ylab = "z post mean",
     main = "z recovery", pch = 20, col = rgb(0,0,0,0.5))
abline(0, 1, col = "red", lty = 2)

par(op); dev.off()
cat("Diagnostics saved to dev/02_pg_matern_linear_diagnostics.pdf\n")
