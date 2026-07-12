# =============================================================================
# Validation V1 — signal decomposition
# =============================================================================
# Three runs of the Woodbury PG+BART sampler, varying the mixing of signal
# between BART and spatial:
#   A) covariate-only:  f_true = tree-structured,  z_true = 0 (sigma_m_sim = 0)
#   B) balanced:        both signals present                  (Proto 2b setup)
#   C) spatial-only:    f_true = 0,                z_true = GRF
#
# We check that the sampler correctly attributes signal:
#   - In A, posterior sigma_m should be small (no spatial signal)
#   - In C, the BART surface should be ~0 (no covariate signal)
#   - In B, both components should be recovered
# =============================================================================

suppressPackageStartupMessages({
  library(BayesLogit)
  library(Matrix)
  library(INLA)
})
INLA::inla.setOption(fmesher.evolution.warn = FALSE)
source("dev/rbart.R")

# -----------------------------------------------------------------------------
# Helpers (copied from Proto 2b)
# -----------------------------------------------------------------------------
matern_corr <- function(d, kappa, nu = 1) {
  out <- (2^(1 - nu)) / gamma(nu) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1; out
}

make_Qpsi <- function(spde, sigma_m, rho) {
  INLA::inla.spde.precision(spde, theta = c(log(rho), log(sigma_m)))
}

margprob_woodbury <- function(yhat, Omega, A, Q_psi) {
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

draw_z_woodbury <- function(m_u, Fchol, A) {
  eta <- rnorm(length(m_u))
  u   <- m_u + as.numeric(Matrix::solve(Fchol, eta, system = "Lt"))
  as.numeric(A %*% u)
}

# Single MCMC run, returns posterior summaries
run_sampler <- function(X, locs, n_i, k_i, f_true, z_true,
                        n_iter = 1200, burn = 400, label = "") {
  n_cl <- nrow(X)

  # Mesh and SPDE
  mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(0.08, 0.25), cutoff = 0.03)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.1, 0.5),
                                          prior.sigma = c(1.0, 0.5))
  A <- as(INLA::inla.spde.make.A(mesh = mesh, loc = locs), "CsparseMatrix")

  # Priors
  m_trees <- 20
  tau2    <- (3 / (2 * sqrt(m_trees)))^2
  alpha_bart <- 0.95; beta_bart <- 2

  d_dim <- 2; rho_0 <- 0.10; sigma_0 <- 2.0; alpha_1 <- 0.05; alpha_2 <- 0.05
  lam1 <- -log(alpha_1) * rho_0^(d_dim/2)
  lam2 <- -log(alpha_2) / sigma_0
  log_pc_prior <- function(rho, sigma_m) {
    log(d_dim/2) + log(lam1) + log(lam2) -
      (d_dim/2 + 1) * log(rho) - lam1 * rho^(-d_dim/2) - lam2 * sigma_m
  }

  # State
  ensemble  <- ensemble_init(m_trees, init_val = 0)
  z         <- rep(0, n_cl)
  sigma_m   <- 0.8
  rho       <- 0.5
  f_current <- rep(0, n_cl)

  # Storage
  sigma_m_draws <- numeric(n_iter)
  rho_draws     <- numeric(n_iter)
  f_draws       <- matrix(NA_real_, n_iter, n_cl)
  z_draws       <- matrix(NA_real_, n_iter, n_cl)
  accept_psi    <- logical(n_iter)

  sd_log <- 0.15
  t0 <- Sys.time()
  for (it in seq_len(n_iter)) {
    eta   <- f_current + z
    Omega <- rpg(num = n_cl, h = n_i, z = eta)
    ybar  <- (k_i - n_i / 2) / Omega

    sweep <- ensemble_sweep(ensemble, X, ybar - z, Omega, tau2,
                            alpha = alpha_bart, beta = beta_bart)
    ensemble  <- sweep$ensemble
    f_current <- sweep$fit

    yhat   <- ybar - f_current
    Q_psi  <- make_Qpsi(spde, sigma_m, rho)
    curr   <- margprob_woodbury(yhat, Omega, A, Q_psi)

    sigma_m_p <- exp(log(sigma_m) + rnorm(1, 0, sd_log))
    rho_p     <- exp(log(rho)     + rnorm(1, 0, sd_log))
    Q_psi_p   <- make_Qpsi(spde, sigma_m_p, rho_p)
    prop      <- margprob_woodbury(yhat, Omega, A, Q_psi_p)
    log_alpha <- (prop$mlik + log_pc_prior(rho_p, sigma_m_p)) -
                 (curr$mlik + log_pc_prior(rho,   sigma_m))   +
                 (log(sigma_m_p) + log(rho_p)) -
                 (log(sigma_m)   + log(rho))
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      sigma_m <- sigma_m_p; rho <- rho_p
      accepted <- prop; accept_psi[it] <- TRUE
    } else { accepted <- curr }

    z <- draw_z_woodbury(accepted$m_u, accepted$Fchol, A)

    sigma_m_draws[it] <- sigma_m
    rho_draws[it]     <- rho
    f_draws[it, ]     <- f_current
    z_draws[it, ]     <- z

    if (it %% 200 == 0) {
      cat(sprintf("  [%s] iter %4d  sigma_m=%.3f  rho=%.3f  RMSE(f)=%.3f  RMSE(z)=%.3f  acc=%.2f\n",
                  label, it, sigma_m, rho,
                  sqrt(mean((f_current - f_true)^2)),
                  sqrt(mean((z       - z_true)^2)),
                  mean(accept_psi[seq_len(it)])))
    }
  }
  t1 <- Sys.time()
  cat(sprintf("  [%s] elapsed: %s\n", label, format(t1 - t0)))

  keep <- (burn + 1):n_iter
  list(
    sigma_m_post = sigma_m_draws[keep],
    rho_post     = rho_draws[keep],
    f_post       = colMeans(f_draws[keep, ]),
    z_post       = colMeans(z_draws[keep, ]),
    f_true       = f_true,
    z_true       = z_true,
    accept_psi   = mean(accept_psi),
    label        = label
  )
}

# -----------------------------------------------------------------------------
# Shared simulation skeleton: locations, X, n_i, k_i  -- vary f_true, z_true
# -----------------------------------------------------------------------------
gen_data <- function(seed, n_cl = 200,
                     covariate_signal = TRUE,
                     spatial_signal   = TRUE,
                     sigma_m_sim = 0.8, rho_sim = 0.3) {
  set.seed(seed)
  nu <- 1
  x1 <- runif(n_cl); x2 <- runif(n_cl)
  X  <- cbind(x1, x2); colnames(X) <- c("x1", "x2")
  locs <- cbind(runif(n_cl), runif(n_cl))
  n_i  <- as.numeric(sample(5:15, n_cl, replace = TRUE))

  if (covariate_signal) {
    f_true <- ifelse(x1 < 0.5, 1.5, ifelse(x2 < 0.5, -1.0, 0))
  } else {
    f_true <- rep(0, n_cl)
  }
  if (spatial_signal) {
    D <- as.matrix(dist(locs))
    kappa_sim <- sqrt(8 * nu) / rho_sim
    Sz <- sigma_m_sim^2 * matern_corr(D, kappa_sim, nu) + 1e-6 * diag(n_cl)
    z_true <- as.vector(t(chol(Sz)) %*% rnorm(n_cl))
  } else {
    z_true <- rep(0, n_cl)
  }
  eta_true <- f_true + z_true
  p_true   <- plogis(eta_true)
  k_i      <- vapply(seq_len(n_cl),
                     function(i) sum(rbinom(n_i[i], 1, p_true[i])),
                     numeric(1))
  list(X = X, locs = locs, n_i = n_i, k_i = k_i,
       f_true = f_true, z_true = z_true,
       sigma_m_sim = if (spatial_signal) sigma_m_sim else 0,
       rho_sim = if (spatial_signal) rho_sim else NA)
}

# -----------------------------------------------------------------------------
# Run the three scenarios
# -----------------------------------------------------------------------------
SEED <- 2024
N_ITER <- 1200; BURN <- 400

cat("\n=== Scenario A: covariate-only (no spatial signal) ===\n")
dat_A <- gen_data(SEED, covariate_signal = TRUE, spatial_signal = FALSE)
res_A <- run_sampler(dat_A$X, dat_A$locs, dat_A$n_i, dat_A$k_i,
                     dat_A$f_true, dat_A$z_true,
                     n_iter = N_ITER, burn = BURN, label = "A:cov-only")

cat("\n=== Scenario B: balanced (covariate + spatial) ===\n")
dat_B <- gen_data(SEED, covariate_signal = TRUE, spatial_signal = TRUE)
res_B <- run_sampler(dat_B$X, dat_B$locs, dat_B$n_i, dat_B$k_i,
                     dat_B$f_true, dat_B$z_true,
                     n_iter = N_ITER, burn = BURN, label = "B:balanced")

cat("\n=== Scenario C: spatial-only (no covariate signal) ===\n")
dat_C <- gen_data(SEED, covariate_signal = FALSE, spatial_signal = TRUE)
res_C <- run_sampler(dat_C$X, dat_C$locs, dat_C$n_i, dat_C$k_i,
                     dat_C$f_true, dat_C$z_true,
                     n_iter = N_ITER, burn = BURN, label = "C:spatial-only")

# -----------------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------------
summarize_run <- function(res, dat) {
  c(scenario      = res$label,
    sigma_m_true  = sprintf("%.2f", dat$sigma_m_sim),
    sigma_m_mean  = sprintf("%.3f", mean(res$sigma_m_post)),
    sigma_m_95CI  = sprintf("[%.2f, %.2f]",
                            quantile(res$sigma_m_post, 0.025),
                            quantile(res$sigma_m_post, 0.975)),
    rho_mean      = sprintf("%.3f", mean(res$rho_post)),
    BART_RMSE     = sprintf("%.3f", sqrt(mean((res$f_post - res$f_true)^2))),
    BART_max_abs  = sprintf("%.3f", max(abs(res$f_post))),
    z_RMSE        = sprintf("%.3f", sqrt(mean((res$z_post - res$z_true)^2))),
    z_max_abs     = sprintf("%.3f", max(abs(res$z_post))),
    psi_acc       = sprintf("%.2f", res$accept_psi))
}

cat("\n=== Decomposition validation summary ===\n")
tbl <- rbind(summarize_run(res_A, dat_A),
             summarize_run(res_B, dat_B),
             summarize_run(res_C, dat_C))
print(tbl, quote = FALSE)

cat("\nKey checks:\n")
cat("- Scenario A (no spatial signal):\n")
cat(sprintf("    posterior sigma_m mean %.3f  (should be small, well below 0.8)\n",
            mean(res_A$sigma_m_post)))
cat(sprintf("    posterior z RMSE %.3f vs zero truth (should be small)\n",
            sqrt(mean(res_A$z_post^2))))
cat("- Scenario C (no covariate signal):\n")
cat(sprintf("    BART max |f_post| %.3f  (should be near 0 across covariate space)\n",
            max(abs(res_C$f_post))))
cat(sprintf("    BART RMSE vs zero truth %.3f\n",
            sqrt(mean(res_C$f_post^2))))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
pdf("dev/05_validation_decomposition.pdf", width = 10, height = 9)
op <- par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

plot_row <- function(res, label) {
  plot(res$sigma_m_post, type = "l", xlab = "iter (post-burn)", ylab = "sigma_m",
       main = paste(label, "sigma_m"))
  abline(h = if (grepl("spatial-only", label) || grepl("balanced", label)) 0.8 else 0,
         col = "red", lty = 2)
  plot(res$f_true, res$f_post, pch = 20, col = rgb(0,0,0,0.5),
       xlab = "f true", ylab = "f post mean", main = paste(label, "BART"))
  abline(0, 1, col = "red", lty = 2)
  plot(res$z_true, res$z_post, pch = 20, col = rgb(0,0,0,0.5),
       xlab = "z true", ylab = "z post mean", main = paste(label, "spatial"))
  abline(0, 1, col = "red", lty = 2)
}
plot_row(res_A, "A cov-only")
plot_row(res_B, "B balanced")
plot_row(res_C, "C spatial-only")
par(op); dev.off()
cat("\nDiagnostics saved to dev/05_validation_decomposition.pdf\n")

saveRDS(list(A = res_A, B = res_B, C = res_C), "dev/05_validation_decomposition.rds")
