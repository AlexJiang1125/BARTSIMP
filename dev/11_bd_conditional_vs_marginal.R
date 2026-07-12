# =============================================================================
# Local smoke test: BD conditional (Path B) vs BD marginal (Path A)
# =============================================================================
# Same simulated data, same priors, same n_iter. Compare:
#   - Posterior summaries (sigma_m, rho, BART surface, spatial field)
#   - psi MH acceptance rate
#   - BD acceptance rate
#   - cor(f, z) and variance-partition diagnostics
#   - Wall time
#
# Purpose: verify Path A runs without errors, gives qualitatively similar
# posteriors to Path B on a simple balanced-scenario dataset, and gauge the
# computational cost ratio.
# =============================================================================

source("dev/bartsimp_pg.R")

set.seed(2026)
N_CL    <- 100    # small for speed
N_ITER  <- 300    # smoke-test length
BURN    <- 100
M_TREES <- 10

cat("Simulating balanced-scenario binary dataset (n=", N_CL, ")...\n", sep = "")
dat <- gen_data(seed = 2026, n_cl = N_CL,
                covariate_signal = TRUE, spatial_signal = TRUE)

# -----------------------------------------------------------------------------
# Path B: BD conditional
# -----------------------------------------------------------------------------
cat("\n--- Path B: BD conditional ---\n")
t0 <- Sys.time()
fit_B <- run_sampler_bartsimp_pg(
  X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
  f_true = dat$f_true, z_true = dat$z_true,
  n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
  mesh_max_edge = c(0.15, 0.40),
  bd_mode = "conditional",
  verbose_every = 100, label = "Path-B"
)
cat(sprintf("Path B elapsed: %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# -----------------------------------------------------------------------------
# Path A: BD marginal
# -----------------------------------------------------------------------------
cat("\n--- Path A: BD marginal ---\n")
t0 <- Sys.time()
fit_A <- run_sampler_bartsimp_pg(
  X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
  f_true = dat$f_true, z_true = dat$z_true,
  n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
  mesh_max_edge = c(0.15, 0.40),
  bd_mode = "marginal",
  verbose_every = 100, label = "Path-A"
)
cat(sprintf("Path A elapsed: %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# -----------------------------------------------------------------------------
# Compare
# -----------------------------------------------------------------------------
summ <- function(fit, dat) {
  keep <- (fit$burn + 1):fit$n_iter
  f_p  <- colMeans(fit$f_draws[keep, , drop = FALSE])
  z_p  <- colMeans(fit$z_draws[keep, , drop = FALSE])
  list(
    sigma_m_mean = mean(fit$sigma_m_draws[keep]),
    sigma_m_ci   = quantile(fit$sigma_m_draws[keep], c(0.025, 0.975)),
    rho_mean     = mean(fit$rho_draws[keep]),
    rho_ci       = quantile(fit$rho_draws[keep], c(0.025, 0.975)),
    f_rmse       = sqrt(mean((f_p - dat$f_true)^2)),
    z_rmse       = sqrt(mean((z_p - dat$z_true)^2)),
    eta_cor      = cor(dat$eta_true, f_p + z_p),
    fz_cor_mean  = mean(fit$fz_cor_draws[keep]),
    var_f_mean   = mean(fit$varpart_draws[keep, "var_f"]),
    var_z_mean   = mean(fit$varpart_draws[keep, "var_z"]),
    psi_acc      = mean(fit$accept_psi),
    elapsed_sec  = fit$elapsed_sec
  )
}

sB <- summ(fit_B, dat); sA <- summ(fit_A, dat)
cat("\n=== Side-by-side recovery ===\n")
cat(sprintf("%-22s %-22s %-22s\n", "metric", "Path B (conditional)", "Path A (marginal)"))
fmt_ci <- function(m, ci) sprintf("%.3f [%.3f, %.3f]", m, ci[1], ci[2])
cat(sprintf("%-22s %-22s %-22s\n", "sigma_m",
            fmt_ci(sB$sigma_m_mean, sB$sigma_m_ci),
            fmt_ci(sA$sigma_m_mean, sA$sigma_m_ci)))
cat(sprintf("%-22s %-22s %-22s\n", "rho",
            fmt_ci(sB$rho_mean, sB$rho_ci),
            fmt_ci(sA$rho_mean, sA$rho_ci)))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "f RMSE",       sB$f_rmse,       sA$f_rmse))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "z RMSE",       sB$z_rmse,       sA$z_rmse))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "eta cor",      sB$eta_cor,      sA$eta_cor))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "cor(f, z)",    sB$fz_cor_mean,  sA$fz_cor_mean))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "var(f) mean",  sB$var_f_mean,   sA$var_f_mean))
cat(sprintf("%-22s %-22.3f %-22.3f\n", "var(z) mean",  sB$var_z_mean,   sA$var_z_mean))
cat(sprintf("%-22s %-22.2f %-22.2f\n", "psi MH acc",   sB$psi_acc,      sA$psi_acc))
cat(sprintf("%-22s %-22.1f %-22.1f\n", "elapsed (s)",  sB$elapsed_sec,  sA$elapsed_sec))
cat(sprintf("%-22s %-22s %-22.2fx\n",  "speed ratio",  "1x", sA$elapsed_sec / sB$elapsed_sec))

cat("\n=== Smoke test PASSED ===\n")
