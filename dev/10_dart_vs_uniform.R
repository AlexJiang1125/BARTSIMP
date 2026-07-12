# =============================================================================
# Validation V6 -- DART sparsity prior vs uniform variable selection
# =============================================================================
# Re-runs the V1 decomposition scenarios with DART enabled. The (f, z)
# identifiability issue observed in scenario C (spatial-only) is the main
# target: DART limits BART's tendency to split on spurious covariates, which
# should reduce the absorbed-spatial-signal artifact.
#
# We also compute identifiability diagnostics at every iteration:
#   - cor(f, z): negative correlation indicates trade-off between components
#   - var(f), var(z), cov(f, z): variance partition of eta = f + z
# =============================================================================

source("dev/bartsimp_pg.R")

SMOKE_TEST <- TRUE   # flip to FALSE for production

if (SMOKE_TEST) {
  N_CL    <- 200
  N_ITER  <- 1500
  BURN    <- 500
  M_TREES <- 20
  SCENARIOS <- c("balanced", "spatial_only")  # focus on the ones at issue
} else {
  N_CL    <- 200
  N_ITER  <- 3000
  BURN    <- 1000
  M_TREES <- 20
  SCENARIOS <- c("cov_only", "balanced", "spatial_only")
}

ALPHA_DART <- 0.5   # Linero default-ish; lower => more sparsity

# Replicate seed used in V1
SEED <- 2024
MESH_EDGE <- c(0.08, 0.25)

summarize <- function(fit, dat) {
  keep <- (fit$burn + 1):fit$n_iter
  f_p <- colMeans(fit$f_draws[keep, , drop = FALSE])
  z_p <- colMeans(fit$z_draws[keep, , drop = FALSE])
  eta_true <- dat$f_true + dat$z_true
  vp <- fit$varpart_draws[keep, , drop = FALSE]
  list(
    label          = fit$label,
    sigma_m_mean   = mean(fit$sigma_m_draws[keep]),
    sigma_m_q025   = quantile(fit$sigma_m_draws[keep], 0.025),
    sigma_m_q975   = quantile(fit$sigma_m_draws[keep], 0.975),
    f_rmse         = sqrt(mean((f_p - dat$f_true)^2)),
    f_max_abs      = max(abs(f_p)),
    z_rmse         = sqrt(mean((z_p - dat$z_true)^2)),
    eta_cor        = cor(eta_true, f_p + z_p),
    fz_cor_mean    = mean(fit$fz_cor_draws[keep]),
    var_f_mean     = mean(vp[, "var_f"]),
    var_z_mean     = mean(vp[, "var_z"]),
    cov_fz_mean    = mean(vp[, "cov_fz"]),
    psi_acc        = mean(fit$accept_psi),
    var_probs_mean = if (!is.null(fit$var_probs_draws))
        colMeans(fit$var_probs_draws[keep, , drop = FALSE]) else NULL
  )
}

# -----------------------------------------------------------------------------
# Run pairs
# -----------------------------------------------------------------------------
results <- list()
for (scen in SCENARIOS) {
  cov_sig <- scen %in% c("cov_only", "balanced")
  spa_sig <- scen %in% c("spatial_only", "balanced")
  dat <- gen_data(SEED, n_cl = N_CL,
                  covariate_signal = cov_sig, spatial_signal = spa_sig)

  cat("\n========================\n")
  cat(sprintf("Scenario %s: cov_sig=%s, spa_sig=%s\n", scen, cov_sig, spa_sig))

  for (use_dart in c(FALSE, TRUE)) {
    tag <- sprintf("%s-%s", scen, if (use_dart) "DART" else "unif")
    cat(sprintf("\n--- %s ---\n", tag))
    fit <- run_sampler_bartsimp_pg(
      X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
      f_true = dat$f_true, z_true = dat$z_true,
      n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
      mesh_max_edge = MESH_EDGE,
      use_dart = use_dart, alpha_dart = ALPHA_DART,
      verbose_every = 300, label = tag
    )
    results[[tag]] <- list(fit = fit, dat = dat, summary = summarize(fit, dat))
  }
}

# -----------------------------------------------------------------------------
# Print summary table
# -----------------------------------------------------------------------------
cat("\n\n============ SUMMARY TABLE ============\n")
rows <- lapply(results, function(r) {
  s <- r$summary
  data.frame(
    label       = s$label,
    sigma_m     = sprintf("%.2f [%.2f,%.2f]", s$sigma_m_mean, s$sigma_m_q025, s$sigma_m_q975),
    f_rmse      = sprintf("%.3f", s$f_rmse),
    f_max_abs   = sprintf("%.2f", s$f_max_abs),
    z_rmse      = sprintf("%.3f", s$z_rmse),
    eta_cor     = sprintf("%.3f", s$eta_cor),
    fz_cor      = sprintf("%.3f", s$fz_cor_mean),
    var_f_z_cov = sprintf("(%.2f, %.2f, %.2f)", s$var_f_mean, s$var_z_mean, s$cov_fz_mean)
  )
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

if (any(vapply(results, function(r) !is.null(r$summary$var_probs_mean), logical(1)))) {
  cat("\nVariable-inclusion probabilities (DART runs only):\n")
  for (tag in names(results)) {
    vp <- results[[tag]]$summary$var_probs_mean
    if (!is.null(vp))
      cat(sprintf("  %s: %s\n", tag, paste(sprintf("%.3f", vp), collapse=" ")))
  }
}

# -----------------------------------------------------------------------------
# Diagnostic plot
# -----------------------------------------------------------------------------
pdf("dev/10_dart_vs_uniform_diagnostics.pdf", width = 9, height = 4 * length(SCENARIOS))
op <- par(mfrow = c(length(SCENARIOS), 3), mar = c(4, 4, 2, 1))
for (scen in SCENARIOS) {
  for (use_dart in c(FALSE, TRUE)) {
    tag <- sprintf("%s-%s", scen, if (use_dart) "DART" else "unif")
    r <- results[[tag]]
    keep <- (r$fit$burn + 1):r$fit$n_iter
    f_p  <- colMeans(r$fit$f_draws[keep, , drop = FALSE])
    plot(r$dat$f_true, f_p, pch = 20, col = rgb(0,0,0,0.4),
         xlab = "f true", ylab = "f post mean",
         main = sprintf("%s BART", tag))
    abline(0, 1, col = "red", lty = 2)
  }
  # third panel: fz_cor trace overlaid
  unif_tag <- sprintf("%s-unif", scen)
  dart_tag <- sprintf("%s-DART", scen)
  fz_unif <- results[[unif_tag]]$fit$fz_cor_draws
  fz_dart <- results[[dart_tag]]$fit$fz_cor_draws
  ylim <- range(c(fz_unif, fz_dart), na.rm = TRUE)
  plot(fz_unif, type = "l", col = "firebrick", lwd = 1.2,
       xlab = "iter", ylab = "cor(f, z)",
       main = sprintf("%s cor(f,z)", scen), ylim = ylim)
  lines(fz_dart, col = "steelblue", lwd = 1.2)
  abline(h = 0, lty = 2, col = "gray")
  legend("topright", legend = c("uniform", "DART"),
         col = c("firebrick", "steelblue"), lwd = 1.2, bty = "n")
}
par(op); dev.off()

# Save
saveRDS(results, "dev/10_dart_vs_uniform.rds")
cat(sprintf("\nResults saved to dev/10_dart_vs_uniform.rds and .pdf\n"))
