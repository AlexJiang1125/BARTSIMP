# =============================================================================
# RSR (linear+interaction basis) vs original BARTSIMP-PG
# =============================================================================
# Tests whether hard-orthogonalizing z against a fixed low-rank basis
#   B0 = [1, x1, x2, x1*x2]
# reduces the (f, z) confounding observed in V1 scenario C (spatial-only),
# without harming recovery in scenarios A (cov-only) and B (balanced).
#
# Metrics reported for each scenario x method:
#   sigma_m posterior (mean + 95% CI)
#   BART surface RMSE vs truth, max|f_post|
#   z RMSE
#   var(f_post), var(z_post)
#   cor(f, z) per-iter mean
#   psi MH acceptance
# =============================================================================

source("dev/bartsimp_pg.R")

SMOKE_TEST <- TRUE   # flip to FALSE for production

if (SMOKE_TEST) {
  N_CL    <- 200
  N_ITER  <- 1500
  BURN    <- 500
  M_TREES <- 20
  SCENARIOS <- c("cov_only", "balanced", "spatial_only")
  RSR_MODES <- c("none", "linear_interaction")
} else {
  N_CL    <- 200
  N_ITER  <- 3000
  BURN    <- 1000
  M_TREES <- 20
  SCENARIOS <- c("cov_only", "balanced", "spatial_only")
  RSR_MODES <- c("none", "linear_interaction")
}

SEED      <- 2024
MESH_EDGE <- c(0.08, 0.25)

summarize <- function(fit, dat) {
  keep <- (fit$burn + 1):fit$n_iter
  f_p  <- colMeans(fit$f_draws[keep, , drop = FALSE])
  z_p  <- colMeans(fit$z_draws[keep, , drop = FALSE])
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
    elapsed_sec    = fit$elapsed_sec
  )
}

results <- list()
for (scen in SCENARIOS) {
  cov_sig <- scen %in% c("cov_only", "balanced")
  spa_sig <- scen %in% c("spatial_only", "balanced")
  dat <- gen_data(SEED, n_cl = N_CL,
                  covariate_signal = cov_sig, spatial_signal = spa_sig)

  cat("\n========================\n")
  cat(sprintf("Scenario %s: cov_sig=%s, spa_sig=%s\n", scen, cov_sig, spa_sig))

  for (rsr in RSR_MODES) {
    tag <- sprintf("%s/%s", scen, rsr)
    cat(sprintf("\n--- %s ---\n", tag))
    fit <- run_sampler_bartsimp_pg(
      X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
      f_true = dat$f_true, z_true = dat$z_true,
      n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
      mesh_max_edge = MESH_EDGE,
      rsr_basis = if (rsr == "none") NULL else rsr,
      verbose_every = 300, label = tag
    )
    results[[tag]] <- list(fit = fit, dat = dat, summary = summarize(fit, dat))
  }
}

# Summary table
cat("\n\n============ RSR vs ORIGINAL ============\n")
rows <- lapply(results, function(r) {
  s <- r$summary
  data.frame(
    label       = s$label,
    sigma_m     = sprintf("%.2f [%.2f,%.2f]", s$sigma_m_mean, s$sigma_m_q025, s$sigma_m_q975),
    f_rmse      = sprintf("%.3f", s$f_rmse),
    f_max_abs   = sprintf("%.2f", s$f_max_abs),
    z_rmse      = sprintf("%.3f", s$z_rmse),
    eta_cor     = sprintf("%.3f", s$eta_cor),
    var_f       = sprintf("%.2f", s$var_f_mean),
    var_z       = sprintf("%.2f", s$var_z_mean),
    fz_cor      = sprintf("%.3f", s$fz_cor_mean)
  )
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

# Side-by-side diagnostic plot
pdf("dev/12_rsr_vs_original_diagnostics.pdf",
    width = 10, height = 4 * length(SCENARIOS))
op <- par(mfrow = c(length(SCENARIOS), 3), mar = c(4, 4, 2, 1))
for (scen in SCENARIOS) {
  tag_none <- sprintf("%s/none", scen)
  tag_rsr  <- sprintf("%s/linear_interaction", scen)
  for (tag in c(tag_none, tag_rsr)) {
    r <- results[[tag]]
    keep <- (r$fit$burn + 1):r$fit$n_iter
    f_p  <- colMeans(r$fit$f_draws[keep, , drop = FALSE])
    plot(r$dat$f_true, f_p, pch = 20, col = rgb(0,0,0,0.4),
         xlab = "f true", ylab = "f post mean",
         main = sprintf("%s BART", tag),
         ylim = c(-3, 3))
    abline(0, 1, col = "red", lty = 2)
  }
  # var(f) trace overlay
  f_none <- results[[tag_none]]$fit$varpart_draws[, "var_f"]
  f_rsr  <- results[[tag_rsr]]$fit$varpart_draws[, "var_f"]
  ylim <- range(c(f_none, f_rsr), na.rm = TRUE)
  plot(f_none, type = "l", col = "firebrick", lwd = 0.8,
       xlab = "iter", ylab = "var(f)", main = sprintf("%s var(f)", scen),
       ylim = ylim)
  lines(f_rsr, col = "steelblue", lwd = 0.8)
  legend("topright", legend = c("original", "RSR"),
         col = c("firebrick", "steelblue"), lwd = 1.5, bty = "n")
}
par(op); dev.off()

saveRDS(results, "dev/12_rsr_vs_original.rds")
cat("\nSaved dev/12_rsr_vs_original.rds and .pdf\n")
