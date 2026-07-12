# =============================================================================
# Validation V2 — multi-seed coverage
# =============================================================================
# For each scenario, run N_REPS replicates with different seeds and measure
# coverage of 95% credible intervals for (sigma_m, rho), and average BART/z
# RMSE across replicates.
#
# SMOKE_TEST mode runs a tiny version to verify the code path end-to-end
# without waiting for full results. Flip SMOKE_TEST <- FALSE to run on a
# cluster.
# =============================================================================

source("dev/bartsimp_pg.R")

SMOKE_TEST <- TRUE   # <-- flip to FALSE for production runs

if (SMOKE_TEST) {
  cat("*** SMOKE TEST MODE ***\n")
  N_CL     <- 100   # small dataset
  N_ITER   <- 150   # short chain
  BURN     <- 50
  N_REPS   <- 2     # bare minimum to test the aggregation path
  SCENARIOS <- c("balanced")   # one scenario only
  M_TREES  <- 5
  MESH_EDGE <- c(0.15, 0.40)
} else {
  N_CL     <- 200
  N_ITER   <- 3000
  BURN     <- 1000
  N_REPS   <- 10
  SCENARIOS <- c("cov_only", "balanced", "spatial_only")
  M_TREES  <- 20
  MESH_EDGE <- c(0.08, 0.25)
}

cat(sprintf("Config: n_cl=%d, n_iter=%d, burn=%d, n_reps=%d, scenarios=[%s]\n",
            N_CL, N_ITER, BURN, N_REPS, paste(SCENARIOS, collapse = ",")))

# -----------------------------------------------------------------------------
# Run one replicate within one scenario
# -----------------------------------------------------------------------------
run_one <- function(scenario, seed) {
  cov_sig <- scenario %in% c("cov_only", "balanced")
  spa_sig <- scenario %in% c("spatial_only", "balanced")
  dat <- gen_data(seed, n_cl = N_CL,
                  covariate_signal = cov_sig,
                  spatial_signal   = spa_sig)
  fit <- run_sampler_bartsimp_pg(
    X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
    f_true = dat$f_true, z_true = dat$z_true,
    n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
    mesh_max_edge = MESH_EDGE,
    verbose_every = 0,
    label = sprintf("%s/seed=%d", scenario, seed)
  )
  summ <- summarize_fit(fit,
                        f_true = dat$f_true, z_true = dat$z_true,
                        eta_true = dat$eta_true,
                        sigma_m_true = dat$sigma_m_sim,
                        rho_true     = dat$rho_sim)
  c(scenario   = scenario, seed = seed,
    sigma_m_truth = dat$sigma_m_sim,
    rho_truth     = ifelse(is.na(dat$rho_sim), NA_real_, dat$rho_sim),
    sigma_m_mean  = summ$sigma_m_mean,
    sigma_m_q025  = unname(summ$sigma_m_q025),
    sigma_m_q975  = unname(summ$sigma_m_q975),
    sigma_m_covered = as.numeric(summ$sigma_m_covered),
    rho_mean     = summ$rho_mean,
    rho_q025     = unname(summ$rho_q025),
    rho_q975     = unname(summ$rho_q975),
    rho_covered  = if (is.na(summ$rho_covered)) NA_real_ else as.numeric(summ$rho_covered),
    f_rmse       = if (is.null(summ$f_rmse))   NA_real_ else summ$f_rmse,
    z_rmse       = if (is.null(summ$z_rmse))   NA_real_ else summ$z_rmse,
    eta_cor      = if (is.null(summ$eta_cor))  NA_real_ else summ$eta_cor,
    psi_acc      = summ$psi_acc,
    elapsed_sec  = summ$elapsed_sec)
}

# -----------------------------------------------------------------------------
# Grid: scenario x seed
# -----------------------------------------------------------------------------
grid <- expand.grid(scenario = SCENARIOS,
                    seed     = seq_len(N_REPS) + 1000,
                    stringsAsFactors = FALSE)
cat(sprintf("Total runs: %d\n\n", nrow(grid)))

t_total <- Sys.time()
rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  cat(sprintf("[%d/%d] scenario=%s seed=%d ...\n",
              i, nrow(grid), grid$scenario[i], grid$seed[i]))
  rows[[i]] <- run_one(grid$scenario[i], grid$seed[i])
}
t_total <- difftime(Sys.time(), t_total, units = "secs")

results <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
# coerce numerics
num_cols <- setdiff(names(results), c("scenario"))
results[num_cols] <- lapply(results[num_cols], as.numeric)

cat(sprintf("\n=== Done in %s ===\n", format(t_total)))
cat("\nPer-run summary:\n")
print(results, row.names = FALSE)

# Aggregate by scenario
cat("\n=== Aggregate by scenario ===\n")
agg <- aggregate(cbind(sigma_m_covered, rho_covered,
                       f_rmse, z_rmse, eta_cor, psi_acc, elapsed_sec)
                 ~ scenario,
                 data = results, FUN = function(x) mean(x, na.rm = TRUE))
print(agg, row.names = FALSE)

# Save
saveRDS(list(results = results, aggregate = agg,
             config = list(N_CL = N_CL, N_ITER = N_ITER, BURN = BURN,
                           N_REPS = N_REPS, SCENARIOS = SCENARIOS,
                           SMOKE_TEST = SMOKE_TEST)),
        sprintf("dev/06_validation_v2_%s.rds",
                if (SMOKE_TEST) "smoke" else "full"))
cat(sprintf("\nResults saved to dev/06_validation_v2_%s.rds\n",
            if (SMOKE_TEST) "smoke" else "full"))
