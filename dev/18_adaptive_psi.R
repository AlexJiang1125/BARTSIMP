# =============================================================================
# Test adaptive psi MH (Robbins-Monro shared scale).
# Compares against the fixed-scale chain on the SAME dataset/seed.
#
# Pass criteria:
#   (a) running acceptance rate settles near psi_target_acc by end of burn-in
#   (b) sigma_m / rho posteriors do not shift substantially (both should be
#       valid samplers; adaptive should just mix better)
#   (c) ESS for sigma_m and rho is HIGHER under the adaptive sampler
#   (d) ACR / RMSE on the grid don't regress
# =============================================================================

source("dev/bartsimp_pg.R")
suppressPackageStartupMessages(library(coda))

set.seed(2026)

# ---- toy config -----------------------------------------------------------
N_CL     <- 150L
N_ITER   <- 1500L
BURN     <-  500L
M_TREES  <-   20L
GRID_DIM <-   20L

# Same balanced-scenario dataset generator as dev/13
sigma_m_baseline <- 0.5
rho_baseline     <- sqrt(8 * 1) / 2.5
surface_f0 <- function(x1, x2)
  ifelse(x1 < 0.5, 3, ifelse(x2 < 0.5, -2, 0))
OMEGA <- 0.5

G    <- GRID_DIM^2
step <- 1 / GRID_DIM
gx   <- (seq_len(GRID_DIM) - 0.5) * step
cell_centers <- expand.grid(lon = gx, lat = gx)
x1_grid <- runif(G); x2_grid <- runif(G)

Dg    <- as.matrix(dist(cell_centers))
kappa <- sqrt(8 * 1) / rho_baseline
Sg    <- sigma_m_baseline^2 * matern_corr(Dg, kappa, 1) + 1e-6 * diag(G)
z_grid <- as.vector(t(chol(Sg)) %*% rnorm(G))
f0_grid  <- surface_f0(x1_grid, x2_grid)
eta_grid <- (1 - OMEGA) * z_grid + OMEGA * f0_grid

cluster_cells <- sample.int(G, N_CL, replace = FALSE)
jitter_x <- runif(N_CL, -step/2, step/2)
jitter_y <- runif(N_CL, -step/2, step/2)
locs <- cbind(
  lon = cell_centers$lon[cluster_cells] + jitter_x,
  lat = cell_centers$lat[cluster_cells] + jitter_y
)
X   <- cbind(x1 = x1_grid[cluster_cells], x2 = x2_grid[cluster_cells])
eta_train <- eta_grid[cluster_cells]
p_train   <- plogis(eta_train)
n_i <- as.numeric(sample(5:15, N_CL, replace = TRUE))
k_i <- vapply(seq_len(N_CL),
              function(i) sum(rbinom(n_i[i], 1, p_train[i])),
              numeric(1))
f_train <- OMEGA * surface_f0(X[, 1], X[, 2])
z_train <- (1 - OMEGA) * z_grid[cluster_cells]
X_grid    <- cbind(x1 = x1_grid, x2 = x2_grid)
locs_grid <- as.matrix(cell_centers[, c("lon", "lat")])

cat(sprintf("Dataset: n_cl=%d, grid=%d cells, omega=%.2f\n",
            N_CL, G, OMEGA))

# ---- run FIXED-scale (adapt off) ------------------------------------------
cat("\n--- Fixed-scale (adapt_psi = FALSE) ---\n")
set.seed(2026)
t0 <- Sys.time()
fit_fixed <- run_sampler_bartsimp_pg(
  X = X, locs = locs, n_i = n_i, k_i = k_i,
  f_true = f_train, z_true = z_train,
  X_pred = X_grid, locs_pred = locs_grid,
  n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
  mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
  bd_mode = "conditional",
  adapt_psi = FALSE,
  verbose_every = 500, label = "fixed"
)
cat(sprintf("Fixed elapsed: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- run ADAPTIVE -----------------------------------------------------------
cat("\n--- Adaptive (adapt_psi = TRUE, target = 0.234) ---\n")
set.seed(2026)
t0 <- Sys.time()
fit_adapt <- run_sampler_bartsimp_pg(
  X = X, locs = locs, n_i = n_i, k_i = k_i,
  f_true = f_train, z_true = z_train,
  X_pred = X_grid, locs_pred = locs_grid,
  n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
  mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
  bd_mode = "conditional",
  adapt_psi = TRUE, psi_target_acc = 0.234, adapt_decay = 0.6,
  verbose_every = 500, label = "adapt"
)
cat(sprintf("Adaptive elapsed: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- diagnostics -----------------------------------------------------------
keep <- (fit_fixed$burn + 1):fit_fixed$n_iter
n_keep <- length(keep)

summ <- function(fit, lbl) {
  keep   <- (fit$burn + 1):fit$n_iter
  acc_post <- mean(fit$accept_psi[keep])
  acc_burn <- mean(fit$accept_psi[seq_len(fit$burn)])
  acc_last200 <- mean(fit$accept_psi[(fit$burn - 199):fit$burn])
  ess_sm <- coda::effectiveSize(coda::mcmc(fit$sigma_m_draws[keep]))
  ess_rh <- coda::effectiveSize(coda::mcmc(fit$rho_draws[keep]))
  cat(sprintf("\n  [%s]\n", lbl))
  cat(sprintf("    acc (burn full)            = %.3f\n", acc_burn))
  cat(sprintf("    acc (last 200 of burn)     = %.3f\n", acc_last200))
  cat(sprintf("    acc (post-burn)            = %.3f\n", acc_post))
  cat(sprintf("    sigma_m mean (truth 0.500) = %.3f\n",
              mean(fit$sigma_m_draws[keep])))
  cat(sprintf("    rho     mean (truth 1.131) = %.3f\n",
              mean(fit$rho_draws[keep])))
  cat(sprintf("    ESS(sigma_m)               = %.1f / %d\n", ess_sm, n_keep))
  cat(sprintf("    ESS(rho)                   = %.1f / %d\n", ess_rh, n_keep))
  cat(sprintf("    ESS/sec (sigma_m)          = %.2f\n",
              ess_sm / fit$elapsed_sec))
  cat(sprintf("    ESS/sec (rho)              = %.2f\n",
              ess_rh / fit$elapsed_sec))
  if (!is.null(fit$log_c_psi_draws)) {
    cat(sprintf("    final c_psi (proposal scaling) = %.3f\n",
                exp(fit$log_c_psi_draws[fit$burn])))
  }
  invisible(list(acc_post = acc_post, ess_sm = ess_sm, ess_rh = ess_rh))
}

cat("\n=== Diagnostics ===")
s_fixed <- summ(fit_fixed, "fixed")
s_adapt <- summ(fit_adapt, "adapt")

# ACR on grid
eval_acr <- function(fit) {
  k    <- (fit$burn + 1):fit$n_iter
  ep   <- fit$eta_pred_draws[k, , drop = FALSE]
  q025 <- apply(ep, 2, quantile, probs = 0.025)
  q975 <- apply(ep, 2, quantile, probs = 0.975)
  list(
    acr   = mean(eta_grid >= q025 & eta_grid <= q975),
    ail   = mean(q975 - q025),
    rmse  = sqrt(mean((colMeans(ep) - eta_grid)^2))
  )
}
cat("\n=== Grid-eta metrics (eta_true = (1-omega)*z* + omega*f0) ===\n")
mf <- eval_acr(fit_fixed); ma <- eval_acr(fit_adapt)
cat(sprintf("  %-7s  ACR=%.3f  AIL=%.3f  RMSE=%.3f\n",
            "fixed", mf$acr, mf$ail, mf$rmse))
cat(sprintf("  %-7s  ACR=%.3f  AIL=%.3f  RMSE=%.3f\n",
            "adapt", ma$acr, ma$ail, ma$rmse))

saveRDS(list(fixed = fit_fixed, adapt = fit_adapt,
             eta_grid = eta_grid),
        "dev/18_adaptive_psi.rds")
cat("\nSaved dev/18_adaptive_psi.rds\n")
