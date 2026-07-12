# =============================================================================
# RSR vs no-RSR head-to-head, paper-aligned generator + adaptive MH.
#
# Same dataset, same seed, both samplers run with adaptive psi MH on.
# The only difference: rsr_basis = "linear_interaction" vs NULL.
#
# Question: does RSR close the sigma_m underestimation gap (target 0.5,
# currently ~0.3 with adaptive MH alone)?
#
# Pass criteria for RSR to be worth deploying:
#   (a) sigma_m mean closer to truth (>= 0.4 ideally)
#   (b) ACR stays in [0.93, 0.99]
#   (c) RMSE(eta) doesn't regress noticeably
# =============================================================================

source("dev/bartsimp_pg.R")
source("dev/rbart.R")

set.seed(2026)

# ---- toy config -----------------------------------------------------------
N_CL     <- 150L
N_ITER   <- 1500L
BURN     <-  500L
M_TREES  <-   20L
GRID_DIM <-   20L

# Paper-aligned generator (same as dev/13_acr_diagnosis.R and dev/18_adaptive_psi.R)
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

# Helper to run + summarize
run_one <- function(rsr_basis, label) {
  cat(sprintf("\n--- %s ---\n", label))
  set.seed(2026)   # same RNG path so chains are directly comparable
  t0 <- Sys.time()
  fit <- run_sampler_bartsimp_pg(
    X = X, locs = locs, n_i = n_i, k_i = k_i,
    f_true = f_train, z_true = z_train,
    X_pred = X_grid, locs_pred = locs_grid,
    n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
    mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
    bd_mode = "conditional",
    adapt_psi = TRUE, psi_target_acc = 0.234,
    rsr_basis = rsr_basis,
    verbose_every = 250, label = label
  )
  cat(sprintf("Elapsed: %.1fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  fit
}

fit_none <- run_one(NULL,                  "no-RSR")
fit_rsr  <- run_one("linear_interaction",  "RSR-linint")

# ---- diagnostics ----------------------------------------------------------
summ <- function(fit, lbl) {
  keep <- (fit$burn + 1):fit$n_iter
  ep   <- fit$eta_pred_draws[keep, , drop = FALSE]
  q025 <- apply(ep, 2, quantile, probs = 0.025)
  q975 <- apply(ep, 2, quantile, probs = 0.975)
  acr  <- mean(eta_grid >= q025 & eta_grid <= q975)
  ail  <- mean(q975 - q025)
  rmse <- sqrt(mean((colMeans(ep) - eta_grid)^2))
  cat(sprintf("\n  [%s]\n", lbl))
  cat(sprintf("    sigma_m mean (truth 0.500) = %.3f\n",
              mean(fit$sigma_m_draws[keep])))
  cat(sprintf("    sigma_m 95%% CI            = [%.3f, %.3f]\n",
              quantile(fit$sigma_m_draws[keep], 0.025),
              quantile(fit$sigma_m_draws[keep], 0.975)))
  cat(sprintf("    rho     mean (truth 1.131) = %.3f\n",
              mean(fit$rho_draws[keep])))
  cat(sprintf("    psi MH acc (post-burn)     = %.3f\n",
              mean(fit$accept_psi[keep])))
  cat(sprintf("    var(f) mean                = %.3f\n",
              mean(fit$varpart_draws[keep, "var_f"])))
  cat(sprintf("    var(z) mean                = %.3f\n",
              mean(fit$varpart_draws[keep, "var_z"])))
  cat(sprintf("    cor(f, z) post mean        = %.3f\n",
              mean(fit$fz_cor_draws[keep])))
  cat(sprintf("    ACR(eta) grid              = %.3f\n", acr))
  cat(sprintf("    AIL(eta) grid              = %.3f\n", ail))
  cat(sprintf("    RMSE(eta) grid             = %.3f\n", rmse))
  invisible(list(acr = acr, ail = ail, rmse = rmse,
                 sigma_m_mean = mean(fit$sigma_m_draws[keep])))
}

cat("\n=================== RSR vs no-RSR ===================")
s_none <- summ(fit_none, "no-RSR")
s_rsr  <- summ(fit_rsr,  "RSR-linint")

cat("\n=== Delta (RSR - no-RSR) ===\n")
cat(sprintf("  sigma_m_mean: %+.3f   (positive = RSR closer to truth 0.5)\n",
            s_rsr$sigma_m_mean - s_none$sigma_m_mean))
cat(sprintf("  ACR         : %+.3f\n", s_rsr$acr - s_none$acr))
cat(sprintf("  AIL         : %+.3f\n", s_rsr$ail - s_none$ail))
cat(sprintf("  RMSE(eta)   : %+.3f\n", s_rsr$rmse - s_none$rmse))

saveRDS(list(fit_none = fit_none, fit_rsr = fit_rsr,
             eta_grid = eta_grid,
             X = X, locs = locs, n_i = n_i, k_i = k_i,
             X_grid = X_grid, locs_grid = locs_grid),
        "dev/19_rsr_test.rds")
cat("\nSaved dev/19_rsr_test.rds\n")
