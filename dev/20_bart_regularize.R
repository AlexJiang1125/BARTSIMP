# =============================================================================
# BART regularization head-to-head: does (m_trees=100, k_bart=3) close the
# sigma_m gap that (m_trees=20, k_bart=2) leaves open?
#
# Same dataset/seed as dev/18 and dev/19. Adaptive MH on in both. Only the
# two BART regularization knobs change.
#
# Pass criteria:
#   (a) sigma_m mean closer to truth (>= 0.4)
#   (b) var(f) closer to truth (~ 0.6 in this balanced scenario)
#   (c) ACR stays in [0.93, 0.99]
#   (d) RMSE(eta) doesn't regress
#
# Cost: m=100 trees is ~5x the work per sweep. Expect 700s vs 150s wall.
# =============================================================================

source("dev/bartsimp_pg.R")
source("dev/rbart.R")

set.seed(2026)

# ---- toy config (same as dev/18, dev/19) ----------------------------------
N_CL     <- 150L
N_ITER   <- 1500L
BURN     <-  500L
GRID_DIM <-   20L

# Paper-aligned generator
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

# truth: var(omega*f0(X_train)) for the n_cl training X's
var_f_truth <- var(OMEGA * surface_f0(X[, 1], X[, 2]))
var_z_truth <- var((1 - OMEGA) * z_grid[cluster_cells])

cat(sprintf("Dataset: n_cl=%d, grid=%d cells, omega=%.2f\n",
            N_CL, G, OMEGA))
cat(sprintf("Truth: var(f) at train = %.3f   var(z) at train = %.3f\n",
            var_f_truth, var_z_truth))

# Helper: same harness, swap m_trees + k_bart
run_one <- function(m_trees, k_bart, label) {
  cat(sprintf("\n--- %s (m=%d, k=%g) ---\n", label, m_trees, k_bart))
  set.seed(2026)
  t0 <- Sys.time()
  fit <- run_sampler_bartsimp_pg(
    X = X, locs = locs, n_i = n_i, k_i = k_i,
    f_true = f_train, z_true = z_train,
    X_pred = X_grid, locs_pred = locs_grid,
    n_iter = N_ITER, burn = BURN,
    m_trees = m_trees, k_bart = k_bart,
    mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
    bd_mode = "conditional",
    adapt_psi = TRUE, psi_target_acc = 0.234,
    verbose_every = 250, label = label
  )
  cat(sprintf("Elapsed: %.1fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  fit
}

fit_default <- run_one(20L, 2,  "m20-k2")
fit_reg     <- run_one(100L, 3, "m100-k3")

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
  cat(sprintf("    var(f) mean (truth %.3f)   = %.3f\n",
              var_f_truth, mean(fit$varpart_draws[keep, "var_f"])))
  cat(sprintf("    var(z) mean (truth %.3f)   = %.3f\n",
              var_z_truth, mean(fit$varpart_draws[keep, "var_z"])))
  cat(sprintf("    cor(f, z) post mean        = %.3f\n",
              mean(fit$fz_cor_draws[keep])))
  cat(sprintf("    ACR(eta) grid              = %.3f\n", acr))
  cat(sprintf("    AIL(eta) grid              = %.3f\n", ail))
  cat(sprintf("    RMSE(eta) grid             = %.3f\n", rmse))
  cat(sprintf("    elapsed                    = %.1fs\n", fit$elapsed_sec))
  invisible(list(
    acr = acr, ail = ail, rmse = rmse,
    sigma_m_mean = mean(fit$sigma_m_draws[keep]),
    var_f_mean   = mean(fit$varpart_draws[keep, "var_f"]),
    var_z_mean   = mean(fit$varpart_draws[keep, "var_z"]),
    elapsed = fit$elapsed_sec
  ))
}

cat("\n=================== BART regularization ===================")
s_def <- summ(fit_default, "m20-k2  (current)")
s_reg <- summ(fit_reg,     "m100-k3 (proposed)")

cat("\n=== Delta (m100-k3 minus m20-k2) ===\n")
cat(sprintf("  sigma_m_mean: %+.3f   (positive = closer to truth 0.5)\n",
            s_reg$sigma_m_mean - s_def$sigma_m_mean))
cat(sprintf("  var(f)      : %+.3f   (negative = closer to truth %.3f)\n",
            s_reg$var_f_mean - s_def$var_f_mean, var_f_truth))
cat(sprintf("  var(z)      : %+.3f\n",
            s_reg$var_z_mean - s_def$var_z_mean))
cat(sprintf("  ACR         : %+.3f\n", s_reg$acr - s_def$acr))
cat(sprintf("  AIL         : %+.3f\n", s_reg$ail - s_def$ail))
cat(sprintf("  RMSE(eta)   : %+.3f\n", s_reg$rmse - s_def$rmse))
cat(sprintf("  wall-clock  : %.1fx slower\n",
            s_reg$elapsed / s_def$elapsed))

saveRDS(list(fit_default = fit_default, fit_reg = fit_reg,
             eta_grid = eta_grid,
             var_f_truth = var_f_truth, var_z_truth = var_z_truth),
        "dev/20_bart_regularize.rds")
cat("\nSaved dev/20_bart_regularize.rds\n")
