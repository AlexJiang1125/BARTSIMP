# =============================================================================
# ACR diagnosis: show the systematic undercoverage from the cluster sim study
# is an interval-construction artifact, not a sampler problem.
#
# The cluster's sim_run_task.R builds grid intervals as
#     eta_q025 = f_grid_point + z_grid_q025
#     eta_q975 = f_grid_point + z_grid_q975
# i.e. it treats BART as known and only propagates spatial uncertainty.
#
# The correct construction uses per-iteration eta = f_iter + z_iter on the grid
# and quantiles those directly. This script runs ONE toy dataset and reports
# ACR both ways so the gap can be measured.
# =============================================================================

source("dev/bartsimp_pg.R")
source("dev/rbart.R")

set.seed(2026)

# ---- toy config (fast, but representative) ---------------------------------
N_CL     <- 150L
N_ITER   <- 1500L
BURN     <-  500L
M_TREES  <-   20L
GRID_DIM <-   20L

# Generate one balanced-scenario dataset PLUS a held-out grid using the cluster
# sim_gen_data convention (jointly drawn Matern GRF on (locs, grid), then split).
# This mirrors what sim_run_task.R does on the cluster -- we just inline it
# here so the script is self-contained.

sigma_m_baseline <- 0.5
rho_baseline     <- sqrt(8 * 1) / 2.5    # nu=1, range as in the paper

# Surface: "tree" step function (matches surface_f0 in sim_gen_data.R)
surface_f0 <- function(x1, x2)
  ifelse(x1 < 0.5, 3, ifelse(x2 < 0.5, -2, 0))

# Use a balanced omega so both BART and the GRF contribute meaningfully.
OMEGA <- 0.5

# Grid in [0,1]^2
G    <- GRID_DIM^2
step <- 1 / GRID_DIM
gx   <- (seq_len(GRID_DIM) - 0.5) * step
cell_centers <- expand.grid(lon = gx, lat = gx)

# Per-cell covariates
x1_grid <- runif(G)
x2_grid <- runif(G)

# Joint Matern GRF on grid cells
Dg    <- as.matrix(dist(cell_centers))
kappa <- sqrt(8 * 1) / rho_baseline
Sg    <- sigma_m_baseline^2 * matern_corr(Dg, kappa, 1) + 1e-6 * diag(G)
z_grid <- as.vector(t(chol(Sg)) %*% rnorm(G))

# Truth on grid
f0_grid  <- surface_f0(x1_grid, x2_grid)
eta_grid <- (1 - OMEGA) * z_grid + OMEGA * f0_grid

# Sample N_CL clusters (one per random cell) with a small jitter inside the cell
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

# Training truths needed by run_sampler_bartsimp_pg's RMSE printout
f_train <- OMEGA * surface_f0(X[, 1], X[, 2])
z_train <- (1 - OMEGA) * z_grid[cluster_cells]

X_grid    <- cbind(x1 = x1_grid, x2 = x2_grid)
locs_grid <- as.matrix(cell_centers[, c("lon", "lat")])

cat(sprintf("Dataset: n_cl=%d, grid=%dx%d=%d cells, omega=%.2f\n",
            N_CL, GRID_DIM, GRID_DIM, G, OMEGA))

# ---- run sampler with per-iter grid prediction enabled ---------------------
cat("\n--- Running sampler (Path B, conditional BD) ---\n")
t0 <- Sys.time()
fit <- run_sampler_bartsimp_pg(
  X = X, locs = locs, n_i = n_i, k_i = k_i,
  f_true = f_train, z_true = z_train,
  X_pred = X_grid, locs_pred = locs_grid,
  n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
  mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
  bd_mode = "conditional",
  verbose_every = 250, label = "acr-toy"
)
cat(sprintf("Elapsed: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- ACR / AIL the two ways ------------------------------------------------
keep <- (fit$burn + 1):fit$n_iter
alpha <- 0.05

eta_true <- eta_grid    # truth on the grid

# (1) BUGGY: cluster sim_run_task.R style
#     f_grid_pt = point estimate from the final BART ensemble
#     z_grid_qXXX = quantiles of z_pred_draws  (proper spatial uncertainty)
#     eta_qXXX = f_grid_pt + z_grid_qXXX
f_grid_pt <- ensemble_predict(fit$ensemble_final, X_grid)
z_pred    <- fit$z_pred_draws[keep, , drop = FALSE]
z_q025    <- apply(z_pred, 2, quantile, probs = 0.025)
z_q975    <- apply(z_pred, 2, quantile, probs = 0.975)
eta_buggy_pt   <- f_grid_pt + colMeans(z_pred)
eta_buggy_q025 <- f_grid_pt + z_q025
eta_buggy_q975 <- f_grid_pt + z_q975

acr_buggy <- mean(eta_true >= eta_buggy_q025 & eta_true <= eta_buggy_q975)
ail_buggy <- mean(eta_buggy_q975 - eta_buggy_q025)
rmse_buggy <- sqrt(mean((eta_buggy_pt - eta_true)^2))

# (2) CORRECT: per-iter eta = f_iter + z_iter on the grid, then quantile
eta_pred <- fit$eta_pred_draws[keep, , drop = FALSE]
eta_correct_pt   <- colMeans(eta_pred)
eta_correct_q025 <- apply(eta_pred, 2, quantile, probs = 0.025)
eta_correct_q975 <- apply(eta_pred, 2, quantile, probs = 0.975)

acr_correct <- mean(eta_true >= eta_correct_q025 &
                    eta_true <= eta_correct_q975)
ail_correct <- mean(eta_correct_q975 - eta_correct_q025)
rmse_correct <- sqrt(mean((eta_correct_pt - eta_true)^2))

# AIS (Gneiting-Raftery interval score) both ways
ais_fun <- function(true, q025, q975, a = alpha) {
  width <- q975 - q025
  mean(width
       + (2/a) * (q025 - true) * (true < q025)
       + (2/a) * (true - q975) * (true > q975))
}
ais_buggy   <- ais_fun(eta_true, eta_buggy_q025,   eta_buggy_q975)
ais_correct <- ais_fun(eta_true, eta_correct_q025, eta_correct_q975)

# Also: BART-only interval (how much f-uncertainty is there?)
f_pred <- fit$f_pred_draws[keep, , drop = FALSE]
f_q025 <- apply(f_pred, 2, quantile, probs = 0.025)
f_q975 <- apply(f_pred, 2, quantile, probs = 0.975)
ail_f  <- mean(f_q975 - f_q025)
ail_z  <- mean(z_q975 - z_q025)

# Coverage decomposition
f_only_acr <- mean(OMEGA * f0_grid >= f_q025 &
                   OMEGA * f0_grid <= f_q975)
z_only_acr <- mean((1 - OMEGA) * z_grid >= z_q025 &
                   (1 - OMEGA) * z_grid <= z_q975)

# ---- discrete-outcome coverage (independent double-check) ------------------
# (a) ACR(p) on grid: monotonic plogis transform of eta; should equal ACR(eta).
p_pred       <- plogis(eta_pred)               # n_keep x G
p_correct_q025 <- apply(p_pred, 2, quantile, probs = 0.025)
p_correct_q975 <- apply(p_pred, 2, quantile, probs = 0.975)
p_true       <- plogis(eta_true)
acr_p_correct <- mean(p_true >= p_correct_q025 & p_true <= p_correct_q975)
ail_p_correct <- mean(p_correct_q975 - p_correct_q025)

# Same on the buggy reconstruction (plogis is monotonic so equal to eta ACR)
p_buggy_q025 <- plogis(eta_buggy_q025)
p_buggy_q975 <- plogis(eta_buggy_q975)
acr_p_buggy  <- mean(p_true >= p_buggy_q025 & p_true <= p_buggy_q975)
ail_p_buggy  <- mean(p_buggy_q975 - p_buggy_q025)

# (b) Posterior predictive ACR of observed k_i at training clusters.
# Per iter, draw k_iter ~ Binomial(n_i, plogis(f_iter + z_iter)) for each
# cluster, then take 2.5/97.5 quantiles across iterations. Coverage of the
# observed k_i tests that the model is calibrated on the data we actually saw,
# INCLUDING binomial sampling variability (independent of the grid checks).
f_train_draws <- fit$f_draws[keep, , drop = FALSE]
z_train_draws <- fit$z_draws[keep, , drop = FALSE]
eta_train_draws <- f_train_draws + z_train_draws
p_train_draws   <- plogis(eta_train_draws)     # n_keep x n_cl

n_keep <- length(keep)
n_i_mat <- matrix(n_i, nrow = n_keep, ncol = N_CL, byrow = TRUE)
k_pp_draws <- matrix(rbinom(n_keep * N_CL,
                            size = as.integer(n_i_mat),
                            prob = p_train_draws),
                     nrow = n_keep, ncol = N_CL)
k_pp_q025 <- apply(k_pp_draws, 2, quantile, probs = 0.025)
k_pp_q975 <- apply(k_pp_draws, 2, quantile, probs = 0.975)
acr_k_pp  <- mean(k_i >= k_pp_q025 & k_i <= k_pp_q975)
ail_k_pp  <- mean(k_pp_q975 - k_pp_q025)

# Latent-p ACR at training clusters (no binomial noise; truth = p_train).
p_train_q025 <- apply(p_train_draws, 2, quantile, probs = 0.025)
p_train_q975 <- apply(p_train_draws, 2, quantile, probs = 0.975)
p_train_true <- plogis(eta_train)
acr_p_train  <- mean(p_train_true >= p_train_q025 & p_train_true <= p_train_q975)
ail_p_train  <- mean(p_train_q975 - p_train_q025)

# ---- report ----------------------------------------------------------------
cat("\n========================================================\n")
cat("ACR diagnosis: omega = ", OMEGA,
    ", surface = tree, n_cl = ", N_CL, "\n", sep = "")
cat("========================================================\n")
cat(sprintf("\n%-25s   %-12s %-12s\n",
            "metric", "buggy(cluster)", "correct(per-iter)"))
cat(sprintf("%-25s   %-12.3f %-12.3f\n",
            "RMSE(eta)",  rmse_buggy, rmse_correct))
cat(sprintf("%-25s   %-12.3f %-12.3f\n",
            "AIL (eta)",  ail_buggy, ail_correct))
cat(sprintf("%-25s   %-12.3f %-12.3f\n",
            "ACR (eta)",  acr_buggy, acr_correct))
cat(sprintf("%-25s   %-12.3f %-12.3f\n",
            "AIS (eta)",  ais_buggy, ais_correct))

cat("\nComponent diagnostics\n")
cat(sprintf("  AIL(f) at grid        = %.3f\n", ail_f))
cat(sprintf("  AIL(z) at grid        = %.3f\n", ail_z))
cat(sprintf("  ACR(f) vs omega*f0    = %.3f\n", f_only_acr))
cat(sprintf("  ACR(z) vs (1-w)*z*    = %.3f\n", z_only_acr))

cat("\nDiscrete-outcome coverage (independent double-check)\n")
cat(sprintf("%-35s   %-12s %-12s\n",
            "metric", "buggy", "correct"))
cat(sprintf("%-35s   %-12.3f %-12.3f\n",
            "ACR(p) on grid",  acr_p_buggy, acr_p_correct))
cat(sprintf("%-35s   %-12.3f %-12.3f\n",
            "AIL(p) on grid",  ail_p_buggy, ail_p_correct))
cat(sprintf("  ACR(p) at training (latent)         = %.3f\n",
            acr_p_train))
cat(sprintf("  AIL(p) at training (latent)         = %.3f\n",
            ail_p_train))
cat(sprintf("  ACR(k) posterior predictive @ train = %.3f   (cov of observed k_i)\n",
            acr_k_pp))
cat(sprintf("  AIL(k) posterior predictive @ train = %.3f\n",
            ail_k_pp))

cat("\nSampler diagnostics (post-burn)\n")
cat(sprintf("  sigma_m posterior mean = %.3f   (truth = %.3f)\n",
            mean(fit$sigma_m_draws[keep]), sigma_m_baseline))
cat(sprintf("  rho     posterior mean = %.3f   (truth = %.3f)\n",
            mean(fit$rho_draws[keep]), rho_baseline))
cat(sprintf("  psi MH acceptance      = %.3f   (target ~ 0.234)\n",
            mean(fit$accept_psi)))
cat(sprintf("  cor(f, z) posterior    = %.3f\n",
            mean(fit$fz_cor_draws[keep])))

saveRDS(
  list(
    config = list(N_CL = N_CL, N_ITER = N_ITER, BURN = BURN,
                  M_TREES = M_TREES, GRID_DIM = GRID_DIM, OMEGA = OMEGA),
    truth = list(eta_grid = eta_grid, f0_grid = f0_grid, z_grid = z_grid,
                 sigma_m = sigma_m_baseline, rho = rho_baseline),
    buggy = list(pt = eta_buggy_pt, q025 = eta_buggy_q025, q975 = eta_buggy_q975,
                 acr = acr_buggy, ail = ail_buggy, ais = ais_buggy,
                 rmse = rmse_buggy),
    correct = list(pt = eta_correct_pt, q025 = eta_correct_q025,
                   q975 = eta_correct_q975,
                   acr = acr_correct, ail = ail_correct, ais = ais_correct,
                   rmse = rmse_correct),
    f_grid_pt = f_grid_pt,
    components = list(ail_f = ail_f, ail_z = ail_z,
                      acr_f = f_only_acr, acr_z = z_only_acr),
    discrete = list(
      acr_p_grid_buggy   = acr_p_buggy,
      acr_p_grid_correct = acr_p_correct,
      ail_p_grid_buggy   = ail_p_buggy,
      ail_p_grid_correct = ail_p_correct,
      acr_p_train        = acr_p_train,
      ail_p_train        = ail_p_train,
      acr_k_pp_train     = acr_k_pp,
      ail_k_pp_train     = ail_k_pp
    ),
    fit_summary = list(
      sigma_m_mean = mean(fit$sigma_m_draws[keep]),
      rho_mean     = mean(fit$rho_draws[keep]),
      psi_acc      = mean(fit$accept_psi),
      fz_cor       = mean(fit$fz_cor_draws[keep]),
      elapsed_sec  = fit$elapsed_sec
    )
  ),
  "dev/13_acr_diagnosis.rds"
)
cat("\nSaved dev/13_acr_diagnosis.rds\n")
