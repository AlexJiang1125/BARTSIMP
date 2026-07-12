# =============================================================================
# Smoke test for the cluster scripts after porting per-iter grid predictions
# (X_pred / locs_pred / eta_pred_draws / f_pred_draws / z_pred_draws).
#
# Runs ONE small task by sourcing the cluster's bartsimp_pg.R + sim_gen_data.R
# directly (the cluster sim_run_task.R hard-codes /users/zj2451/... output
# paths, so we mirror just the parts that exercise the new code path).
# =============================================================================

cluster_R <- "/Users/alexziyujiang/Documents/GitHub/BARTSIMP_glm/R"
source(file.path(cluster_R, "bartsimp_pg.R"))
source(file.path(cluster_R, "sim_gen_data.R"))

set.seed(20260529)

dat <- sim_gen_dataset(seed = 20260529, surface = "tree", omega = 0.5,
                       n_cl = 150L, grid_dim = 15L)

X_grid    <- as.matrix(dat$grid[, c("x1", "x2")])
locs_grid <- as.matrix(dat$grid[, c("lon", "lat")])

cat(sprintf("Cluster smoke: n_cl=%d grid=%d cells (15x15)\n",
            nrow(dat$X), nrow(dat$grid)))

t0 <- Sys.time()
fit <- run_sampler_bartsimp_pg(
  X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
  f_true = dat$f_train, z_true = dat$z_train,
  X_pred = X_grid, locs_pred = locs_grid,
  n_iter = 800L, burn = 300L, m_trees = 20L,
  mesh_max_edge = c(0.08, 0.25), mesh_cutoff = 0.03,
  bd_mode = "conditional",
  verbose_every = 200L, label = "cluster-smoke"
)
cat(sprintf("Elapsed: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

stopifnot(!is.null(fit$eta_pred_draws),
          !is.null(fit$f_pred_draws),
          !is.null(fit$z_pred_draws),
          all(dim(fit$eta_pred_draws) == c(800L, nrow(X_grid))))

keep <- (fit$burn + 1):fit$n_iter

# Replicate the new sim_run_task.R metric block (CORRECT + BUGGY + DISCRETE)
eta_pred      <- fit$eta_pred_draws[keep, , drop = FALSE]
eta_grid_pt   <- colMeans(eta_pred)
eta_grid_q025 <- apply(eta_pred, 2, quantile, probs = 0.025)
eta_grid_q975 <- apply(eta_pred, 2, quantile, probs = 0.975)

f_grid_pt    <- ensemble_predict(fit$ensemble_final, X_grid)
z_pred       <- fit$z_pred_draws[keep, , drop = FALSE]
z_grid_q025  <- apply(z_pred, 2, quantile, probs = 0.025)
z_grid_q975  <- apply(z_pred, 2, quantile, probs = 0.975)
eta_buggy_pt   <- f_grid_pt + colMeans(z_pred)
eta_buggy_q025 <- f_grid_pt + z_grid_q025
eta_buggy_q975 <- f_grid_pt + z_grid_q975

eta_true <- dat$grid$eta_true
acr      <- mean(eta_true >= eta_grid_q025  & eta_true <= eta_grid_q975)
acr_b    <- mean(eta_true >= eta_buggy_q025 & eta_true <= eta_buggy_q975)

# Discrete: ACR(p) on grid + posterior predictive ACR(k) at training
p_true     <- dat$grid$p_true
p_grid_q025 <- plogis(eta_grid_q025)
p_grid_q975 <- plogis(eta_grid_q975)
acr_p_grid  <- mean(p_true >= p_grid_q025 & p_true <= p_grid_q975)

n_keep <- length(keep)
n_cl <- length(dat$n_i)
p_train_draws <- plogis(fit$f_draws[keep, , drop = FALSE] +
                        fit$z_draws[keep, , drop = FALSE])
n_i_mat <- matrix(dat$n_i, nrow = n_keep, ncol = n_cl, byrow = TRUE)
k_pp <- matrix(rbinom(n_keep * n_cl, as.integer(n_i_mat), p_train_draws),
               nrow = n_keep, ncol = n_cl)
k_pp_q025 <- apply(k_pp, 2, quantile, probs = 0.025)
k_pp_q975 <- apply(k_pp, 2, quantile, probs = 0.975)
acr_k_pp  <- mean(dat$k_i >= k_pp_q025 & dat$k_i <= k_pp_q975)

cat(sprintf("\n  CORRECT  ACR(eta)=%.3f\n", acr))
cat(sprintf("  BUGGY    ACR(eta)=%.3f\n", acr_b))
cat(sprintf("  DISCRETE ACR(p_grid)=%.3f  ACR(k_pp_train)=%.3f\n",
            acr_p_grid, acr_k_pp))

cat("\nSmoke test PASSED.\n")
