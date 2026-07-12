# =============================================================================
# Validation V3 — predictive comparison on held-out test set
# =============================================================================
# Methods compared:
#   M1: BARTSIMP-PG (our R prototype)
#   M2: SPDE-binomial    -- INLA linear-in-covariates + spatial
#   M3: SPDE0-binomial   -- INLA intercept-only + spatial
#   M4: BART-no-spatial  -- PG-augmented BART with no spatial random field
#
# Metrics on held-out test clusters:
#   - eta RMSE (linear predictor scale)
#   - p RMSE (probability scale)
#   - Brier score on success proportions k_test / n_test
#   - log-loss (binomial)
#
# SMOKE_TEST -> tiny dataset + few iterations to verify the pipelines run.
# =============================================================================

source("dev/bartsimp_pg.R")

SMOKE_TEST <- TRUE   # flip to FALSE for production

if (SMOKE_TEST) {
  cat("*** SMOKE TEST MODE ***\n")
  N_TRAIN <- 80
  N_TEST  <- 20
  N_ITER  <- 150
  BURN    <- 50
  M_TREES <- 5
  MESH_EDGE <- c(0.15, 0.40)
  METHODS <- c("BARTSIMP_PG", "SPDE_binom", "SPDE0_binom", "BART_only")
} else {
  N_TRAIN <- 200
  N_TEST  <- 50
  N_ITER  <- 3000
  BURN    <- 1000
  M_TREES <- 20
  MESH_EDGE <- c(0.08, 0.25)
  METHODS <- c("BARTSIMP_PG", "SPDE_binom", "SPDE0_binom", "BART_only")
}

cat(sprintf("Config: n_train=%d, n_test=%d, n_iter=%d, m_trees=%d\n",
            N_TRAIN, N_TEST, N_ITER, M_TREES))
cat("Methods:", paste(METHODS, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Data: balanced scenario, simulate train + test together so test eta_true known
# -----------------------------------------------------------------------------
set.seed(2026)
dat <- gen_data(seed = 2026, n_cl = N_TRAIN + N_TEST,
                covariate_signal = TRUE, spatial_signal = TRUE)
train_idx <- seq_len(N_TRAIN)
test_idx  <- (N_TRAIN + 1):(N_TRAIN + N_TEST)

dat_train <- list(X     = dat$X[train_idx, , drop = FALSE],
                  locs  = dat$locs[train_idx, , drop = FALSE],
                  n_i   = dat$n_i[train_idx],
                  k_i   = dat$k_i[train_idx],
                  f_true = dat$f_true[train_idx],
                  z_true = dat$z_true[train_idx])
dat_test  <- list(X        = dat$X[test_idx, , drop = FALSE],
                  locs     = dat$locs[test_idx, , drop = FALSE],
                  n_i      = dat$n_i[test_idx],
                  k_i      = dat$k_i[test_idx],
                  f_true   = dat$f_true[test_idx],
                  z_true   = dat$z_true[test_idx],
                  eta_true = dat$eta_true[test_idx],
                  p_true   = dat$p_true[test_idx])

# -----------------------------------------------------------------------------
# M1: BARTSIMP-PG
# -----------------------------------------------------------------------------
fit_bartsimp_pg <- function(dat_train, dat_test) {
  fit <- run_sampler_bartsimp_pg(
    X = dat_train$X, locs = dat_train$locs,
    n_i = dat_train$n_i, k_i = dat_train$k_i,
    n_iter = N_ITER, burn = BURN, m_trees = M_TREES,
    mesh_max_edge = MESH_EDGE, verbose_every = 0,
    label = "M1-BARTSIMP-PG"
  )
  # Predict at test locations: f from final ensemble, z from posterior-mean u
  A_test <- as(INLA::inla.spde.make.A(mesh = fit$mesh, loc = dat_test$locs),
               "CsparseMatrix")
  keep <- (fit$burn + 1):fit$n_iter
  u_mean <- colMeans(fit$u_draws[keep, , drop = FALSE])
  z_test <- as.numeric(A_test %*% u_mean)
  f_test <- ensemble_predict(fit$ensemble_final, dat_test$X)
  list(eta_test = f_test + z_test, elapsed_sec = fit$elapsed_sec,
       f_test = f_test, z_test = z_test)
}

# -----------------------------------------------------------------------------
# M2/M3 helper: INLA binomial fit with optional covariates
# -----------------------------------------------------------------------------
fit_inla_binom <- function(dat_train, dat_test, include_covariates = TRUE) {
  loc_all <- rbind(dat_train$locs, dat_test$locs)
  mesh <- INLA::inla.mesh.2d(loc = loc_all, max.edge = MESH_EDGE, cutoff = 0.03)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.1, 0.5),
                                          prior.sigma = c(1.0, 0.5))
  A_train <- INLA::inla.spde.make.A(mesh = mesh, loc = dat_train$locs)
  A_test  <- INLA::inla.spde.make.A(mesh = mesh, loc = dat_test$locs)

  if (include_covariates) {
    eff_train <- list(i = seq_len(spde$n.spde),
                      data.frame(intercept = rep(1, nrow(dat_train$X)),
                                 x1 = dat_train$X[, "x1"],
                                 x2 = dat_train$X[, "x2"]))
    eff_test <- list(i = seq_len(spde$n.spde),
                     data.frame(intercept = rep(1, nrow(dat_test$X)),
                                x1 = dat_test$X[, "x1"],
                                x2 = dat_test$X[, "x2"]))
    formula <- y ~ -1 + intercept + x1 + x2 + f(i, model = spde)
  } else {
    eff_train <- list(i = seq_len(spde$n.spde),
                      data.frame(intercept = rep(1, nrow(dat_train$X))))
    eff_test <- list(i = seq_len(spde$n.spde),
                     data.frame(intercept = rep(1, nrow(dat_test$X))))
    formula <- y ~ -1 + intercept + f(i, model = spde)
  }

  stk_train <- INLA::inla.stack(
    data    = list(y = dat_train$k_i, n_trials = dat_train$n_i),
    A       = list(A_train, 1),
    effects = eff_train,
    tag     = "train"
  )
  stk_test <- INLA::inla.stack(
    data    = list(y = rep(NA, nrow(dat_test$X)),
                   n_trials = rep(NA, nrow(dat_test$X))),
    A       = list(A_test, 1),
    effects = eff_test,
    tag     = "test"
  )
  stk <- INLA::inla.stack(stk_train, stk_test)

  t0 <- Sys.time()
  res <- INLA::inla(
    formula,
    family = "binomial",
    Ntrials = INLA::inla.stack.data(stk)$n_trials,
    data    = INLA::inla.stack.data(stk),
    control.predictor = list(A = INLA::inla.stack.A(stk),
                              compute = TRUE, link = 1),
    control.compute   = list(dic = FALSE, waic = FALSE)
  )
  t1 <- Sys.time()
  idx_test <- INLA::inla.stack.index(stk, "test")$data
  eta_test <- res$summary.linear.predictor[idx_test, "mean"]
  list(eta_test = eta_test,
       elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")))
}

# -----------------------------------------------------------------------------
# M4: BART with no spatial component (PG-augmented logistic BART)
# Simply run our sampler but force sigma_m near zero (sampler will still
# update; an alternative is to clamp z=0 directly).
# We do a minimal pure-BART-with-PG: no spatial draws, fixed z=0.
# -----------------------------------------------------------------------------
fit_bart_only <- function(dat_train, dat_test) {
  n_cl <- nrow(dat_train$X)
  tau2 <- (3 / (2 * sqrt(M_TREES)))^2
  ensemble  <- ensemble_init(M_TREES, init_val = 0)
  f_current <- rep(0, n_cl)

  t0 <- Sys.time()
  for (it in seq_len(N_ITER)) {
    eta   <- f_current             # z fixed at 0
    Omega <- rpg(num = n_cl, h = dat_train$n_i, z = eta)
    ybar  <- (dat_train$k_i - dat_train$n_i / 2) / Omega
    sw <- ensemble_sweep(ensemble, dat_train$X, ybar, Omega, tau2)
    ensemble  <- sw$ensemble
    f_current <- sw$fit
  }
  t1 <- Sys.time()
  f_test <- ensemble_predict(ensemble, dat_test$X)
  list(eta_test = f_test,
       elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")))
}

# -----------------------------------------------------------------------------
# Run all selected methods + collect metrics
# -----------------------------------------------------------------------------
metrics <- function(eta_pred, eta_true, p_true, k, n) {
  p_pred <- plogis(eta_pred)
  list(
    eta_rmse  = sqrt(mean((eta_pred - eta_true)^2)),
    p_rmse    = sqrt(mean((p_pred  - p_true)^2)),
    brier     = mean((p_pred - (k / n))^2),
    log_loss  = -mean(k * log(pmax(p_pred, 1e-10)) +
                      (n - k) * log(pmax(1 - p_pred, 1e-10)))
  )
}

results <- list()
for (m in METHODS) {
  cat(sprintf("\n--- Fitting %s ---\n", m))
  out <- tryCatch({
    switch(m,
      BARTSIMP_PG = fit_bartsimp_pg(dat_train, dat_test),
      SPDE_binom  = fit_inla_binom(dat_train, dat_test, include_covariates = TRUE),
      SPDE0_binom = fit_inla_binom(dat_train, dat_test, include_covariates = FALSE),
      BART_only   = fit_bart_only(dat_train, dat_test)
    )
  }, error = function(e) {
    cat(sprintf("  [%s] FAILED: %s\n", m, conditionMessage(e)))
    NULL
  })

  if (!is.null(out)) {
    mt <- metrics(out$eta_test, dat_test$eta_true, dat_test$p_true,
                  dat_test$k_i, dat_test$n_i)
    cat(sprintf("  eta_rmse=%.3f  p_rmse=%.3f  Brier=%.4f  log_loss=%.3f  time=%.2fs\n",
                mt$eta_rmse, mt$p_rmse, mt$brier, mt$log_loss, out$elapsed_sec))
    results[[m]] <- c(method = m, mt, elapsed_sec = out$elapsed_sec)
  }
}

# Pretty table
tbl <- do.call(rbind, lapply(results, function(r) {
  data.frame(method      = r$method,
             eta_rmse    = r$eta_rmse,
             p_rmse      = r$p_rmse,
             brier       = r$brier,
             log_loss    = r$log_loss,
             elapsed_sec = r$elapsed_sec)
}))
cat("\n=== V3 comparison summary ===\n")
print(tbl, row.names = FALSE)

saveRDS(list(results = results, tbl = tbl,
             config = list(N_TRAIN = N_TRAIN, N_TEST = N_TEST,
                           N_ITER = N_ITER, METHODS = METHODS,
                           SMOKE_TEST = SMOKE_TEST)),
        sprintf("dev/07_validation_v3_%s.rds",
                if (SMOKE_TEST) "smoke" else "full"))
cat(sprintf("\nSaved to dev/07_validation_v3_%s.rds\n",
            if (SMOKE_TEST) "smoke" else "full"))
