# =============================================================================
# SPDE-binomial: linear-in-covariates binomial spatial model via INLA
#
# eta_i = beta0 + beta1 * x1_i + beta2 * x2_i + z_i,
# z ~ Matern via SPDE,
# k_i ~ Binomial(n_i, plogis(eta_i)).
#
# Posterior summaries on the supplied grid (X_pred, locs_pred) are obtained
# by stacking grid rows as held-out predictions in the INLA call.
# =============================================================================

# INLA is loaded lazily at call time so sourcing this file is harmless when
# INLA isn't installed (other competing methods can still run).

fit_spde_binomial <- function(X_train, locs_train, n_i, k_i,
                              X_pred,  locs_pred,
                              n_iter = NULL, burn = NULL,
                              mesh_max_edge = c(0.08, 0.25),
                              mesh_cutoff   = 0.03,
                              prior_range_l = 0.10,
                              prior_range_p = 0.5,
                              prior_sigma_u = 1.0,
                              prior_sigma_p = 0.5,
                              num_threads   = "1:1",
                              verbose = TRUE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("Package 'INLA' is required for fit_spde_binomial() but not installed.")
  stopifnot(length(n_i) == nrow(X_train),
            length(k_i) == nrow(X_train),
            ncol(X_train) == 2L, ncol(X_pred) == 2L)

  n_train <- nrow(X_train); n_pred <- nrow(X_pred)
  if (verbose)
    cat(sprintf("[spde_binomial] n_train=%d, n_pred=%d\n", n_train, n_pred))

  # Mesh built on training locations only (consistent with how BARTSIMP-PG
  # builds its mesh from training locs).
  mesh <- INLA::inla.mesh.2d(loc = locs_train,
                             max.edge = mesh_max_edge,
                             cutoff   = mesh_cutoff)
  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(prior_range_l, prior_range_p),
    prior.sigma = c(prior_sigma_u, prior_sigma_p)
  )

  A_train <- INLA::inla.spde.make.A(mesh = mesh, loc = locs_train)
  A_pred  <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(locs_pred))

  stk_train <- INLA::inla.stack(
    tag = "train",
    data = list(y = k_i, Ntrials = n_i),
    A    = list(A_train, 1),
    effects = list(
      list(spatial = seq_len(mesh$n)),
      list(intercept = rep(1, n_train),
           x1 = X_train[, 1],
           x2 = X_train[, 2])
    )
  )
  stk_pred <- INLA::inla.stack(
    tag = "pred",
    data = list(y = rep(NA_real_, n_pred), Ntrials = rep(1, n_pred)),
    A    = list(A_pred, 1),
    effects = list(
      list(spatial = seq_len(mesh$n)),
      list(intercept = rep(1, n_pred),
           x1 = X_pred[, 1],
           x2 = X_pred[, 2])
    )
  )
  stk <- INLA::inla.stack(stk_train, stk_pred)

  formula <- y ~ -1 + intercept + x1 + x2 + f(spatial, model = spde)

  t0 <- Sys.time()
  fit <- INLA::inla(
    formula,
    family       = "binomial",
    data         = INLA::inla.stack.data(stk),
    Ntrials      = INLA::inla.stack.data(stk)$Ntrials,
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE,
                             link = 1),
    control.compute   = list(return.marginals.predictor = TRUE),
    num.threads  = num_threads,
    verbose      = verbose
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Extract predictions on the grid (pred rows of the stack)
  idx_pred <- INLA::inla.stack.index(stk, tag = "pred")$data
  eta_marg <- fit$marginals.linear.predictor[idx_pred]
  eta_mean <- vapply(eta_marg, function(m)
    INLA::inla.emarginal(function(x) x, m), numeric(1))
  eta_q025 <- vapply(eta_marg, function(m)
    INLA::inla.qmarginal(0.025, m), numeric(1))
  eta_q975 <- vapply(eta_marg, function(m)
    INLA::inla.qmarginal(0.975, m), numeric(1))

  # Probability scale: transform marginals through plogis
  p_mean <- vapply(eta_marg, function(m) {
    mt <- INLA::inla.tmarginal(plogis, m)
    INLA::inla.emarginal(function(x) x, mt)
  }, numeric(1))
  p_q025 <- vapply(eta_marg, function(m) {
    mt <- INLA::inla.tmarginal(plogis, m)
    INLA::inla.qmarginal(0.025, mt)
  }, numeric(1))
  p_q975 <- vapply(eta_marg, function(m) {
    mt <- INLA::inla.tmarginal(plogis, m)
    INLA::inla.qmarginal(0.975, mt)
  }, numeric(1))

  list(
    eta_pred_mean  = eta_mean,
    eta_pred_q025  = eta_q025,
    eta_pred_q975  = eta_q975,
    p_pred_mean    = p_mean,
    p_pred_q025    = p_q025,
    p_pred_q975    = p_q975,
    eta_pred_draws = NULL,   # INLA: marginals, not draws
    elapsed_sec    = elapsed,
    method         = "spde_binomial",
    inla_fit       = fit
  )
}
