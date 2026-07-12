# =============================================================================
# BARTSIMP fitter for the Section 5 simulation
# =============================================================================
# Wraps the two-stage BARTSIMP workflow behind the same contract used by the
# competitor methods in competitors.R:
#
#   list(point = <length |G|>, lower = <length |G|>, upper = <length |G|>)
#
# Stage 1  bartsimp()          : MCMC over the BART tree ensemble and the
#                                Matern-SPDE spatial hyperparameters
#                                (kappa = sqrt(8 nu)/rho, sigma_m, sigma_e).
#                                Returns per-draw tree predictions at train and
#                                grid locations, the hyperparameter draws, and
#                                the internal SPDE mesh.
# Stage 2  bartsimp_predict()  : for each draw, recovers the Matern field from
#                                the tree residuals (hyperparameters fixed at the
#                                drawn kappa/sigma_m) and projects it to the grid,
#                                giving posterior draws of the mean field
#                                f(x_g) = tree(x_g) + w(s_g).
#                                predict_method = "woodbury" (default) uses the
#                                closed-form GMRF posterior mean (sparse Cholesky)
#                                instead of refitting INLA per draw: numerically
#                                identical (~1e-7) but orders of magnitude faster.
#                                predict_method = "inla" restores the per-draw
#                                INLA fit.
#
# ORDERING CONTRACT
#   Both bartsimp() and bartsimp_predict() sort training rows by order(s1, s2)
#   internally. We pre-sort `obs` by the same key so both internal sorts are
#   identity permutations and the returned tree-draw columns stay aligned with
#   the (sorted) response we hand to stage 2. Test/grid rows are never reordered,
#   so output columns align with `grid` row order and hence with grid$f_true.
#
# The prediction target is the mean field (no observation noise), matching the
# competitors and how the paper evaluates over the grid G.
# =============================================================================

fit_bartsimp <- function(obs, grid,
                         ntree = 50L, ndpost = 1000L, nskip = 1L,
                         nwarmup = 1000L, alpha = 0.05, seed = 1L,
                         use_sampler_sigma = TRUE,
                         predict_method = c("woodbury", "inla"),
                         cov_names = c("x1", "x2"),
                         return_draws = FALSE,
                         verbose = FALSE, ...) {
  if (!requireNamespace("BARTSIMP", quietly = TRUE))
    stop("Package 'BARTSIMP' is required for fit_bartsimp().")
  predict_method <- match.arg(predict_method)

  # The Woodbury closed form needs the per-draw observation noise SD, so it must
  # come from the sampler. Fall back to INLA (which can estimate the noise) if
  # the caller has explicitly turned that off.
  if (predict_method == "woodbury" && !use_sampler_sigma) {
    warning("predict_method = 'woodbury' requires the sampler noise SD; ",
            "setting use_sampler_sigma = TRUE.", call. = FALSE)
    use_sampler_sigma <- TRUE
  }

  # Pre-sort so bartsimp()'s and bartsimp_predict()'s internal order(s1, s2)
  # sorts are both no-ops, keeping tree-draw columns aligned with the response.
  ord <- order(obs$s1, obs$s2)
  obs <- obs[ord, , drop = FALSE]

  fit <- BARTSIMP::bartsimp(
    x.train = as.matrix(obs[, cov_names, drop = FALSE]),
    y.train = obs$y,
    s1 = obs$s1, s2 = obs$s2,
    x.test  = as.matrix(grid[, cov_names, drop = FALSE]),
    s1.test = grid$s1, s2.test = grid$s2,
    ntree = ntree, ndpost = ndpost, nskip = nskip,
    nwarmup = nwarmup, seed = seed, ...
  )

  pred <- BARTSIMP::bartsimp_predict(
    y_train          = obs$y,
    tree_draws_train = fit$yhat.train,
    tree_draws_test  = fit$yhat.test,
    s_train          = obs[, c("s1", "s2")],
    s_test           = grid[, c("s1", "s2")],
    kappas           = fit$kappas,
    sigmams          = fit$sigmams,
    sigma_obs        = if (use_sampler_sigma) fit$sigma else NULL,
    mesh             = fit$mesh,
    method           = predict_method,
    verbose          = verbose
  )

  draws <- pred$prediction_mean_test          # ndpost x |G| posterior of f(x_g)
  out <- list(
    point = colMeans(draws),
    lower = apply(draws, 2, quantile, probs = alpha / 2,     names = FALSE),
    upper = apply(draws, 2, quantile, probs = 1 - alpha / 2, names = FALSE),
    fit   = fit                                # retained for diagnostics
  )
  # Areal aggregation (Section 6) needs the full cell-level posterior, since an
  # area's credible interval must come from the areal posterior, not from
  # combining cell-level intervals. Return the draws only when asked.
  if (return_draws) out$draws <- draws
  out
}
