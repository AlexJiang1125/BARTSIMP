# =============================================================================
# BART-logistic (no spatial component)
#
# Fits BART::lbart on per-trial binary outcomes obtained by expanding the
# cluster-level (k_i, n_i) counts into n_i Bernoulli trials per cluster.
# Then predicts f(x) at the supplied grid (X_pred). locs_pred is ignored
# (BART-logistic has no spatial component).
# =============================================================================

# BART package is loaded lazily at call time so that sourcing this file is
# harmless even when the BART CRAN package isn't installed (the runner can
# then dispatch to the other 3 methods without aborting).

fit_bart_logistic <- function(X_train, locs_train, n_i, k_i,
                              X_pred,  locs_pred,
                              n_iter = 1500, burn = 500,
                              ntree = 50L,
                              verbose = TRUE) {
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package 'BART' is required for fit_bart_logistic() but not installed. ",
         "Run install.packages('BART') (or the cluster install_deps.R), or ",
         "drop 'bart_logistic' from the methods list in sim_design.R.")
  }
  stopifnot(length(n_i) == nrow(X_train),
            length(k_i) == nrow(X_train))

  # Expand to per-trial binary outcomes
  rep_idx <- rep(seq_len(nrow(X_train)), times = n_i)
  X_expanded <- X_train[rep_idx, , drop = FALSE]
  y_expanded <- unlist(lapply(seq_along(n_i), function(i) {
    c(rep(1L, k_i[i]), rep(0L, n_i[i] - k_i[i]))
  }))
  stopifnot(length(y_expanded) == sum(n_i),
            nrow(X_expanded) == sum(n_i))

  if (verbose)
    cat(sprintf("[bart_logistic] n_trial=%d (n_cl=%d), n_pred=%d, ntree=%d, ndpost=%d\n",
                length(y_expanded), nrow(X_train), nrow(X_pred), ntree, n_iter - burn))

  t0 <- Sys.time()
  fit <- BART::lbart(
    x.train = X_expanded, y.train = y_expanded,
    x.test  = X_pred,
    ntree   = ntree,
    ndpost  = n_iter - burn,
    nskip   = burn,
    printevery = if (verbose) max((n_iter %/% 6L), 50L) else n_iter + 1L
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # yhat.test is matrix [ndpost, n_pred] of f(x) on the LATENT (logit) scale
  eta_draws <- fit$yhat.test
  p_draws   <- plogis(eta_draws)

  list(
    eta_pred_mean  = colMeans(eta_draws),
    eta_pred_q025  = apply(eta_draws, 2, quantile, probs = 0.025),
    eta_pred_q975  = apply(eta_draws, 2, quantile, probs = 0.975),
    p_pred_mean    = colMeans(p_draws),
    p_pred_q025    = apply(p_draws, 2, quantile, probs = 0.025),
    p_pred_q975    = apply(p_draws, 2, quantile, probs = 0.975),
    eta_pred_draws = eta_draws,
    elapsed_sec    = elapsed,
    method         = "bart_logistic",
    n_iter         = n_iter, burn = burn, ntree = ntree
  )
}
