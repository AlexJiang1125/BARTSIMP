# =============================================================================
# BARTSIMP-PG wrapper exposing the common competing-methods interface.
# Same signature as fit_bart_logistic / fit_spde_binomial / fit_spde0_binomial.
# =============================================================================

# Caller is expected to have already sourced dev/bartsimp_pg.R (which sources
# dev/rbart.R) -- this wrapper just calls run_sampler_bartsimp_pg with the
# common-interface argument names and reshapes the output.

fit_bartsimp_pg <- function(X_train, locs_train, n_i, k_i,
                            X_pred,  locs_pred,
                            n_iter  = 1500, burn = 500,
                            m_trees = 100L, k_bart = 3,
                            mesh_max_edge = c(0.08, 0.25),
                            mesh_cutoff   = 0.03,
                            bd_mode  = "conditional",
                            adapt_psi = TRUE,
                            verbose  = TRUE) {

  t0 <- Sys.time()
  fit <- run_sampler_bartsimp_pg(
    X = X_train, locs = locs_train, n_i = n_i, k_i = k_i,
    X_pred = X_pred, locs_pred = locs_pred,
    n_iter = n_iter, burn = burn,
    m_trees = m_trees, k_bart = k_bart,
    mesh_max_edge = mesh_max_edge, mesh_cutoff = mesh_cutoff,
    bd_mode = bd_mode, adapt_psi = adapt_psi,
    verbose_every = if (verbose) max(n_iter %/% 6L, 50L) else 0L,
    label = "bartsimp_pg"
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  keep <- (fit$burn + 1):fit$n_iter
  eta_draws <- fit$eta_pred_draws[keep, , drop = FALSE]
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
    method         = "bartsimp_pg",
    fit            = fit
  )
}
