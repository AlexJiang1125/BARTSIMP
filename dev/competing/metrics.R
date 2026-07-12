# =============================================================================
# V3 metrics (paper §4.3.3): RMSE, Brier, log-loss, AIL/ACR/AIS at alpha=0.05.
#
# Inputs are the standardized output schema from the competing wrappers,
# plus ground truth eta_true and p_true at the pred locations.
# Optionally also: observed (k_obs, n_obs) at the pred locations -- used for
# Brier and log-loss, and posterior-predictive ACR(k).
# =============================================================================

ais_interval <- function(true, q025, q975, alpha = 0.05) {
  width <- q975 - q025
  mean(width
       + (2 / alpha) * (q025 - true) * (true < q025)
       + (2 / alpha) * (true - q975) * (true > q975))
}

v3_metrics <- function(fit, eta_true, p_true,
                       k_obs = NULL, n_obs = NULL,
                       alpha = 0.05) {
  stopifnot(length(eta_true) == length(fit$eta_pred_mean),
            length(p_true)   == length(fit$p_pred_mean))

  rmse_eta <- sqrt(mean((fit$eta_pred_mean - eta_true)^2))
  rmse_p   <- sqrt(mean((fit$p_pred_mean   - p_true)^2))

  # eta-scale credible interval metrics
  acr_eta <- mean(eta_true >= fit$eta_pred_q025 & eta_true <= fit$eta_pred_q975)
  ail_eta <- mean(fit$eta_pred_q975 - fit$eta_pred_q025)
  ais_eta <- ais_interval(eta_true, fit$eta_pred_q025, fit$eta_pred_q975, alpha)

  # p-scale credible interval metrics
  acr_p   <- mean(p_true >= fit$p_pred_q025 & p_true <= fit$p_pred_q975)
  ail_p   <- mean(fit$p_pred_q975 - fit$p_pred_q025)

  out <- list(
    method   = fit$method,
    elapsed_sec = fit$elapsed_sec,
    rmse_eta = rmse_eta,
    rmse_p   = rmse_p,
    acr_eta  = acr_eta,
    ail_eta  = ail_eta,
    ais_eta  = ais_eta,
    acr_p    = acr_p,
    ail_p    = ail_p
  )

  # Brier + log-loss (per-trial scoring rules) IF held-out trials are supplied
  if (!is.null(k_obs) && !is.null(n_obs)) {
    stopifnot(length(k_obs) == length(p_true),
              length(n_obs) == length(p_true))
    # Brier: per-trial expected squared error.
    # For each cluster with n_i trials and k_i successes:
    #   Brier_i = (1/n_i) * (k_i * (1 - p_hat_i)^2 + (n_i - k_i) * p_hat_i^2)
    p_hat <- fit$p_pred_mean
    brier <- mean( (k_obs * (1 - p_hat)^2 + (n_obs - k_obs) * p_hat^2) / n_obs )
    # log-loss: clamp p to avoid log(0)
    p_clamp <- pmin(pmax(p_hat, 1e-15), 1 - 1e-15)
    log_loss <- -mean( (k_obs * log(p_clamp)
                       + (n_obs - k_obs) * log(1 - p_clamp)) / n_obs )
    out$brier    <- brier
    out$log_loss <- log_loss
  }
  out
}

# Print a single-row summary
v3_print_row <- function(m) {
  cat(sprintf("  %-15s  RMSE(eta)=%.3f  RMSE(p)=%.3f  ACR=%.3f  AIL=%.3f  AIS=%.2f",
              m$method, m$rmse_eta, m$rmse_p,
              m$acr_eta, m$ail_eta, m$ais_eta))
  if (!is.null(m$brier))    cat(sprintf("  Brier=%.4f", m$brier))
  if (!is.null(m$log_loss)) cat(sprintf("  log_loss=%.4f", m$log_loss))
  cat(sprintf("  t=%.1fs\n", m$elapsed_sec))
}
