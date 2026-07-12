# =============================================================================
# Section 5 performance criteria (paper Section 5.2)
# =============================================================================
# All metrics are evaluated over the gridded surface G (|G| cells) against the
# true field f(x_g). A method supplies, per grid cell, a point prediction and a
# (1 - alpha) prediction interval [L, U]. The paper uses alpha = 0.05.
#
#   RMSE  point-prediction accuracy.
#         The paper writes  (1/|G|) sum (fhat - f)^2  (i.e. a mean squared
#         error). We report the conventional root-mean-squared error
#         sqrt(mean((fhat - f)^2)); set `sqrt = FALSE` to match the paper's
#         written formula exactly.
#   AIL   average interval length,  mean(U - L).
#   ACR   average coverage rate,    mean( f in [L, U] ).   Nominal = 1 - alpha.
#   AIS   average interval score (Gneiting & Raftery, 2007), lower is better:
#         (U - L) + (2/alpha)(L - f) 1{f < L} + (2/alpha)(f - U) 1{f > U}.
# =============================================================================

rmse <- function(point, f_true, sqrt = TRUE) {
  v <- mean((point - f_true)^2)
  if (sqrt) base::sqrt(v) else v
}

# Interval metrics from explicit lower/upper bounds.
interval_metrics <- function(point, lower, upper, f_true, alpha = 0.05,
                             rmse_sqrt = TRUE) {
  stopifnot(length(point) == length(f_true),
            length(lower) == length(f_true),
            length(upper) == length(f_true))
  ail <- mean(upper - lower)
  acr <- mean(f_true >= lower & f_true <= upper)
  ais <- mean((upper - lower) +
              (2 / alpha) * (lower - f_true) * (f_true < lower) +
              (2 / alpha) * (f_true - upper) * (f_true > upper))
  c(RMSE = rmse(point, f_true, sqrt = rmse_sqrt),
    AIL = ail, ACR = acr, AIS = ais)
}

# Convenience wrapper: compute all four metrics from a matrix of posterior
# draws (rows = draws, cols = grid cells). Point = posterior mean; interval =
# central (1 - alpha) posterior quantiles.
metrics_from_draws <- function(draws, f_true, alpha = 0.05, rmse_sqrt = TRUE) {
  draws <- as.matrix(draws)
  stopifnot(ncol(draws) == length(f_true))
  point <- colMeans(draws)
  lower <- apply(draws, 2, quantile, probs = alpha / 2,     names = FALSE)
  upper <- apply(draws, 2, quantile, probs = 1 - alpha / 2, names = FALSE)
  interval_metrics(point, lower, upper, f_true, alpha = alpha,
                   rmse_sqrt = rmse_sqrt)
}
