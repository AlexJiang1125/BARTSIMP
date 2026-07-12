# =============================================================================
# Section 5 competitor methods
# =============================================================================
# Each fitter trains on the 250-cluster observation data and predicts the true
# field f(x_g) at all |G| grid cells, returning a common structure:
#
#   list(point = <length |G|>, lower = <length |G|>, upper = <length |G|>)
#
# so results feed directly into interval_metrics() (see metrics.R). Intervals
# are central (1 - alpha) intervals; the prediction target is the mean field
# (no observation noise), consistent with how BARTSIMP predicts f.
#
# Methods (paper Section 5.1):
#   BART   : standard BART on covariates only, no spatial field (dbarts).
#   SPDE   : Matern-SPDE GMRF spatial field + linear covariate effects (INLA).
#   SPDE0  : Matern-SPDE GMRF spatial field, intercept only, no covariates.
#
# BARTSIMP itself is fit in the vignette/driver via bartsimp() +
# bartsimp_predict(); it is not included here to keep this file package-free
# except for dbarts/INLA.
# =============================================================================

# ---- BART (covariate-only) ---------------------------------------------------
# cov_names: covariate columns to use (default the Section 5 pair c("x1","x2");
# Section 6 passes its six geospatial covariates).
fit_bart <- function(obs, grid, ntree = 20L, ndpost = 2000L, nskip = 2000L,
                     alpha = 0.05, seed = 1L, cov_names = c("x1", "x2")) {
  if (!requireNamespace("dbarts", quietly = TRUE))
    stop("Package 'dbarts' is required for fit_bart().")
  set.seed(seed)
  fit <- dbarts::bart(
    x.train = as.matrix(obs[, cov_names, drop = FALSE]),
    y.train = obs$y,
    x.test  = as.matrix(grid[, cov_names, drop = FALSE]),
    ntree = ntree, ndpost = ndpost, nskip = nskip,
    keeptrees = FALSE, verbose = FALSE
  )
  draws <- fit$yhat.test                      # ndpost x |G| posterior of f(x)
  list(
    point = colMeans(draws),
    lower = apply(draws, 2, quantile, probs = alpha / 2,     names = FALSE),
    upper = apply(draws, 2, quantile, probs = 1 - alpha / 2, names = FALSE)
  )
}

# ---- SPDE family (shared engine) --------------------------------------------
# covariates = TRUE  -> SPDE  (intercept + linear covariates + spatial field)
# covariates = FALSE -> SPDE0 (intercept + spatial field)
# cov_names: which covariate columns enter the linear predictor when
# covariates = TRUE (default the Section 5 pair; Section 6 passes its six).
fit_spde_family <- function(obs, grid, covariates = TRUE, alpha = 0.05,
                            cov_names = c("x1", "x2"),
                            max.edge = c(0.05, 0.2), cutoff = 0.02,
                            prior.range = c(0.3, 0.5), prior.sigma = c(1, 0.5)) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("Package 'INLA' is required for the SPDE methods.")

  loc_train <- as.matrix(obs[, c("s1", "s2")])
  loc_grid  <- as.matrix(grid[, c("s1", "s2")])

  mesh <- INLA::inla.mesh.2d(loc = loc_train, max.edge = max.edge,
                             cutoff = cutoff)
  spde <- INLA::inla.spde2.pcmatern(mesh = mesh,
                                    prior.range = prior.range,
                                    prior.sigma = prior.sigma)

  A_train <- INLA::inla.spde.make.A(mesh, loc = loc_train)
  A_grid  <- INLA::inla.spde.make.A(mesh, loc = loc_grid)
  idx <- seq_len(spde$n.spde)

  # Intercept must be full-length so the identity-mapped effect block matches
  # the response rows in inla.stack (a scalar 1 would give a 1-row effect).
  eff_train <- data.frame(intercept = rep(1, nrow(obs)))
  eff_grid  <- data.frame(intercept = rep(1, nrow(grid)))
  if (covariates) {
    eff_train <- cbind(eff_train, obs[, cov_names, drop = FALSE])
    eff_grid  <- cbind(eff_grid,  grid[, cov_names, drop = FALSE])
  }

  stk_est <- INLA::inla.stack(
    data = list(y = obs$y), A = list(A_train, 1),
    effects = list(i = idx, eff_train), tag = "est")
  stk_pred <- INLA::inla.stack(
    data = list(y = NA), A = list(A_grid, 1),
    effects = list(i = idx, eff_grid), tag = "pred")
  stk <- INLA::inla.stack(stk_est, stk_pred)

  rhs <- if (covariates)
    paste(c("intercept", cov_names), collapse = " + ")
  else
    "intercept"
  form <- stats::as.formula(paste0("y ~ -1 + ", rhs, " + f(i, model = spde)"))

  res <- INLA::inla(
    form, data = INLA::inla.stack.data(stk), family = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
    verbose = FALSE)

  ii <- INLA::inla.stack.index(stk, tag = "pred")$data
  lp <- res$summary.linear.predictor[ii, ]
  qlo <- paste0(alpha / 2, "quant")
  qhi <- paste0(1 - alpha / 2, "quant")
  list(
    point = lp[, "mean"],
    lower = lp[, if (qlo %in% names(lp)) qlo else "0.025quant"],
    upper = lp[, if (qhi %in% names(lp)) qhi else "0.975quant"]
  )
}

fit_spde  <- function(obs, grid, ...) fit_spde_family(obs, grid, covariates = TRUE,  ...)
fit_spde0 <- function(obs, grid, ...) fit_spde_family(obs, grid, covariates = FALSE, ...)
