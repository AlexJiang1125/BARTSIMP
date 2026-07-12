# =============================================================================
# Section 6 cross-validation exercise (synthetic Kenya-DHS-like data)
# =============================================================================
# Reproduces the *structure* of the paper's CV comparison (Section 6.1): split
# the clusters (enumeration areas) into train/test with a stratified design,
# fit BARTSIMP / BART / SPDE / SPDE0 on the training children, predict at the
# held-out cluster locations, and score with RMSE / AIL / ACR / AIS.
#
# EVALUATION TARGET. We predict the mean WHZ surface and score it against the
# held-out CLUSTER-MEAN WHZ (the average of the test children in each cluster).
# Averaging over a cluster's children removes most of the child-level noise, so
# the mean-field intervals returned by the fitters are the right object to score
# (the residual noise in a cluster mean is sigma_e^2 / n_j, small for n_j ~ 20).
# The paper instead scores per-child with fully predictive intervals; that only
# changes interval width by a roughly constant sigma_e inflation and does not
# change the ranking of the methods, which is what this exercise demonstrates.
#
# Sourcing run_scenario.R already pulls in metrics/competitors/bartsimp_fit;
# this file additionally needs dgp_section6.R and areal.R.
# =============================================================================

# Build the per-test-cluster target table: one row per held-out cluster with its
# location, covariate values, and the mean WHZ of its test children.
.make_test_targets <- function(obs_test, cov_names, id_col = "cluster") {
  agg_y <- tapply(obs_test$y, obs_test[[id_col]], mean)
  first <- !duplicated(obs_test[[id_col]])
  base  <- obs_test[first, c(id_col, "s1", "s2", cov_names), drop = FALSE]
  base  <- base[match(names(agg_y), base[[id_col]]), , drop = FALSE]
  base$y <- as.numeric(agg_y)
  rownames(base) <- NULL
  base
}

# One CV fold: returns a tidy data frame of metrics for the requested methods.
#
#   sim      : output of simulate_section6().
#   methods  : subset of c("BARTSIMP","BART","SPDE","SPDE0").
#   seed     : split seed (also seeds the fitters).
#   *_args   : extra args forwarded to each fitter.
run_cv_fold_section6 <- function(sim,
                                 methods = c("BARTSIMP", "BART", "SPDE", "SPDE0"),
                                 prop_train = 0.8, alpha = 0.05, seed = 1L,
                                 cov_names = SECTION6_COVARIATES,
                                 bartsimp_args = list(), bart_args = list(),
                                 spde_args = list(), spde0_args = list(),
                                 verbose = FALSE) {
  methods <- match.arg(methods, several.ok = TRUE)
  sp <- stratified_cluster_split(sim$clusters, prop_train = prop_train, seed = seed)
  obs_tr <- subset_by_cluster(sim$obs, sp$train)
  obs_te <- subset_by_cluster(sim$obs, sp$test)
  target <- .make_test_targets(obs_te, cov_names)   # has $y = cluster-mean WHZ

  fitters <- list(
    BARTSIMP = list(fn = fit_bartsimp,
                    args = c(bartsimp_args, list(alpha = alpha, cov_names = cov_names))),
    BART     = list(fn = fit_bart,
                    args = c(bart_args,     list(alpha = alpha, cov_names = cov_names))),
    SPDE     = list(fn = fit_spde,
                    args = c(spde_args,     list(alpha = alpha, cov_names = cov_names))),
    SPDE0    = list(fn = fit_spde0,
                    args = c(spde0_args,    list(alpha = alpha)))
  )

  rows <- lapply(methods, function(mn) {
    if (verbose) message(sprintf("  [%s] fold seed=%s ...", mn, seed))
    t0 <- Sys.time()
    fit <- tryCatch(do.call(fitters[[mn]]$fn, c(list(obs_tr, target), fitters[[mn]]$args)),
                    error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(fit, "error")) {
      warning(sprintf("[%s] failed: %s", mn, conditionMessage(fit)), call. = FALSE)
      m <- c(RMSE = NA, AIL = NA, ACR = NA, AIS = NA)
    } else {
      m <- interval_metrics(fit$point, fit$lower, fit$upper, target$y, alpha = alpha)
    }
    data.frame(fold = seed, method = mn,
               RMSE = m["RMSE"], AIL = m["AIL"], ACR = m["ACR"], AIS = m["AIS"],
               n_train = length(sp$train), n_test = length(sp$test),
               seconds = secs, row.names = NULL)
  })
  do.call(rbind, rows)
}

# Repeat the CV fold over several split seeds and stack the results (paper: 10
# repeats). Each call re-splits the same simulated dataset.
run_cv_section6 <- function(sim, reps = 1:5,
                            methods = c("BARTSIMP", "BART", "SPDE", "SPDE0"),
                            prop_train = 0.8, alpha = 0.05,
                            cov_names = SECTION6_COVARIATES,
                            bartsimp_args = list(), bart_args = list(),
                            spde_args = list(), spde0_args = list(),
                            verbose = TRUE) {
  res <- lapply(reps, function(s) {
    if (verbose) message(sprintf("[CV fold %d]", s))
    run_cv_fold_section6(sim, methods = methods, prop_train = prop_train,
                         alpha = alpha, seed = s, cov_names = cov_names,
                         bartsimp_args = bartsimp_args, bart_args = bart_args,
                         spde_args = spde_args, spde0_args = spde0_args,
                         verbose = verbose)
  })
  do.call(rbind, res)
}

# Aggregate CV folds to per-method means (+ SD), paper-table layout.
summarize_cv_section6 <- function(results) {
  agg_mean <- aggregate(cbind(RMSE, AIL, ACR, AIS) ~ method, data = results, FUN = mean)
  agg_sd   <- aggregate(cbind(RMSE, AIL, ACR, AIS) ~ method, data = results, FUN = sd)
  names(agg_sd)[-1] <- paste0(names(agg_sd)[-1], "_sd")
  merge(agg_mean, agg_sd, by = "method")
}
