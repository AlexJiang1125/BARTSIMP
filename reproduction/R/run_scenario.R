# =============================================================================
# Section 5 scenario runner
# =============================================================================
# Simulates one dataset for a (omega, surface) scenario, fits the requested
# methods, and returns a tidy data frame of the four performance criteria
# (RMSE, AIL, ACR, AIS) evaluated over the grid G against the true field.
#
# Sourcing this file pulls in all of the Section 5 building blocks, so a driver
# script only needs:  source("reproduction/R/run_scenario.R").
# =============================================================================

# Resolve the directory holding the Section 5 building blocks. Honor an explicit
# `.repro_R` set by the caller (the driver and vignette do this); otherwise try
# the Rscript --file path, then any active source() frame, then a getwd() guess.
.repro_R <- if (exists(".repro_R", inherits = FALSE)) .repro_R else {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  from_source <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { from_source <- of; break }
  }
  if (length(f))               dirname(normalizePath(f))
  else if (!is.null(from_source)) dirname(normalizePath(from_source))
  else if (dir.exists("reproduction/R")) normalizePath("reproduction/R")
  else if (dir.exists("R"))    normalizePath("R")
  else file.path(getwd(), "reproduction", "R")
}
source(file.path(.repro_R, "dgp_section5.R"))
source(file.path(.repro_R, "metrics.R"))
source(file.path(.repro_R, "competitors.R"))
source(file.path(.repro_R, "bartsimp_fit.R"))

# Fit one method and return the four metrics as a named numeric vector.
# `fitter` returns list(point, lower, upper); metrics are vs grid$f_true.
# A fitter failure (e.g. an occasional INLA non-convergence on a sparse
# replicate) is caught and reported as NA metrics so a long sweep is not lost.
.eval_method <- function(fitter, obs, grid, alpha, args = list()) {
  fit <- tryCatch(do.call(fitter, c(list(obs, grid), args)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    warning(sprintf("fitter failed: %s", conditionMessage(fit)), call. = FALSE)
    return(c(RMSE = NA_real_, AIL = NA_real_, ACR = NA_real_, AIS = NA_real_))
  }
  interval_metrics(fit$point, fit$lower, fit$upper, grid$f_true, alpha = alpha)
}

# Run one scenario across the requested methods.
#
#   omega, surface  : scenario definition (see simulate_section5()).
#   seed            : replicate seed (redraws covariate + spatial surfaces).
#   methods         : any subset of c("BARTSIMP","BART","SPDE","SPDE0").
#   n_side, n_clusters : grid + sampling design (paper: 50, 250).
#   *_args          : extra args forwarded to each fitter (e.g. ndpost).
#
# Returns a data frame: replicate, omega, surface, method, RMSE, AIL, ACR, AIS.
run_scenario <- function(omega, surface, seed,
                         methods = c("BARTSIMP", "BART", "SPDE", "SPDE0"),
                         n_side = 50L, n_clusters = 250L,
                         alpha = 0.05,
                         bartsimp_args = list(), bart_args = list(),
                         spde_args = list(), spde0_args = list(),
                         verbose = TRUE) {
  methods <- match.arg(methods, several.ok = TRUE)
  sim <- simulate_section5(omega = omega, surface = surface,
                           n_side = n_side, n_clusters = n_clusters, seed = seed)
  obs <- sim$obs; grid <- sim$grid

  fitters <- list(
    BARTSIMP = list(fn = fit_bartsimp, args = c(bartsimp_args, list(alpha = alpha))),
    BART     = list(fn = fit_bart,     args = c(bart_args,     list(alpha = alpha))),
    SPDE     = list(fn = fit_spde,     args = c(spde_args,     list(alpha = alpha))),
    SPDE0    = list(fn = fit_spde0,    args = c(spde0_args,    list(alpha = alpha)))
  )

  rows <- lapply(methods, function(mn) {
    if (verbose) message(sprintf("  [%s] omega=%.2f surface=%s seed=%s ...",
                                 mn, omega, surface, seed))
    t0 <- Sys.time()
    m <- .eval_method(fitters[[mn]]$fn, obs, grid, alpha, fitters[[mn]]$args)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.frame(replicate = seed, omega = omega, surface = surface, method = mn,
               RMSE = m["RMSE"], AIL = m["AIL"], ACR = m["ACR"], AIS = m["AIS"],
               seconds = secs, row.names = NULL)
  })
  do.call(rbind, rows)
}

# Convenience: run a grid of scenarios x replicates and stack the results.
#   omegas    : vector of covariate-signal shares (paper: c(1,.8,.5,.2,0)).
#   surfaces  : covariate surfaces to sweep.
#   reps      : integer replicate seeds (e.g. 1:20).
run_section5 <- function(omegas = c(1, 0.8, 0.5, 0.2, 0),
                         surfaces = c("tree", "linear", "smooth"),
                         reps = 1:5,
                         methods = c("BARTSIMP", "BART", "SPDE", "SPDE0"),
                         n_side = 50L, n_clusters = 250L, alpha = 0.05,
                         bartsimp_args = list(), bart_args = list(),
                         spde_args = list(), spde0_args = list(),
                         verbose = TRUE) {
  grid_scen <- expand.grid(seed = reps, omega = omegas, surface = surfaces,
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(grid_scen)), function(i) {
    sc <- grid_scen[i, ]
    if (verbose) message(sprintf("[scenario %d/%d] omega=%.2f surface=%s rep=%d",
                                 i, nrow(grid_scen), sc$omega, sc$surface, sc$seed))
    run_scenario(omega = sc$omega, surface = sc$surface, seed = sc$seed,
                 methods = methods, n_side = n_side, n_clusters = n_clusters,
                 alpha = alpha, bartsimp_args = bartsimp_args, bart_args = bart_args,
                 spde_args = spde_args, spde0_args = spde0_args, verbose = verbose)
  })
  do.call(rbind, res)
}

# Aggregate replicate-level results to scenario means (paper-style table).
summarize_section5 <- function(results) {
  agg <- aggregate(cbind(RMSE, AIL, ACR, AIS) ~ omega + surface + method,
                   data = results, FUN = mean)
  agg[order(agg$surface, -agg$omega, agg$method), ]
}
