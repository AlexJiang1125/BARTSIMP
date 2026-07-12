# Smoke test for dgp_section5.R, metrics.R, competitors.R
# Run with framework R: Rscript reproduction/R/_smoke_test.R
# Locate this script's own directory so it runs from anywhere.
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
here <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(here, "dgp_section5.R"))
source(file.path(here, "metrics.R"))
source(file.path(here, "competitors.R"))

set.seed(2025)
# small grid for speed; still exercises the full |G|-length contract
sim <- simulate_section5(omega = 0.5, surface = "tree", n_side = 20L,
                         n_clusters = 120L, seed = 7L)
obs  <- sim$obs
grid <- sim$grid
G <- nrow(grid)
cat(sprintf("grid cells |G| = %d ; obs rows = %d\n", G, nrow(obs)))

check_fit <- function(name, fit) {
  ok <- is.list(fit) &&
    all(c("point", "lower", "upper") %in% names(fit)) &&
    length(fit$point) == G &&
    length(fit$lower) == G &&
    length(fit$upper) == G &&
    all(is.finite(fit$point)) &&
    all(fit$lower <= fit$upper)
  cat(sprintf("[%s] length/contract ok = %s\n", name, ok))
  if (!ok) stop(sprintf("%s failed contract", name))
  m <- interval_metrics(fit$point, fit$lower, fit$upper, grid$f_true)
  cat(sprintf("  RMSE=%.3f AIL=%.3f ACR=%.3f AIS=%.3f\n",
              m["RMSE"], m["AIL"], m["ACR"], m["AIS"]))
  invisible(m)
}

cat("\n-- fit_bart --\n")
fb <- fit_bart(obs, grid, ntree = 20L, ndpost = 200L, nskip = 200L, seed = 1L)
check_fit("BART", fb)

cat("\n-- fit_spde --\n")
fs <- fit_spde(obs, grid)
check_fit("SPDE", fs)

cat("\n-- fit_spde0 --\n")
fs0 <- fit_spde0(obs, grid)
check_fit("SPDE0", fs0)

cat("\n-- metrics_from_draws sanity --\n")
# fabricate draws centered on f_true to confirm ACR ~ nominal
draws <- matrix(rnorm(50 * G, mean = rep(grid$f_true, each = 50), sd = 0.5),
                nrow = 50, ncol = G)
md <- metrics_from_draws(draws, grid$f_true)
cat(sprintf("  RMSE=%.3f AIL=%.3f ACR=%.3f AIS=%.3f\n",
            md["RMSE"], md["AIL"], md["ACR"], md["AIS"]))

cat("\nALL SMOKE TESTS PASSED\n")
