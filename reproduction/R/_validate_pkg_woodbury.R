# =============================================================================
# Validate the integrated method = "woodbury" branch of BARTSIMP::bartsimp_predict
# against the default method = "inla" path, using devtools::load_all() so the
# in-tree R source (with the new `method` argument) is exercised directly.
#
# Run from the package root:
#   Rscript reproduction/R/_validate_pkg_woodbury.R
# =============================================================================
suppressWarnings(suppressMessages({
  library(devtools)
  if (requireNamespace("INLA", quietly = TRUE))
    try(INLA::inla.setOption(fmesher.evolution.warn = FALSE), silent = TRUE)
}))

pkg_root <- "/Users/alexziyujiang/Documents/GitHub/BARTSIMP"
message("Loading package from: ", pkg_root)
devtools::load_all(pkg_root, quiet = TRUE)

set.seed(11)

# ---- small synthetic prediction problem --------------------------------------
n_train <- 120L
n_test  <- 200L
ndpost  <- 12L

s_train <- data.frame(s1 = runif(n_train), s2 = runif(n_train))
s_test  <- expand.grid(s1 = seq(0.05, 0.95, length.out = 20),
                       s2 = seq(0.05, 0.95, length.out = 10))
s_test  <- data.frame(s1 = s_test$s1, s2 = s_test$s2)

mesh <- INLA::inla.mesh.2d(loc = as.matrix(s_train),
                           max.edge = c(0.1, 0.4), cutoff = 0.02)

y_train <- rnorm(n_train)
tree_draws_train <- matrix(rnorm(ndpost * n_train), nrow = ndpost)
tree_draws_test  <- matrix(rnorm(ndpost * n_test),  nrow = ndpost)

kappas    <- runif(ndpost, 1.5, 3.5)
sigmams   <- runif(ndpost, 0.4, 0.9)
sigma_obs <- runif(ndpost, 0.8, 1.2)

# ---- INLA path (reference) ---------------------------------------------------
t0 <- Sys.time()
res_inla <- bartsimp_predict(
  y_train = y_train,
  tree_draws_train = tree_draws_train,
  tree_draws_test = tree_draws_test,
  s_train = s_train, s_test = s_test,
  kappas = kappas, sigmams = sigmams, sigma_obs = sigma_obs,
  mesh = mesh, method = "inla", verbose = FALSE
)
t_inla <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Woodbury path (fast) ----------------------------------------------------
t0 <- Sys.time()
res_wb <- bartsimp_predict(
  y_train = y_train,
  tree_draws_train = tree_draws_train,
  tree_draws_test = tree_draws_test,
  s_train = s_train, s_test = s_test,
  kappas = kappas, sigmams = sigmams, sigma_obs = sigma_obs,
  mesh = mesh, method = "woodbury", verbose = FALSE
)
t_wb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- compare -----------------------------------------------------------------
sm_i <- res_inla$spatial_mean_test
sm_w <- res_wb$spatial_mean_test
max_abs <- max(abs(sm_i - sm_w))
rel     <- max_abs / max(1e-12, max(abs(sm_i)))
corr    <- cor(as.numeric(sm_i), as.numeric(sm_w))

cat(sprintf("\n==== integrated woodbury vs inla ====\n"))
cat(sprintf("spatial_mean dims      : %d x %d\n", nrow(sm_i), ncol(sm_i)))
cat(sprintf("max |INLA - Woodbury|  : %.3e\n", max_abs))
cat(sprintf("relative max diff      : %.3e\n", rel))
cat(sprintf("correlation            : %.8f\n", corr))
cat(sprintf("time INLA / Woodbury   : %.2fs / %.2fs  (%.0fx)\n",
            t_inla, t_wb, t_inla / max(t_wb, 1e-6)))

pm_max <- max(abs(res_inla$prediction_mean_test - res_wb$prediction_mean_test))
cat(sprintf("max |pred_mean diff|   : %.3e\n", pm_max))

ok <- (max_abs < 1e-5) && (corr > 0.9999)
cat(sprintf("\nRESULT: %s\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
