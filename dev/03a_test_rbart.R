# =============================================================================
# Smoke test for dev/rbart.R
# =============================================================================
# Plain heteroskedastic Gaussian regression, NO PG, NO spatial.
#   y_i ~ N(f(x_i), sigma_i^2)
#   precision Omega_i = 1/sigma_i^2 (known)
# We use a tree-structured truth so BART has a fighting chance to find it:
#   f(x) =  3   if x1 < 0.5
#          -2   if x1 >= 0.5 and x2 < 0.5
#           0   if x1 >= 0.5 and x2 >= 0.5
# If sum-of-trees doesn't recover this surface, the BART module has bugs.
# =============================================================================

source("dev/rbart.R")
set.seed(11)

# Simulate
n  <- 400
x1 <- runif(n); x2 <- runif(n)
X  <- cbind(x1, x2)
f_true <- ifelse(x1 < 0.5, 3,
          ifelse(x2 < 0.5, -2, 0))
sigma  <- 0.5
y      <- f_true + rnorm(n, 0, sigma)
Omega  <- rep(1 / sigma^2, n)

# Prior: BART default-ish. tau chosen so 2*tau*sqrt(m) ~ (y_max - y_min) / 2
m <- 50
tau <- ((max(y) - min(y)) / (2 * 2 * sqrt(m)))  # k = 2
tau2 <- tau^2
cat("tau =", round(tau, 3), "\n")

# Initialize
ensemble <- ensemble_init(m, init_val = mean(y))

# Run
n_iter <- 600
burn   <- 200
fits   <- matrix(NA_real_, n_iter, n)
n_leaves <- numeric(n_iter)

t0 <- Sys.time()
for (it in seq_len(n_iter)) {
  out <- ensemble_sweep(ensemble, X, y, Omega, tau2)
  ensemble <- out$ensemble
  fits[it, ] <- out$fit
  n_leaves[it] <- sum(vapply(ensemble, function(t) sum(t$is_leaf), integer(1)))
  if (it %% 100 == 0) {
    rmse <- sqrt(mean((out$fit - f_true)^2))
    cat(sprintf("iter %4d  acc=%3d/%d  total_leaves=%4d  RMSE(f)=%.3f\n",
                it, out$accepts, m, n_leaves[it], rmse))
  }
}
t1 <- Sys.time()
cat("Elapsed:", format(t1 - t0), "\n")

# Posterior summaries
keep <- (burn + 1):n_iter
f_post <- colMeans(fits[keep, ])
rmse <- sqrt(mean((f_post - f_true)^2))
cat("\nPosterior mean RMSE vs truth:", round(rmse, 3), "\n")
cat("cor(f_post, f_true):", round(cor(f_post, f_true), 3), "\n")

# Coverage on a coarse grid
gridx <- seq(0.05, 0.95, by = 0.1)
grid <- expand.grid(x1 = gridx, x2 = gridx)
Xgrid <- as.matrix(grid)
fgrid_true <- ifelse(Xgrid[,1] < 0.5, 3,
              ifelse(Xgrid[,2] < 0.5, -2, 0))
# predict ensemble samples at grid
fgrid <- matrix(NA_real_, length(keep), nrow(Xgrid))
for (i in seq_along(keep)) {
  # we don't store ensembles per draw to save mem; recompute from final
  # ensemble would not give samples. Skip credible-band check, use point estimate only.
  fgrid[i, ] <- NA
}
# Instead just save current-state predictions at grid + look at point estimate
fgrid_post <- ensemble_predict(ensemble, Xgrid)
cat("Grid RMSE (final state):", round(sqrt(mean((fgrid_post - fgrid_true)^2)), 3), "\n")

# Diagnostic plot
pdf("dev/03a_test_rbart_diagnostics.pdf", width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(f_true, f_post, pch = 20, col = rgb(0,0,0,0.4),
     xlab = "f true", ylab = "f post mean", main = "BART recovery (training)")
abline(0, 1, col = "red", lty = 2)
plot(n_leaves, type = "l", xlab = "iter", ylab = "total leaves",
     main = "Tree complexity")
plot(rowMeans((fits - matrix(f_true, n_iter, n, byrow = TRUE))^2)^0.5,
     type = "l", xlab = "iter", ylab = "RMSE vs truth",
     main = "RMSE trace")
plot(Xgrid, col = ifelse(fgrid_post > 1.5, "red",
                  ifelse(fgrid_post < -1,  "blue", "gray")),
     pch = 15, cex = 2, xlab = "x1", ylab = "x2",
     main = "Final predicted surface")
abline(v = 0.5, lty = 3); abline(h = 0.5, lty = 3)
par(op); dev.off()
cat("Diagnostics saved to dev/03a_test_rbart_diagnostics.pdf\n")
