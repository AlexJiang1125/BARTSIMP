# =============================================================================
# Verify a native pure-R `make_Qpsi_native` matches INLA::inla.spde.precision,
# and benchmark the two against each other.
#
# Goal: eliminate the one INLA call sitting inside the MCMC hot loop
# (current sampler invokes INLA::inla.spde.precision 2x per iteration).
# The setup-time calls (mesh, spde, A) stay -- they're one-shot per fit.
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
})

# ---- build an spde object on a representative mesh -------------------------
set.seed(2026)
n_cl <- 250
locs <- cbind(runif(n_cl), runif(n_cl))
mesh <- INLA::inla.mesh.2d(loc = locs,
                           max.edge = c(0.08, 0.25),
                           cutoff   = 0.03)
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = c(0.1, 0.5),
                                  prior.sigma = c(1.0, 0.5))
cat(sprintf("Mesh n = %d  (sparse Q dimension %d x %d)\n",
            mesh$n, mesh$n, mesh$n))

# spde$param.inla stores the precomputed sparse blocks
M0 <- as(spde$param.inla$M0, "CsparseMatrix")
M1 <- as(spde$param.inla$M1, "CsparseMatrix")
M2 <- as(spde$param.inla$M2, "CsparseMatrix")
cat(sprintf("nnz: M0=%d  M1=%d  M2=%d\n",
            length(M0@x), length(M1@x), length(M2@x)))

# ---- native implementation -------------------------------------------------
# Lindgren et al. 2011, Matern smoothness nu=1 in dimension d=2:
#   kappa = sqrt(8 nu) / rho
#   sigma^2 = Gamma(nu) / (Gamma(nu+d/2) (4 pi)^(d/2) kappa^(2 nu) tau^2)
#           = 1 / (4 pi kappa^2 tau^2)             (for nu=1, d=2)
#   tau^2 = 1 / (4 pi kappa^2 sigma^2)
# and Q = tau^2 ( kappa^4 M0 + 2 kappa^2 M1 + M2 )
make_Qpsi_native <- function(M0, M1, M2, sigma_m, rho, nu = 1, d = 2) {
  kappa <- sqrt(8 * nu) / rho
  tau2  <- 1 / (4 * pi * kappa^2 * sigma_m^2)
  tau2 * (kappa^4 * M0 + 2 * kappa^2 * M1 + M2)
}

# ---- numerical match against INLA ------------------------------------------
test_grid <- expand.grid(
  sigma_m = c(0.1, 0.5, 1.0, 2.0),
  rho     = c(0.05, 0.2, 0.5, 1.0, 2.0)
)
cat("\n=== Q_native vs Q_inla, frobenius rel error ===\n")
errs <- numeric(nrow(test_grid))
for (i in seq_len(nrow(test_grid))) {
  sg <- test_grid$sigma_m[i]; rh <- test_grid$rho[i]
  Q_inla   <- INLA::inla.spde.precision(spde, theta = c(log(rh), log(sg)))
  Q_native <- make_Qpsi_native(M0, M1, M2, sg, rh)
  d_norm <- sqrt(sum((Q_inla - Q_native)@x^2))
  q_norm <- sqrt(sum(Q_inla@x^2))
  errs[i] <- d_norm / q_norm
  cat(sprintf("  sigma=%.2f rho=%.2f  rel||Q_inla-Q_native||/||Q_inla|| = %.3e\n",
              sg, rh, errs[i]))
}
stopifnot(max(errs) < 1e-10)
cat("\nAll Q's match INLA to <1e-10.\n")

# ---- benchmark --------------------------------------------------------------
N_BENCH <- 1000L
sg <- 0.5; rh <- 0.3
t1 <- system.time(for (k in seq_len(N_BENCH))
  Q1 <- INLA::inla.spde.precision(spde, theta = c(log(rh), log(sg))))
t2 <- system.time(for (k in seq_len(N_BENCH))
  Q2 <- make_Qpsi_native(M0, M1, M2, sg, rh))

cat(sprintf("\n=== %d calls (mesh n=%d) ===\n", N_BENCH, mesh$n))
cat(sprintf("  INLA::inla.spde.precision : %.3f s  (%.1f us / call)\n",
            t1["elapsed"], 1e6 * t1["elapsed"] / N_BENCH))
cat(sprintf("  make_Qpsi_native          : %.3f s  (%.1f us / call)\n",
            t2["elapsed"], 1e6 * t2["elapsed"] / N_BENCH))
cat(sprintf("  speedup                   : %.1fx\n",
            t1["elapsed"] / t2["elapsed"]))

# ---- impact on a realistic MCMC: 2 calls per iter, 1500 iter ---------------
n_iter <- 1500L
cat(sprintf("\nProjected savings for 1 chain x %d iters (2 Q calls/iter):\n", n_iter))
cat(sprintf("  INLA path  : %.2f s\n", 2 * n_iter * t1["elapsed"] / N_BENCH))
cat(sprintf("  native     : %.2f s\n", 2 * n_iter * t2["elapsed"] / N_BENCH))
cat(sprintf("  saved      : %.2f s\n",
            2 * n_iter * (t1["elapsed"] - t2["elapsed"]) / N_BENCH))
