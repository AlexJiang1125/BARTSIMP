# =============================================================================
# Validation V5 -- scalability test for marginal-likelihood evaluation
# =============================================================================
# Two ways to evaluate log p(yhat | psi, Omega) under the PG-augmented model
#   yhat_i | z_i ~ N(z_i, Omega_i^{-1})
#   z ~ N(0, Sigma_z),  Sigma_z = sigma_m^2 Matern correlation
#
# (a) PG-Exact   -- form Sigma_y = Sigma_z + diag(1/Omega) (dense n x n),
#                   Cholesky, plug into Gaussian log-density. O(n^3).
# (b) PG-Woodbury -- replace Sigma_z by A Q_psi^{-1} A', use Woodbury identity
#                    on the data precision. One sparse Cholesky of Q_u. Cost
#                    is roughly O(n_mesh^{1.5}) for the sparse solve.
#
# This mirrors the scalability comparison in Section S.6 of the published
# BARTSIMP supplementary material (which compared INLA-SPDE vs dense exact
# under a Gaussian model).
# =============================================================================

source("dev/bartsimp_pg.R")

if (!exists("SMOKE_TEST")) SMOKE_TEST <- TRUE

if (SMOKE_TEST) {
  N_VALS   <- c(100, 200, 500)
  N_REPS   <- 3
  N_CALLS  <- 5
} else {
  N_VALS   <- c(100, 200, 500, 1000, 2000)
  N_REPS   <- 5
  N_CALLS  <- 10
}

cat(sprintf("Config: n_vals=[%s], n_reps=%d, n_calls=%d\n",
            paste(N_VALS, collapse=","), N_REPS, N_CALLS))

# -----------------------------------------------------------------------------
# Dense (PG-Exact) marginal log-likelihood
#   Sigma_y = Sigma_z(sigma_m, rho) + diag(1/Omega)
#   log p(yhat) = -0.5 [n log(2pi) + log|Sigma_y| + yhat' Sigma_y^{-1} yhat]
# -----------------------------------------------------------------------------
margprob_pg_exact <- function(yhat, Omega, D, sigma_m, rho, nu = 1) {
  kappa  <- sqrt(8 * nu) / rho
  S_y    <- sigma_m^2 * matern_corr(D, kappa, nu)
  diag(S_y) <- diag(S_y) + 1 / Omega
  L      <- chol(S_y)                     # dense Cholesky
  logdet <- 2 * sum(log(diag(L)))
  q      <- backsolve(L, yhat, transpose = TRUE)
  quad   <- sum(q * q)
  N      <- length(yhat)
  -0.5 * (N * log(2 * pi) + logdet + quad)
}

# -----------------------------------------------------------------------------
# Setup helper: simulate balanced-scenario data + mesh objects at given n
# -----------------------------------------------------------------------------
make_setup <- function(n_cl, seed) {
  dat <- gen_data(seed, n_cl = n_cl,
                  covariate_signal = TRUE, spatial_signal = TRUE)
  # PG-style yhat: use the true linear predictor for clean comparison
  Omega <- BayesLogit::rpg(num = n_cl, h = dat$n_i, z = dat$eta_true)
  ybar  <- (dat$k_i - dat$n_i / 2) / Omega
  yhat  <- ybar - dat$f_true
  D     <- as.matrix(dist(dat$locs))

  mesh <- INLA::inla.mesh.2d(loc = dat$locs,
                             max.edge = c(0.08, 0.25), cutoff = 0.03)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.1, 0.5),
                                          prior.sigma = c(1.0, 0.5))
  A <- as(INLA::inla.spde.make.A(mesh = mesh, loc = dat$locs),
          "CsparseMatrix")
  list(yhat = yhat, Omega = Omega, D = D, A = A, spde = spde, n_mesh = mesh$n)
}

# -----------------------------------------------------------------------------
# Time-helper: average per-call wall time over n_calls
# -----------------------------------------------------------------------------
time_per_call <- function(fn, n_calls) {
  fn()  # warm-up
  t0 <- Sys.time()
  for (i in seq_len(n_calls)) fn()
  as.numeric(difftime(Sys.time(), t0, units = "secs")) / n_calls
}

# -----------------------------------------------------------------------------
# Run the grid
# -----------------------------------------------------------------------------
out_rows <- list()
for (n_cl in N_VALS) {
  cat(sprintf("\n=== n_cl = %d ===\n", n_cl))
  t_exact <- numeric(N_REPS)
  t_wood  <- numeric(N_REPS)
  mesh_n  <- NA_integer_

  for (rep in seq_len(N_REPS)) {
    s <- make_setup(n_cl, seed = 1000 + rep)
    mesh_n <- s$n_mesh
    f_exact <- function() margprob_pg_exact(s$yhat, s$Omega, s$D,
                                            sigma_m = 0.8, rho = 0.3)
    Q_psi   <- make_Qpsi(s$spde, sigma_m = 0.8, rho = 0.3)
    f_wood  <- function() margprob_woodbury(s$yhat, s$Omega, s$A, Q_psi)

    t_exact[rep] <- time_per_call(f_exact, N_CALLS)
    t_wood[rep]  <- time_per_call(f_wood,  N_CALLS)
  }

  cat(sprintf("  exact   : mean=%.4f s   sd=%.4f s\n",
              mean(t_exact), sd(t_exact)))
  cat(sprintf("  woodbury: mean=%.4f s   sd=%.4f s   (mesh=%d)\n",
              mean(t_wood), sd(t_wood), mesh_n))
  cat(sprintf("  speedup = exact/wood = %.2fx\n",
              mean(t_exact) / mean(t_wood)))

  out_rows[[length(out_rows) + 1]] <- data.frame(
    n_cl       = n_cl,
    n_mesh     = mesh_n,
    exact_mean = mean(t_exact),
    exact_sd   = sd(t_exact),
    wood_mean  = mean(t_wood),
    wood_sd    = sd(t_wood),
    speedup    = mean(t_exact) / mean(t_wood)
  )
}

results <- do.call(rbind, out_rows)
cat("\n=== Summary ===\n")
print(results, row.names = FALSE)

# Save raw + plot
saveRDS(results,
        sprintf("dev/09_scalability_%s.rds", if (SMOKE_TEST) "smoke" else "full"))

pdf(sprintf("dev/09_scalability_%s.pdf", if (SMOKE_TEST) "smoke" else "full"),
    width = 7, height = 5)
yrange <- range(c(results$exact_mean, results$wood_mean))
plot(results$n_cl, results$exact_mean,
     log = "xy", type = "b", pch = 16, lwd = 2, col = "firebrick",
     xlab = "n (number of clusters)", ylab = "time per evaluation (s, log scale)",
     ylim = yrange,
     main = "Marginal-likelihood evaluation: PG-Exact vs PG-Woodbury")
lines(results$n_cl, results$wood_mean, type = "b", pch = 17, lwd = 2,
      col = "steelblue")
legend("topleft", legend = c("PG-Exact (dense Matérn)", "PG-Woodbury (SPDE)"),
       col = c("firebrick", "steelblue"), pch = c(16, 17), lwd = 2, bty = "n")
grid()
dev.off()

cat(sprintf("\nResults saved to dev/09_scalability_%s.rds and .pdf\n",
            if (SMOKE_TEST) "smoke" else "full"))
