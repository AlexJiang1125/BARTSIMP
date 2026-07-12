# =============================================================================
# Benchmark: replace `Matrix::determinant(Q, logarithm=TRUE)` calls with
# log-det readouts from a cached Cholesky factor.
#
# In `margprob_woodbury`, we already have Fchol = Cholesky(Qu). Computing
# log|Qu| via Matrix::determinant(Qu, log=TRUE) re-factorizes (LU). We should
# instead read it off Fchol's diagonal. Same idea for log|Q_psi|.
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
})

set.seed(2026)
n_cl <- 250
locs <- cbind(runif(n_cl), runif(n_cl))
mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(0.08, 0.25), cutoff = 0.03)
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = c(0.1, 0.5),
                                  prior.sigma = c(1.0, 0.5))
A     <- as(INLA::inla.spde.make.A(mesh = mesh, loc = locs), "CsparseMatrix")
Omega <- runif(n_cl, 0.5, 2.0)
yhat  <- rnorm(n_cl)
Q_psi <- INLA::inla.spde.precision(spde, theta = c(log(0.3), log(0.5)))

# ---- OLD: log-dets via Matrix::determinant ---------------------------------
margprob_old <- function(yhat, Omega, A, Q_psi) {
  AtWA  <- Matrix::crossprod(A, as.vector(Omega) * A)
  Qu    <- Matrix::forceSymmetric(AtWA + Q_psi)
  Fchol <- Matrix::Cholesky(Qu, LDL = FALSE, super = FALSE, perm = FALSE)
  rhs   <- as.numeric(Matrix::crossprod(A, Omega * yhat))
  m_u   <- as.numeric(Matrix::solve(Fchol, rhs, system = "A"))
  quad  <- sum(Omega * yhat^2) - sum(rhs * m_u)
  logdet_W    <- sum(log(Omega))
  logdet_Qpsi <- as.numeric(Matrix::determinant(Q_psi, logarithm = TRUE)$modulus)
  logdet_Qu   <- as.numeric(Matrix::determinant(Qu,    logarithm = TRUE)$modulus)
  N <- length(yhat)
  mlik <- 0.5 * (-N * log(2 * pi) + logdet_W + logdet_Qpsi - logdet_Qu - quad)
  list(mlik = mlik, m_u = m_u, Fchol = Fchol)
}

# ---- NEW: reuse the Cholesky factor we already have for log|Q_u| -----------
# CHMfactor's `determinant` reads the pre-computed diagonal directly --
# O(n), no materialization. Leave the Q_psi log-det alone (it's already fast).
margprob_new <- function(yhat, Omega, A, Q_psi) {
  AtWA  <- Matrix::crossprod(A, as.vector(Omega) * A)
  Qu    <- Matrix::forceSymmetric(AtWA + Q_psi)
  Fchol <- Matrix::Cholesky(Qu, LDL = FALSE, super = FALSE, perm = FALSE)
  rhs   <- as.numeric(Matrix::crossprod(A, Omega * yhat))
  m_u   <- as.numeric(Matrix::solve(Fchol, rhs, system = "A"))
  quad  <- sum(Omega * yhat^2) - sum(rhs * m_u)
  logdet_W    <- sum(log(Omega))
  logdet_Qpsi <- as.numeric(Matrix::determinant(Q_psi, logarithm = TRUE)$modulus)
  # log|Q_u|: read off the Cholesky we already have
  logdet_Qu   <- 2 * as.numeric(Matrix::determinant(Fchol, logarithm = TRUE)$modulus)
  N <- length(yhat)
  mlik <- 0.5 * (-N * log(2 * pi) + logdet_W + logdet_Qpsi - logdet_Qu - quad)
  list(mlik = mlik, m_u = m_u, Fchol = Fchol)
}

# ---- Equivalence ------------------------------------------------------------
r_old <- margprob_old(yhat, Omega, A, Q_psi)
r_new <- margprob_new(yhat, Omega, A, Q_psi)
cat(sprintf("mlik old = %.12f\nmlik new = %.12f\ndiff     = %.3e\n",
            r_old$mlik, r_new$mlik, abs(r_old$mlik - r_new$mlik)))
stopifnot(abs(r_old$mlik - r_new$mlik) < 1e-9)
stopifnot(max(abs(r_old$m_u - r_new$m_u)) < 1e-9)
cat("mlik and m_u identical to ~1e-9.\n")

# ---- Benchmark --------------------------------------------------------------
N_B <- 1000L
t_old <- system.time(for (i in seq_len(N_B)) r <- margprob_old(yhat, Omega, A, Q_psi))
t_new <- system.time(for (i in seq_len(N_B)) r <- margprob_new(yhat, Omega, A, Q_psi))

cat(sprintf("\n=== %d calls (mesh n=%d) ===\n", N_B, mesh$n))
cat(sprintf("  margprob_old (determinant) : %.3f s  (%.1f us/call)\n",
            t_old["elapsed"], 1e6 * t_old["elapsed"] / N_B))
cat(sprintf("  margprob_new (Cholesky)    : %.3f s  (%.1f us/call)\n",
            t_new["elapsed"], 1e6 * t_new["elapsed"] / N_B))
cat(sprintf("  speedup                    : %.2fx\n",
            t_old["elapsed"] / t_new["elapsed"]))

# Projection for one chain: 2 margprob calls / iter * 1500 iter
n_iter <- 1500L
cat(sprintf("\nProjected savings, 1 chain x %d iter (2 calls/iter):\n", n_iter))
cat(sprintf("  old : %.2f s\n", 2 * n_iter * t_old["elapsed"] / N_B))
cat(sprintf("  new : %.2f s\n", 2 * n_iter * t_new["elapsed"] / N_B))
cat(sprintf("  saved %.2f s (%.1f%%)\n",
            2 * n_iter * (t_old["elapsed"] - t_new["elapsed"]) / N_B,
            100 * (1 - t_new["elapsed"] / t_old["elapsed"])))
