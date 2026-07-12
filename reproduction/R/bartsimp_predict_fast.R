# =============================================================================
# Fast (Woodbury) spatial predictor for BARTSIMP  -- opt-in, drop-in
# =============================================================================
# bartsimp_predict() (in the package) refits an INLA SPDE model *per posterior
# draw* to recover the spatial residual field. With the Matern hyperparameters
# held FIXED at each draw's (kappa, sigma_m), that INLA fit is just the posterior
# mean of a Gaussian Markov random field, which has a closed form. This routine
# computes it directly via a sparse Cholesky (the "Woodbury" form), giving the
# same per-draw spatial mean at a fraction of the cost and with no INLA call in
# the loop.
#
# For draw m, with residual r = y_train - tree_draws_train[m, ], per-observation
# precision Omega = 1/sigma_obs[m]^2, projection A (obs -> mesh nodes), and SPDE
# precision Q_psi(kappa_m, sigma_m):
#
#   Q_u = A' Omega A + Q_psi           (sparse, n_mesh x n_mesh)
#   m_u = Q_u^{-1} A' Omega r           (GMRF posterior mean of the field weights)
#   spatial_mean_test = A_test %*% m_u
#
# This matches INLA::inla(...)$summary.fitted.values[pred, "mean"] for the
# intercept-free field model  y ~ -1 + f(i, model = spde)  with fixed
# hyperparameters. Validated against the INLA path in
# reproduction/R/_validate_fast_predict.R.
#
# Contract and argument names mirror BARTSIMP::bartsimp_predict() so this is a
# drop-in replacement, EXCEPT sigma_obs is required here (the closed form needs
# the observation precision; the INLA path could estimate it instead).
# =============================================================================

# Sparse SPDE precision for fixed (kappa, sigma_m). rho = sqrt(8 nu)/kappa is the
# Matern range; theta convention is (log(range), log(sigma_m)).
.qpsi_from_kappa <- function(spde, kappa, sigma_m, nu = 1) {
  rho <- sqrt(8 * nu) / kappa
  INLA::inla.spde.precision(spde, theta = c(log(rho), log(sigma_m)))
}

# GMRF posterior mean of the mesh weights, m_u = Q_u^{-1} A' Omega r.
.woodbury_mu <- function(r, Omega, A, Q_psi) {
  AtWA <- Matrix::crossprod(A, Omega * A)
  Qu   <- Matrix::forceSymmetric(AtWA + Q_psi)
  Fch  <- Matrix::Cholesky(Qu, LDL = FALSE, super = FALSE, perm = TRUE)
  rhs  <- as.numeric(Matrix::crossprod(A, Omega * r))
  as.numeric(Matrix::solve(Fch, rhs, system = "A"))
}

bartsimp_predict_woodbury <- function(y_train,
                                      tree_draws_train,
                                      tree_draws_test,
                                      s_train,
                                      s_test,
                                      kappas,
                                      sigmams,
                                      sigma_obs,
                                      mesh,
                                      nu = 1,
                                      verbose = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("Package 'INLA' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.")
  if (missing(sigma_obs) || is.null(sigma_obs))
    stop("bartsimp_predict_woodbury() requires sigma_obs (per-draw noise SD).")

  tree_draws_train <- as.matrix(tree_draws_train)
  tree_draws_test  <- as.matrix(tree_draws_test)
  s_train <- as.data.frame(s_train); colnames(s_train) <- c("s1", "s2")
  s_test  <- as.data.frame(s_test);  colnames(s_test)  <- c("s1", "s2")

  # Match bartsimp_predict()'s internal sort: order training by (s1, s2) and
  # apply the same permutation to the tree-draw columns.
  ord_train <- order(s_train$s1, s_train$s2)
  y_train          <- y_train[ord_train]
  s_train          <- s_train[ord_train, , drop = FALSE]
  tree_draws_train <- tree_draws_train[, ord_train, drop = FALSE]

  ndpost <- nrow(tree_draws_train)
  n_test <- nrow(s_test)
  stopifnot(length(y_train) == ncol(tree_draws_train),
            nrow(tree_draws_test) == ndpost,
            ncol(tree_draws_test) == n_test,
            length(kappas) == ndpost, length(sigmams) == ndpost,
            length(sigma_obs) == ndpost)

  # SPDE object: priors are irrelevant because we set theta (the hyperparameters)
  # explicitly per draw; only the mesh-dependent finite-element structure is used.
  spde <- INLA::inla.spde2.pcmatern(mesh = mesh,
                                    prior.range = c(0.5, 0.5),
                                    prior.sigma = c(0.5, 0.5))
  A_train <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(s_train))
  A_test  <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(s_test))

  spatial_mean_test <- matrix(NA_real_, nrow = ndpost, ncol = n_test)
  for (m in seq_len(ndpost)) {
    if (verbose && (m %% 50 == 0 || m == 1))
      message("Woodbury draw ", m, " / ", ndpost)
    r_m     <- y_train - tree_draws_train[m, ]
    Omega_m <- rep(1 / sigma_obs[m]^2, length(r_m))
    Q_psi   <- .qpsi_from_kappa(spde, kappas[m], sigmams[m], nu = nu)
    m_u     <- .woodbury_mu(r_m, Omega_m, A_train, Q_psi)
    spatial_mean_test[m, ] <- as.numeric(A_test %*% m_u)
  }

  list(
    spatial_mean_test    = spatial_mean_test,
    prediction_mean_test = spatial_mean_test + tree_draws_test,
    s_train = s_train,
    s_test  = s_test
  )
}
