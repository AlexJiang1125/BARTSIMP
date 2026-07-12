# =============================================================================
# Section 5 simulation: data-generating process (DGP)
# =============================================================================
# Reproduces the simulation design of Section 5 ("Simulation Experiments") of
#
#   Jiang, A. Z., & Wakefield, J. (2025). BARTSIMP: Flexible spatial covariate
#   modeling and prediction using Bayesian Additive Regression Trees.
#   Spatial and Spatio-temporal Epidemiology, 55, 100757.
#
# Design (paper Section 5.1):
#   * 50 x 50 grid over [0,1]^2; |G| = 2500 cells.
#   * Two covariates x1, x2 ~ Uniform(0,1), drawn independently per cell.
#   * Baseline spatial field z* ~ GRF, Matern(kappa = 2.5, nu = 1, sigma_m^2 = 0.5).
#     kappa = 2.5 corresponds to range rho = sqrt(8 nu)/kappa ~= 1.13.
#   * Covariate surface f0(x): "tree", "linear", or "smooth" (see below).
#   * True field:  f(x_g) = (1 - omega) z*_g + omega f0(x_g).
#     omega in {1, 0.8, 0.5, 0.2, 0} = scenarios 1..5 (covariate-signal share).
#   * 250 clusters: sample 250 grid cells, one point uniformly within each cell.
#   * Per cluster, n_j ~ Uniform{5,...,10} observations.
#   * y_ij = f(cell_j) + N(0, sigma_e^2 = 1).
#
# These functions have no dependency on the BARTSIMP package; they only need
# base R (the GRF uses a dense Matern Cholesky, which is fine for 2500 cells).
# =============================================================================

# Matern correlation with the paper's kappa parameterization (nu fixed at 1).
# C(0) = 1; C(d) = 2^(1-nu)/Gamma(nu) * (kappa d)^nu * K_nu(kappa d).
matern_corr <- function(d, kappa = 2.5, nu = 1) {
  out <- (2^(1 - nu) / gamma(nu)) * (kappa * d)^nu * besselK(kappa * d, nu)
  out[d == 0] <- 1
  out
}

# Draw one baseline spatial field z* on the given coordinates.
#   coords: n x 2 matrix/data.frame of locations.
#   Returns a length-n numeric vector with marginal variance sigma2_m.
simulate_grf <- function(coords, kappa = 2.5, nu = 1, sigma2_m = 0.5,
                         jitter = 1e-8) {
  coords <- as.matrix(coords)
  D <- as.matrix(dist(coords))
  Sigma <- sigma2_m * matern_corr(D, kappa = kappa, nu = nu)
  # small nugget on the diagonal for numerical positive-definiteness
  diag(Sigma) <- diag(Sigma) + jitter
  L <- chol(Sigma)                       # upper-triangular, Sigma = L' L
  as.numeric(crossprod(L, rnorm(nrow(coords))))
}

# Deterministic covariate surface f0(x1, x2), paper Section 5.1.
covariate_surface <- function(x1, x2, type = c("tree", "linear", "smooth")) {
  type <- match.arg(type)
  switch(type,
    tree = ifelse(x1 < 0.5, 3,
                  ifelse(x2 >= 0.5, 0, -2)),
    linear = 2 * x1 - x2,
    smooth = sin(2 * pi * x1) + cos(2 * pi * x2)
  )
}

# Generate one simulated dataset for a given scenario.
#
#   omega    : covariate-signal share in {1, 0.8, 0.5, 0.2, 0}.
#   surface  : "tree" | "linear" | "smooth".
#   n_side   : grid is n_side x n_side (default 50 -> |G| = 2500).
#   n_clusters, obs_range, sigma2_e, kappa, nu, sigma2_m : as in the paper.
#   seed     : RNG seed; covariate + spatial surfaces are redrawn each call
#              (the paper does NOT hold them fixed across replicates).
#
# Returns a list:
#   $obs   : per-observation training data (one row per observation) with
#            columns cluster, cell, s1, s2, x1, x2, y.
#   $grid  : the full 2500-cell surface with columns cell, s1, s2, x1, x2,
#            z_star, f0, f_true  (f_true is the prediction target).
#   $params: the scenario settings, for provenance.
simulate_section5 <- function(omega,
                              surface = c("tree", "linear", "smooth"),
                              n_side = 50L,
                              n_clusters = 250L,
                              obs_range = 5:10,
                              sigma2_e = 1,
                              kappa = 2.5, nu = 1, sigma2_m = 0.5,
                              seed = NULL) {
  surface <- match.arg(surface)
  if (!is.null(seed)) set.seed(seed)

  # ---- grid cell centers over [0,1]^2 (cell midpoints) ----
  mid <- (seq_len(n_side) - 0.5) / n_side
  half <- 0.5 / n_side
  grid <- expand.grid(s1 = mid, s2 = mid)
  n_cells <- nrow(grid)                       # 2500 for n_side = 50
  grid$cell <- seq_len(n_cells)

  # ---- per-cell covariates ~ U(0,1) ----
  grid$x1 <- runif(n_cells)
  grid$x2 <- runif(n_cells)

  # ---- baseline spatial field and covariate surface ----
  grid$z_star <- simulate_grf(grid[, c("s1", "s2")],
                              kappa = kappa, nu = nu, sigma2_m = sigma2_m)
  grid$f0 <- covariate_surface(grid$x1, grid$x2, type = surface)
  grid$f_true <- (1 - omega) * grid$z_star + omega * grid$f0

  # ---- sample 250 clusters: one point uniformly inside each chosen cell ----
  sel <- sample.int(n_cells, size = n_clusters, replace = FALSE)
  n_j <- sample(obs_range, n_clusters, replace = TRUE)

  cell_of <- rep(sel, times = n_j)
  # jitter the cluster location uniformly within its cell
  s1_pt <- rep(grid$s1[sel] + runif(n_clusters, -half, half), times = n_j)
  s2_pt <- rep(grid$s2[sel] + runif(n_clusters, -half, half), times = n_j)

  obs <- data.frame(
    cluster = rep(seq_len(n_clusters), times = n_j),
    cell    = cell_of,
    s1      = s1_pt,
    s2      = s2_pt,
    x1      = grid$x1[cell_of],
    x2      = grid$x2[cell_of]
  )
  obs$y <- grid$f_true[cell_of] + rnorm(nrow(obs), 0, sqrt(sigma2_e))

  list(
    obs = obs,
    grid = grid[, c("cell", "s1", "s2", "x1", "x2", "z_star", "f0", "f_true")],
    params = list(omega = omega, surface = surface, n_side = n_side,
                  n_clusters = n_clusters, sigma2_e = sigma2_e,
                  kappa = kappa, nu = nu, sigma2_m = sigma2_m, seed = seed)
  )
}
