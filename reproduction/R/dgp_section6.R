# =============================================================================
# Section 6 application: SYNTHETIC Kenya-DHS-like data-generating process
# =============================================================================
# The empirical analysis in Section 6 ("Application") of
#
#   Jiang, A. Z., & Wakefield, J. (2025). BARTSIMP: Flexible spatial covariate
#   modeling and prediction using Bayesian Additive Regression Trees.
#   Spatial and Spatio-temporal Epidemiology, 55, 100757.
#
# uses child wasting (weight-for-height Z-score, WHZ) from the 2014 Kenya DHS,
# with six geospatial covariates and population-weighted aggregation to
# administrative areas. The DHS microdata are NOT redistributable: you must
# register with the DHS Program and download them yourself (see
# reproduction/R/covariates_kenya.R for data access, and the Section 6
# vignette for the full recipe).
#
# To make the *workflow* reproducible without the restricted data, this file
# generates a SYNTHETIC dataset with the same STRUCTURE as the real one:
#
#   * clusters (enumeration areas) at point locations over a Kenya-like box;
#   * several household observations per cluster, sharing the cluster's
#     location and covariate values (the paper's within-unit assumption);
#   * a stratified design: strata = admin-1 area x urban/rural (the paper has
#     92 strata = 47 counties x urban/rural);
#   * six named covariates with realistic spatial smoothness;
#   * an under-five population-density weight per cell (for areal aggregation);
#   * nested administrative labels (admin1, admin2) for small-area estimation.
#
# It reuses the Matern GRF sampler and helpers from dgp_section5.R. Everything
# here is synthetic --- it reproduces the pipeline, not the paper's numbers.
# =============================================================================

if (!exists("simulate_grf", mode = "function")) {
  .this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  source(file.path(.this_dir, "dgp_section5.R"))
}

# The six geospatial covariates used in the paper (Table of Section 2).
SECTION6_COVARIATES <- c("pop_density", "night_light", "ndvi",
                         "temp", "precip", "access_city")

# A smooth, spatially structured covariate surface on [0,1]^2, built from a few
# random Fourier modes so neighbouring cells share values like a real raster.
# Returned on the raw z-scale; the generator rescales/kinks per covariate.
.smooth_surface <- function(s1, s2, n_modes = 4L) {
  val <- numeric(length(s1))
  for (k in seq_len(n_modes)) {
    a1 <- runif(1, 0.5, 3); a2 <- runif(1, 0.5, 3)
    ph <- runif(1, 0, 2 * pi); wt <- rnorm(1)
    val <- val + wt * sin(2 * pi * (a1 * s1 + a2 * s2) + ph)
  }
  as.numeric(scale(val))
}

# Build the six covariates over the grid. Each is a smooth spatial surface,
# transformed to a plausible marginal shape/scale (positive & right-skewed for
# density/lights/access; bounded for ndvi; roughly Gaussian for climate). These
# are stand-ins; with the real data you would extract raster values at the
# grid/cluster coordinates instead (see covariates_kenya.R).
.make_covariates <- function(s1, s2) {
  z <- lapply(SECTION6_COVARIATES, function(nm) .smooth_surface(s1, s2))
  names(z) <- SECTION6_COVARIATES
  data.frame(
    pop_density = exp(1.2 * z$pop_density),                 # >0, right-skewed
    night_light = pmax(0, exp(0.9 * z$night_light) - 0.6),  # >=0, many small
    ndvi        = plogis(1.1 * z$ndvi),                     # in (0,1)
    temp        = 24 + 4 * z$temp,                          # deg C-ish
    precip      = pmax(0, 70 + 25 * z$precip),              # mm-ish, >=0
    access_city = exp(1.0 * z$access_city)                  # >0 travel time
  )
}

# Standardize covariates to comparable scales for the covariate effect g(x).
# (BART is scale-invariant, but this keeps the true signal interpretable.)
.z <- function(x) as.numeric(scale(x))

# The true covariate effect g(x) on WHZ. Sign/shape chosen to echo the paper's
# partial-dependence findings: WHZ increases with population density and
# precipitation, decreases with temperature and access-to-city travel time.
.covariate_effect <- function(cov) {
  0.45 * .z(log(cov$pop_density)) +
    0.25 * .z(cov$precip) -
    0.35 * .z(cov$temp) -
    0.30 * .z(log1p(cov$access_city)) +
    0.10 * .z(cov$ndvi)
}

# Generate one synthetic Kenya-DHS-like dataset.
#
#   n_side       : prediction grid is n_side x n_side over the Kenya-like box.
#   n_clusters   : number of sampled enumeration areas (paper: 1584).
#   obs_range    : households sampled per cluster (paper: ~25); WHZ per child.
#   n_admin1     : number of admin-1 areas (paper Kenya: 47 counties).
#   admin2_per   : admin-2 areas per admin-1 (paper Kenya: ~290 total / 47 ~ 6).
#   omega        : covariate-signal share, as in Section 5. omega = 0.7 gives a
#                  covariate-driven-but-spatially-correlated surface like WHZ.
#   sigma2_e     : observation (child-level) noise variance.
#   kappa,nu,sigma2_m : Matern spatial-field hyperparameters.
#   urban_frac   : fraction of clusters that are urban (oversampled in DHS).
#   seed         : RNG seed.
#
# Returns a list:
#   $obs    : one row per child (cluster, admin1, admin2, stratum, urban, s1, s2,
#             the six covariates, and WHZ = y).
#   $grid   : the prediction surface, one row per cell (cell, s1, s2, admin1,
#             admin2, pop_u5 = under-5 population weight, the six covariates,
#             z_star, g_cov, f_true = target WHZ surface).
#   $clusters : one row per cluster (design/location table).
#   $params : provenance.
simulate_section6 <- function(n_side = 40L,
                              n_clusters = 400L,
                              obs_range = 15:25,
                              n_admin1 = 12L,
                              admin2_per = 6L,
                              omega = 0.7,
                              sigma2_e = 1,
                              kappa = 3.5, nu = 1, sigma2_m = 0.6,
                              urban_frac = 0.35,
                              whz_shift = -0.6,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ---- prediction grid over a Kenya-like unit box ----------------------------
  mid  <- (seq_len(n_side) - 0.5) / n_side
  half <- 0.5 / n_side
  grid <- expand.grid(s1 = mid, s2 = mid)
  n_cells <- nrow(grid)
  grid$cell <- seq_len(n_cells)

  # ---- nested administrative partition (rectangular tiling of the box) -------
  # admin1: a near-square grid of n_admin1 tiles; admin2: split each into
  # admin2_per sub-tiles. Labels are assigned by which tile a cell falls in.
  a1_cols <- ceiling(sqrt(n_admin1)); a1_rows <- ceiling(n_admin1 / a1_cols)
  a1_col  <- pmin(a1_cols, floor(grid$s1 * a1_cols) + 1L)
  a1_row  <- pmin(a1_rows, floor(grid$s2 * a1_rows) + 1L)
  grid$admin1 <- pmin(n_admin1, (a1_row - 1L) * a1_cols + a1_col)

  a2_cols <- ceiling(sqrt(admin2_per)); a2_rows <- ceiling(admin2_per / a2_cols)
  # position within the admin1 tile, in [0,1)
  fx <- (grid$s1 * a1_cols) %% 1; fy <- (grid$s2 * a1_rows) %% 1
  a2_col <- pmin(a2_cols, floor(fx * a2_cols) + 1L)
  a2_row <- pmin(a2_rows, floor(fy * a2_rows) + 1L)
  a2_within <- pmin(admin2_per, (a2_row - 1L) * a2_cols + a2_col)
  grid$admin2 <- (grid$admin1 - 1L) * admin2_per + a2_within

  # ---- covariates, under-5 population weight, spatial field, truth -----------
  cov_grid <- .make_covariates(grid$s1, grid$s2)
  grid <- cbind(grid, cov_grid)
  # under-5 population density weight: tied to (but noisier than) pop_density
  grid$pop_u5 <- cov_grid$pop_density * exp(rnorm(n_cells, 0, 0.3))

  grid$z_star <- simulate_grf(grid[, c("s1", "s2")],
                              kappa = kappa, nu = nu, sigma2_m = sigma2_m)
  g_cov <- .covariate_effect(cov_grid)
  grid$g_cov  <- g_cov
  # blend covariate signal and spatial field, then shift onto the WHZ scale
  grid$f_true <- whz_shift + (1 - omega) * grid$z_star + omega * g_cov

  # ---- sample clusters (enumeration areas) -----------------------------------
  sel  <- sample.int(n_cells, size = n_clusters, replace = TRUE)
  urban <- rbinom(n_clusters, 1, urban_frac)
  n_j  <- sample(obs_range, n_clusters, replace = TRUE)

  clusters <- data.frame(
    cluster = seq_len(n_clusters),
    cell    = sel,
    admin1  = grid$admin1[sel],
    admin2  = grid$admin2[sel],
    urban   = urban,
    # stratum = admin1 x urban/rural (mirrors the paper's 92 strata)
    stratum = paste0("a1_", grid$admin1[sel], ifelse(urban == 1, "_U", "_R")),
    s1      = grid$s1[sel] + runif(n_clusters, -half, half),
    s2      = grid$s2[sel] + runif(n_clusters, -half, half)
  )
  clusters <- cbind(clusters, cov_grid[sel, , drop = FALSE])

  # ---- expand to per-child observations --------------------------------------
  idx <- rep(seq_len(n_clusters), times = n_j)
  obs <- data.frame(
    cluster = clusters$cluster[idx],
    admin1  = clusters$admin1[idx],
    admin2  = clusters$admin2[idx],
    stratum = clusters$stratum[idx],
    urban   = clusters$urban[idx],
    s1      = clusters$s1[idx],
    s2      = clusters$s2[idx]
  )
  obs <- cbind(obs, cov_grid[clusters$cell[idx], , drop = FALSE])
  # per-child WHZ: true surface at the cluster's cell + child-level noise
  obs$y <- grid$f_true[clusters$cell[idx]] + rnorm(nrow(obs), 0, sqrt(sigma2_e))

  list(
    obs = obs,
    grid = grid[, c("cell", "s1", "s2", "admin1", "admin2", "pop_u5",
                    SECTION6_COVARIATES, "z_star", "g_cov", "f_true")],
    clusters = clusters,
    params = list(n_side = n_side, n_clusters = n_clusters,
                  n_admin1 = n_admin1, admin2_per = admin2_per, omega = omega,
                  sigma2_e = sigma2_e, kappa = kappa, nu = nu,
                  sigma2_m = sigma2_m, urban_frac = urban_frac,
                  whz_shift = whz_shift, seed = seed)
  )
}
