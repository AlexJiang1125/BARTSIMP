# =============================================================================
# Local 4-way head-to-head: BARTSIMP-PG vs BART-logistic, SPDE-binomial,
# SPDE0-binomial. Same toy dataset (paper-aligned generator, balanced
# scenario at omega=0.5) used in dev/13, dev/18, dev/19, dev/20.
#
# Each method is evaluated on the same X_grid / locs_grid against the same
# eta_true / p_true ground truth using the V3 metrics in
# dev/competing/metrics.R.
# =============================================================================

# Load BARTSIMP-PG sampler (and rbart.R via its auto-loader)
source("dev/bartsimp_pg.R")
source("dev/rbart.R")

# Load competing-method wrappers
source("dev/competing/bart_logistic.R")
source("dev/competing/spde_binomial.R")
source("dev/competing/spde0_binomial.R")
source("dev/competing/bartsimp_pg_wrap.R")
source("dev/competing/metrics.R")

set.seed(2026)

# ---- toy config (same as dev/18, dev/19, dev/20) --------------------------
N_CL     <- 150L
N_ITER   <- 1500L
BURN     <-  500L
GRID_DIM <-   20L
OMEGA    <- 0.5

# Paper-aligned generator
sigma_m_baseline <- 0.5
rho_baseline     <- sqrt(8 * 1) / 2.5
surface_f0 <- function(x1, x2)
  ifelse(x1 < 0.5, 3, ifelse(x2 < 0.5, -2, 0))

G    <- GRID_DIM^2
step <- 1 / GRID_DIM
gx   <- (seq_len(GRID_DIM) - 0.5) * step
cell_centers <- expand.grid(lon = gx, lat = gx)
x1_grid <- runif(G); x2_grid <- runif(G)

Dg    <- as.matrix(dist(cell_centers))
kappa <- sqrt(8 * 1) / rho_baseline
Sg    <- sigma_m_baseline^2 * matern_corr(Dg, kappa, 1) + 1e-6 * diag(G)
z_grid   <- as.vector(t(chol(Sg)) %*% rnorm(G))
f0_grid  <- surface_f0(x1_grid, x2_grid)
eta_grid <- (1 - OMEGA) * z_grid + OMEGA * f0_grid
p_grid   <- plogis(eta_grid)

cluster_cells <- sample.int(G, N_CL, replace = FALSE)
jitter_x <- runif(N_CL, -step/2, step/2)
jitter_y <- runif(N_CL, -step/2, step/2)
locs <- cbind(
  lon = cell_centers$lon[cluster_cells] + jitter_x,
  lat = cell_centers$lat[cluster_cells] + jitter_y
)
X   <- cbind(x1 = x1_grid[cluster_cells], x2 = x2_grid[cluster_cells])
eta_train <- eta_grid[cluster_cells]
p_train   <- plogis(eta_train)
n_i <- as.numeric(sample(5:15, N_CL, replace = TRUE))
k_i <- vapply(seq_len(N_CL),
              function(i) sum(rbinom(n_i[i], 1, p_train[i])),
              numeric(1))

X_grid    <- cbind(x1 = x1_grid, x2 = x2_grid)
locs_grid <- as.matrix(cell_centers[, c("lon", "lat")])

cat(sprintf("Dataset: n_cl=%d, grid=%d cells, omega=%.2f\n",
            N_CL, G, OMEGA))

# ---- run all four methods -------------------------------------------------
results <- list()

cat("\n=== BARTSIMP-PG (Path B + adaptive MH + m=100, k=3) ===\n")
results$bartsimp_pg <- fit_bartsimp_pg(
  X_train = X, locs_train = locs, n_i = n_i, k_i = k_i,
  X_pred  = X_grid, locs_pred = locs_grid,
  n_iter  = N_ITER, burn = BURN,
  m_trees = 100L, k_bart = 3,
  bd_mode = "conditional", adapt_psi = TRUE,
  verbose = TRUE
)

cat("\n=== BART-logistic (no spatial) ===\n")
results$bart_logistic <- fit_bart_logistic(
  X_train = X, locs_train = locs, n_i = n_i, k_i = k_i,
  X_pred  = X_grid, locs_pred = locs_grid,
  n_iter  = N_ITER, burn = BURN,
  ntree   = 50L, verbose = TRUE
)

cat("\n=== SPDE-binomial (linear + Matern via INLA) ===\n")
results$spde_binomial <- fit_spde_binomial(
  X_train = X, locs_train = locs, n_i = n_i, k_i = k_i,
  X_pred  = X_grid, locs_pred = locs_grid,
  verbose = FALSE
)

cat("\n=== SPDE0-binomial (intercept + Matern via INLA) ===\n")
results$spde0_binomial <- fit_spde0_binomial(
  X_train = X, locs_train = locs, n_i = n_i, k_i = k_i,
  X_pred  = X_grid, locs_pred = locs_grid,
  verbose = FALSE
)

# ---- V3 metrics on the grid -----------------------------------------------
# (no held-out k_obs / n_obs at the grid -- we have only the latent truth
# eta_grid / p_grid since the grid is a deterministic eval surface, not a
# held-out validation sample. Brier / log-loss are reported NA accordingly.)
cat("\n=== V3 metrics on the grid (truth = eta_grid, p_grid) ===\n")
mlist <- lapply(results, v3_metrics,
                eta_true = eta_grid, p_true = p_grid)
for (m in mlist) v3_print_row(m)

# ---- summary table ---------------------------------------------------------
tab <- do.call(rbind, lapply(mlist, function(m) {
  data.frame(
    method   = m$method,
    elapsed  = sprintf("%.1f", m$elapsed_sec),
    rmse_eta = sprintf("%.3f", m$rmse_eta),
    rmse_p   = sprintf("%.3f", m$rmse_p),
    acr_eta  = sprintf("%.3f", m$acr_eta),
    ail_eta  = sprintf("%.3f", m$ail_eta),
    ais_eta  = sprintf("%.2f", m$ais_eta),
    acr_p    = sprintf("%.3f", m$acr_p),
    ail_p    = sprintf("%.3f", m$ail_p),
    stringsAsFactors = FALSE
  )
}))
cat("\n=== Summary table ===\n")
print(tab, row.names = FALSE)

saveRDS(list(
  config = list(N_CL = N_CL, OMEGA = OMEGA, GRID_DIM = GRID_DIM,
                N_ITER = N_ITER, BURN = BURN),
  truth = list(eta_grid = eta_grid, p_grid = p_grid),
  results = results,
  metrics = mlist
), "dev/21_competing_test.rds")
cat("\nSaved dev/21_competing_test.rds\n")
