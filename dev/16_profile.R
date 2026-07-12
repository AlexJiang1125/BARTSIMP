# =============================================================================
# Profile the sampler to identify the actual hot spots.
# =============================================================================

source("dev/bartsimp_pg.R")

set.seed(2026)
dat <- gen_data(seed = 2026, n_cl = 200,
                covariate_signal = TRUE, spatial_signal = TRUE)

prof_file <- tempfile(fileext = ".out")
Rprof(prof_file, interval = 0.01, line.profiling = TRUE,
      memory.profiling = FALSE)
fit <- run_sampler_bartsimp_pg(
  X = dat$X, locs = dat$locs, n_i = dat$n_i, k_i = dat$k_i,
  f_true = dat$f_true, z_true = dat$z_true,
  n_iter = 600L, burn = 200L, m_trees = 20L,
  mesh_max_edge = c(0.08, 0.25),
  bd_mode = "conditional",
  verbose_every = 0, label = "prof"
)
Rprof(NULL)

cat(sprintf("\nElapsed: %.1fs\n", fit$elapsed_sec))
cat("\n=== Top by self-time ===\n")
print(head(summaryRprof(prof_file)$by.self, 20))
cat("\n=== Top by total-time ===\n")
print(head(summaryRprof(prof_file)$by.total, 20))
