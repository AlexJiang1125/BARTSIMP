# Quick smoke test for cv_section6.R on synthetic Section 6 data.
suppressMessages(library(INLA))   # bartsimp() sampler needs INLA on the search path
.repro_R <- normalizePath("R")
source(file.path(.repro_R, "run_scenario.R"))   # metrics/competitors/bartsimp_fit
source(file.path(.repro_R, "dgp_section5.R"))
source(file.path(.repro_R, "dgp_section6.R"))
source(file.path(.repro_R, "areal.R"))
source(file.path(.repro_R, "cv_section6.R"))

sim <- simulate_section6(n_side = 24L, n_clusters = 200L, obs_range = 15:25,
                         n_admin1 = 9L, admin2_per = 4L, seed = 20260712L)
cat(sprintf("sim: %d children, %d clusters, %d grid cells, %d strata\n",
            nrow(sim$obs), nrow(sim$clusters), nrow(sim$grid),
            length(unique(sim$clusters$stratum))))

t0 <- Sys.time()
res <- run_cv_fold_section6(
  sim, methods = c("BARTSIMP", "BART", "SPDE", "SPDE0"),
  prop_train = 0.8, alpha = 0.05, seed = 1L,
  bartsimp_args = list(ntree = 30L, ndpost = 120L, nskip = 1L, nwarmup = 120L),
  verbose = TRUE)
cat(sprintf("\nfold wall time: %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(res, row.names = FALSE)
