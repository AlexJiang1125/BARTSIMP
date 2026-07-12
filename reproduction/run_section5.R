#!/usr/bin/env Rscript
# =============================================================================
# Section 5 simulation driver
# =============================================================================
# Runs the Section 5 experiment (BARTSIMP vs BART / SPDE / SPDE0) over a grid of
# scenarios and replicates, then writes the replicate-level and scenario-mean
# tables to reproduction/results/.
#
# Usage (from the package root, with the framework R that has INLA installed):
#
#   Rscript reproduction/run_section5.R                # quick demo config
#   Rscript reproduction/run_section5.R --paper        # full paper-scale sweep
#
# NOTE ON RUNTIME
#   BARTSIMP dominates the cost (per-draw INLA in bartsimp_predict + INLA inside
#   the C++ sampler). The quick config below is meant to finish in minutes; the
#   --paper config (n_side=50, 250 clusters, ndpost=1000, 20 reps, all surfaces)
#   is a long multi-hour run. See reproduction/vignette_section5.Rmd for the
#   opt-in Woodbury fast sampler that removes the per-draw INLA bottleneck.
# =============================================================================

suppressWarnings(suppressMessages({
  if (!requireNamespace("BARTSIMP", quietly = TRUE))
    stop("Install the BARTSIMP package first (see reproduction/README.md).")
  if (requireNamespace("INLA", quietly = TRUE))
    try(INLA::inla.setOption(fmesher.evolution.warn = FALSE), silent = TRUE)
}))

# locate reproduction/R relative to this script
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
root <- if (length(.f)) dirname(normalizePath(.f)) else file.path(getwd(), "reproduction")
.repro_R <- file.path(root, "R")
source(file.path(.repro_R, "run_scenario.R"))

paper <- "--paper" %in% commandArgs(trailingOnly = TRUE)

if (paper) {
  cfg <- list(
    omegas = c(1, 0.8, 0.5, 0.2, 0), surfaces = c("tree", "linear", "smooth"),
    reps = 1:20, n_side = 50L, n_clusters = 250L,
    bartsimp_args = list(ntree = 50L, ndpost = 1000L, nwarmup = 1000L),
    bart_args = list(ntree = 50L, ndpost = 2000L, nskip = 2000L)
  )
  message("Running FULL paper-scale sweep (this takes many hours).")
} else {
  cfg <- list(
    omegas = c(1, 0.5, 0), surfaces = "tree",
    reps = 1:2, n_side = 25L, n_clusters = 150L,
    bartsimp_args = list(ntree = 40L, ndpost = 200L, nwarmup = 200L),
    bart_args = list(ntree = 40L, ndpost = 1000L, nskip = 1000L)
  )
  message("Running QUICK demo config (pass --paper for the full sweep).")
}

t0 <- Sys.time()
results <- run_section5(
  omegas = cfg$omegas, surfaces = cfg$surfaces, reps = cfg$reps,
  n_side = cfg$n_side, n_clusters = cfg$n_clusters,
  bartsimp_args = cfg$bartsimp_args, bart_args = cfg$bart_args, verbose = TRUE
)
message(sprintf("Sweep finished in %.1f min.",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

out_dir <- file.path(root, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (paper) "paper" else "demo"
write.csv(results, file.path(out_dir, sprintf("section5_%s_replicates.csv", tag)),
          row.names = FALSE)
summ <- summarize_section5(results)
write.csv(summ, file.path(out_dir, sprintf("section5_%s_summary.csv", tag)),
          row.names = FALSE)

cat("\n==== scenario-mean summary ====\n")
print(summ, digits = 3, row.names = FALSE)
cat(sprintf("\nWrote results to %s\n", out_dir))
