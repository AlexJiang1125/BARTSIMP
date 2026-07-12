#!/usr/bin/env Rscript
# smoke_paired.R — single-rep smoke test of the paired A/B:
#   BARTSIMP-PG (d0=0, legacy)  vs  BARTSIMP-PG-SD (d0=2, depth-floor)
# on ONE dgp5 dataset, matched truth+data. Confirms wiring + rough timing +
# that the two arms produce DIFFERENT admin posteriors (knob has an effect).
suppressPackageStartupMessages({ library(optparse) })
CLUSTER <- "/Users/alexziyujiang/Documents/Claude/Projects/BARTSIMP_glm/bartsimp_cluster"
setwd(CLUSTER)
b07 <- "experiments/sim_thm_verification/sim_07/R"
source(file.path(b07, "dgp.R"))
source(file.path(b07, "methods_inla.R"))
source(file.path(b07, "methods_bartsimp.R"))
source(file.path(b07, "coverage.R"))
source("experiments/sim_dgp5_coverage/R/dgp5_core.R")

N <- 500L; A <- 20L; n_iter <- 600L; burn <- 200L; m_trees <- 20L; S <- 500L
d0 <- as.integer(ceiling(sqrt(log(N))))   # spec floor; log(500)=6.21 -> sqrt=2.49 -> 3
cat(sprintf("ceil(sqrt(log(%d))) = %d\n", N, d0))

truth <- make_truth5(alpha_0 = -0.90, seed_truth = 4242L)
admin_struct <- make_admin_structure(A = A, seed = 1L)
set.seed(90000L + 1L)
design   <- draw_design5(N)
admin_id <- assign_admin(design$locs, admin_struct)
tr       <- apply_truth5(design, truth)

score_one <- function(p_post, label) {
  ad <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.20)
  ad95 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.05)
  rmse <- sqrt(mean((ad$p_a_post_mean - ad$p_a_true)^2))
  bias <- mean(ad$p_a_post_mean - ad$p_a_true)
  psd  <- mean(ad$p_a_post_sd)
  cat(sprintf("%-18s cov80=%.3f cov95=%.3f rmse=%.4f bias=%+.4f post_sd=%.4f ratio=%.3f\n",
              label, mean(ad$covered), mean(ad95$covered), rmse, bias, psd, rmse/psd))
  invisible(list(cov80 = mean(ad$covered), bias = bias, post_sd = psd))
}

cat("\n-- BARTSIMP-PG (d0=0) --\n")
set.seed(424242L)
t0 <- Sys.time()
p_leg <- fit_bartsimp_pg(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                         S = S, n_iter = n_iter, burn = burn, m_trees = m_trees)
cat(sprintf("  elapsed %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
r_leg <- score_one(p_leg, "BARTSIMP-PG")

cat("\n-- BARTSIMP-PG-SD (d0=", d0, ") --\n", sep = "")
set.seed(424242L)
t0 <- Sys.time()
p_sd <- fit_bartsimp_pg_sd(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                           S = S, n_iter = n_iter, burn = burn, m_trees = m_trees, d0 = d0)
cat(sprintf("  elapsed %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
r_sd <- score_one(p_sd, "BARTSIMP-PG-SD")

cat(sprintf("\nPAIRED DELTA (SD - legacy): dcov80=%+.3f  dbias=%+.4f  dpost_sd=%+.4f\n",
            r_sd$cov80 - r_leg$cov80, r_sd$bias - r_leg$bias, r_sd$post_sd - r_leg$post_sd))
cat(sprintf("arms differ: %s\n", !isTRUE(all.equal(p_leg, p_sd))))
