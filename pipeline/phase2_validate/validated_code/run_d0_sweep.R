#!/usr/bin/env Rscript
# run_d0_sweep.R — SENSITIVITY + ABLATION stress test for the depth-floor.
# Sweeps d0 in {0,1,2,3,4} on matched-seed dgp5 designs (truth drawn once, data
# per rep). d0 = 0 IS the ablation (no floor = legacy BARTSIMP-PG); d0 in {1..4}
# are the floored variants. Reports admin cov80 / cov95 / bias / post_sd / ratio
# and the pixel-latent guardrail per d0, so we can read (a) the SAFE OPERATING
# RANGE of d0 and (b) whether the advantage is a broad plateau or a narrow sweet
# spot. Only the BARTSIMP family runs (no roster) to keep the sweep affordable.
suppressPackageStartupMessages({ library(optparse) })
opts <- parse_args(OptionParser(option_list = list(
  make_option("--reps",    default = 16L,  type = "integer"),
  make_option("--N",       default = 500L, type = "integer"),
  make_option("--A",       default = 20L,  type = "integer"),
  make_option("--n_draws", default = 500L, type = "integer"),
  make_option("--n_iter",  default = 600L, type = "integer"),
  make_option("--burn",    default = 200L, type = "integer"),
  make_option("--m_trees", default = 20L,  type = "integer"),
  make_option("--d0_grid", default = "0,1,2,3,4"),
  make_option("--base_seed",  default = 90000L, type = "integer"),
  make_option("--seed_truth", default = 4242L,  type = "integer"),
  make_option("--alpha0",     default = -0.90,  type = "double"),
  make_option("--tag",        default = "d0sweep"),
  make_option("--out_dir",
              default = "/Users/alexziyujiang/Documents/GitHub/BARTSIMP/pipeline/phase2_validate/validated_results")
)))
CLUSTER <- "/Users/alexziyujiang/Documents/Claude/Projects/BARTSIMP_glm/bartsimp_cluster"
setwd(CLUSTER)
b07 <- "experiments/sim_thm_verification/sim_07/R"
source(file.path(b07, "dgp.R")); source(file.path(b07, "methods_inla.R"))
source(file.path(b07, "methods_bartsimp.R")); source(file.path(b07, "coverage.R"))
source("experiments/sim_dgp5_coverage/R/spatial_ml_competitors.R")
.carf <- "experiments/sim_dgp5_coverage/R/vendor/CAR-Forest.R"; if (file.exists(.carf)) source(.carf)
source("experiments/sim_dgp5_coverage/R/design_culture_competitors.R")
source("experiments/sim_dgp5_coverage/R/dgp5_core.R")
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
d0_grid <- as.integer(strsplit(opts$d0_grid, ",")[[1]])

compute_pixel_latent_cov <- function(p_post, p_star, alpha = 0.20) {
  N <- ncol(p_post); cov_i <- numeric(N)
  for (i in seq_len(N)) {
    iv <- quantile(p_post[, i], probs = c(alpha/2, 1 - alpha/2), names = FALSE, na.rm = TRUE)
    cov_i[i] <- as.integer(p_star[i] >= iv[1] & p_star[i] <= iv[2])
  }
  mean(cov_i)
}
truth        <- make_truth5(alpha_0 = opts$alpha0, seed_truth = opts$seed_truth)
admin_struct <- make_admin_structure(A = opts$A, seed = 1L)

per_rep <- list()
for (r in seq_len(opts$reps)) {
  set.seed(opts$base_seed + r)
  design   <- draw_design5(opts$N)
  admin_id <- assign_admin(design$locs, admin_struct)
  tr       <- apply_truth5(design, truth)
  for (d0 in d0_grid) {
    set.seed(opts$base_seed * 7L + r)        # matched sampler seed across d0
    p_post <- tryCatch({
      if (d0 == 0L)
        fit_bartsimp_pg(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                        S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn, m_trees = opts$m_trees)
      else
        fit_bartsimp_pg_sd(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                           S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn,
                           m_trees = opts$m_trees, d0 = d0)
    }, error = function(e) { cat(sprintf("  [d0=%d rep %d ERR] %s\n", d0, r, conditionMessage(e))); NULL })
    if (is.null(p_post)) next
    ad80 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.20)
    ad95 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.05)
    rmse <- sqrt(mean((ad80$p_a_post_mean - ad80$p_a_true)^2))
    per_rep[[length(per_rep)+1L]] <- data.frame(
      d0 = d0, rep = r, admin_cov80 = ad80$mean_cov, admin_cov95 = ad95$mean_cov,
      admin_rmse = rmse, admin_bias = mean(ad80$p_a_post_mean - ad80$p_a_true),
      admin_post_sd = mean(ad80$p_a_post_sd), admin_ratio = rmse / mean(ad80$p_a_post_sd),
      pixel_cov80 = compute_pixel_latent_cov(p_post, tr$p_star, 0.20),
      stringsAsFactors = FALSE)
  }
  pr <- do.call(rbind, per_rep)
  write.csv(pr, file.path(opts$out_dir, sprintf("dgp5_%s_per_rep.csv", opts$tag)), row.names = FALSE)
  cat(sprintf("[rep %2d/%d] done (%d d0 values)\n", r, opts$reps, length(d0_grid)))
}
per_rep <- do.call(rbind, per_rep)
mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
agg <- do.call(rbind, lapply(sort(unique(per_rep$d0)), function(dd) {
  d <- per_rep[per_rep$d0 == dd, ]
  data.frame(d0 = dd, n_rep = nrow(d),
    admin_cov80 = mean(d$admin_cov80), admin_cov80_se = mc_se(d$admin_cov80),
    admin_cov95 = mean(d$admin_cov95), admin_cov95_se = mc_se(d$admin_cov95),
    admin_bias = mean(d$admin_bias), admin_abs_bias = mean(abs(d$admin_bias)),
    admin_post_sd = mean(d$admin_post_sd), admin_ratio = mean(d$admin_ratio),
    admin_rmse = mean(d$admin_rmse),
    pixel_cov80 = mean(d$pixel_cov80), pixel_cov80_se = mc_se(d$pixel_cov80),
    stringsAsFactors = FALSE)
}))
num <- sapply(agg, is.numeric); agg[num] <- lapply(agg[num], round, 5)
write.csv(agg, file.path(opts$out_dir, sprintf("dgp5_%s_aggregated.csv", opts$tag)), row.names = FALSE)
cat("\n==== d0 SWEEP (sensitivity + ablation) ====\n"); print(agg, row.names = FALSE)
