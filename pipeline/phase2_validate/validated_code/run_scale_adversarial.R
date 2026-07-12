#!/usr/bin/env Rscript
# run_scale_adversarial.R — SCALE + ADVERSARIAL stress tests, paired legacy vs SD.
#
# SCALE: aggregation-scale sweep A in {5,20,40} (the "calibrated at every scale"
#        property; the spec requires ratio in [0.95,1.10] at every A, and cov80
#        improvement to hold across A).
# ADVERSARIAL: spatial-amplitude crank sigma_gp in {0.5,1.0,1.5} (model stress:
#        a stronger spatial field is harder for the tree mean to leave alone;
#        does the depth-floor degrade MORE or LESS than legacy as the field grows?)
#
# Each (condition) runs the paired BARTSIMP-PG (d0=0) vs BARTSIMP-PG-SD on matched
# seeds. Reduced reps (stress-test budget).
suppressPackageStartupMessages({ library(optparse) })
opts <- parse_args(OptionParser(option_list = list(
  make_option("--reps",    default = 12L,  type = "integer"),
  make_option("--N",       default = 500L, type = "integer"),
  make_option("--n_draws", default = 500L, type = "integer"),
  make_option("--n_iter",  default = 600L, type = "integer"),
  make_option("--burn",    default = 200L, type = "integer"),
  make_option("--m_trees", default = 20L,  type = "integer"),
  make_option("--mode",    default = "scale", help = "scale | adversarial"),
  make_option("--base_seed",  default = 90000L, type = "integer"),
  make_option("--seed_truth", default = 4242L,  type = "integer"),
  make_option("--alpha0",     default = -0.90,  type = "double"),
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

# truth with a tunable spatial-field amplitude (adversarial knob). sigma_gp=0.5
# reproduces make_truth5 exactly.
make_truth5_sigma <- function(alpha_0, seed_truth, sigma_gp = 0.5, N_field = 8000L) {
  set.seed(seed_truth)
  locs_dense <- cbind(s1 = runif(N_field), s2 = runif(N_field))
  g_field <- draw_matern_field(locs_dense, sigma = sigma_gp, kappa = 4, seed = seed_truth + 1L)
  list(alpha_0 = alpha_0, f_truth = f_truth5, g_field = g_field, name = "DGP5")
}
N <- opts$N
d0 <- as.integer(ceiling(sqrt(log(N))))

if (opts$mode == "scale") {
  conditions <- data.frame(A = c(5L,20L,40L), sigma_gp = 0.5)
} else {
  conditions <- data.frame(A = 20L, sigma_gp = c(0.5,1.0,1.5))
}

compute_pixel_latent_cov <- function(p_post, p_star, alpha = 0.20) {
  N <- ncol(p_post); cov_i <- numeric(N)
  for (i in seq_len(N)) {
    iv <- quantile(p_post[, i], probs = c(alpha/2, 1 - alpha/2), names = FALSE, na.rm = TRUE)
    cov_i[i] <- as.integer(p_star[i] >= iv[1] & p_star[i] <= iv[2])
  }
  mean(cov_i)
}
score <- function(p_post, design, admin_id, tr) {
  ad80 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.20)
  ad95 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.05)
  rmse <- sqrt(mean((ad80$p_a_post_mean - ad80$p_a_true)^2))
  list(cov80 = ad80$mean_cov, cov95 = ad95$mean_cov, rmse = rmse,
       bias = mean(ad80$p_a_post_mean - ad80$p_a_true),
       post_sd = mean(ad80$p_a_post_sd), ratio = rmse / mean(ad80$p_a_post_sd),
       pixel80 = compute_pixel_latent_cov(p_post, tr$p_star, 0.20))
}

per_rep <- list()
for (ci in seq_len(nrow(conditions))) {
  A <- conditions$A[ci]; sg <- conditions$sigma_gp[ci]
  truth <- make_truth5_sigma(opts$alpha0, opts$seed_truth, sigma_gp = sg)
  admin_struct <- make_admin_structure(A = A, seed = 1L)
  cat(sprintf("== condition %d/%d: A=%d sigma_gp=%.1f d0=%d ==\n", ci, nrow(conditions), A, sg, d0))
  for (r in seq_len(opts$reps)) {
    set.seed(opts$base_seed + r)
    design   <- draw_design5(N)
    admin_id <- assign_admin(design$locs, admin_struct)
    tr       <- apply_truth5(design, truth)
    for (arm in c("legacy","SD")) {
      set.seed(opts$base_seed * 7L + r)   # matched sampler seed across arms
      p_post <- tryCatch({
        if (arm == "legacy")
          fit_bartsimp_pg(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                          S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn, m_trees = opts$m_trees)
        else
          fit_bartsimp_pg_sd(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                             S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn,
                             m_trees = opts$m_trees, d0 = d0)
      }, error = function(e) { cat(sprintf("   [%s A=%d sg=%.1f r%d ERR] %s\n", arm, A, sg, r, conditionMessage(e))); NULL })
      if (is.null(p_post)) next
      s <- score(p_post, design, admin_id, tr)
      per_rep[[length(per_rep)+1L]] <- data.frame(
        mode = opts$mode, A = A, sigma_gp = sg, arm = arm,
        d0 = if (arm=="SD") d0 else 0L, rep = r,
        admin_cov80 = s$cov80, admin_cov95 = s$cov95, admin_rmse = s$rmse,
        admin_bias = s$bias, admin_post_sd = s$post_sd, admin_ratio = s$ratio,
        pixel_cov80 = s$pixel80, stringsAsFactors = FALSE)
    }
    pr <- do.call(rbind, per_rep)
    write.csv(pr, file.path(opts$out_dir, sprintf("dgp5_%s_per_rep.csv", opts$mode)), row.names = FALSE)
  }
}
per_rep <- do.call(rbind, per_rep)
mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
agg <- do.call(rbind, by(per_rep, list(per_rep$A, per_rep$sigma_gp, per_rep$arm), function(d) {
  data.frame(mode = d$mode[1], A = d$A[1], sigma_gp = d$sigma_gp[1], arm = d$arm[1],
    n_rep = nrow(d),
    admin_cov80 = mean(d$admin_cov80), admin_cov80_se = mc_se(d$admin_cov80),
    admin_cov95 = mean(d$admin_cov95), admin_bias = mean(d$admin_bias),
    admin_abs_bias = mean(abs(d$admin_bias)),
    admin_post_sd = mean(d$admin_post_sd), admin_ratio = mean(d$admin_ratio),
    pixel_cov80 = mean(d$pixel_cov80), stringsAsFactors = FALSE)
}))
num <- sapply(agg, is.numeric); agg[num] <- lapply(agg[num], round, 5)
agg <- agg[order(agg$A, agg$sigma_gp, agg$arm), ]
write.csv(agg, file.path(opts$out_dir, sprintf("dgp5_%s_aggregated.csv", opts$mode)), row.names = FALSE)
cat(sprintf("\n==== %s stress (paired legacy vs SD) ====\n", toupper(opts$mode)))
print(agg, row.names = FALSE)
