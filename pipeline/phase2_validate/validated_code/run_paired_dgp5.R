#!/usr/bin/env Rscript
# run_paired_dgp5.R — PRIMARY validation: paired A/B of the semi-dense depth-floor.
#
#   BARTSIMP-PG     (d0 = 0, legacy)        <-- comparator (the arm to beat)
#   BARTSIMP-PG-SD  (d0 = ceil(sqrt(log N))) <-- proposed
#
# plus the core-8 roster carried on the SAME matched-seed designs. Matched seeds:
# truth drawn once, data drawn per rep with set.seed(base_seed + r); each method
# fit under its own fixed sampler seed so the paired delta isolates the d0 knob.
#
# Writes per-rep + aggregated CSVs (method, metric, value, mc_se) to
# pipeline/phase2_validate/validated_results/.
#
# RUN (Framework R, from anywhere; the script setwd()s to the cluster dir):
#   env -u MallocNanoZone /Library/Frameworks/R.framework/Resources/bin/Rscript \
#       run_paired_dgp5.R --reps 30 --A 20
suppressPackageStartupMessages({ library(optparse) })

opts <- parse_args(OptionParser(option_list = list(
  make_option("--reps",    default = 30L,  type = "integer"),
  make_option("--N",       default = 500L, type = "integer"),
  make_option("--A",       default = 20L,  type = "integer"),
  make_option("--n_draws", default = 500L, type = "integer"),
  make_option("--n_iter",  default = 600L, type = "integer"),
  make_option("--burn",    default = 200L, type = "integer"),
  make_option("--m_trees", default = 20L,  type = "integer"),
  make_option("--d0",      default = -1L,  type = "integer",
              help = "depth-floor for SD arm; -1 => ceil(sqrt(log N))"),
  make_option("--methods", default = "GLM,GLMM-iCAR,GLMM-SPDE,Vanilla-BART,BARTSIMP-PG,BARTSIMP-PG-SD",
              help = "comma-separated roster"),
  make_option("--base_seed",  default = 90000L, type = "integer"),
  make_option("--seed_truth", default = 4242L,  type = "integer"),
  make_option("--alpha0",     default = -0.90,  type = "double"),
  make_option("--tag",        default = "primary"),
  make_option("--out_dir",
              default = "/Users/alexziyujiang/Documents/GitHub/BARTSIMP/pipeline/phase2_validate/validated_results")
)))

CLUSTER <- "/Users/alexziyujiang/Documents/Claude/Projects/BARTSIMP_glm/bartsimp_cluster"
setwd(CLUSTER)
b07 <- "experiments/sim_thm_verification/sim_07/R"
source(file.path(b07, "dgp.R"))
source(file.path(b07, "methods_inla.R"))
source(file.path(b07, "methods_bartsimp.R"))
source(file.path(b07, "coverage.R"))
source("experiments/sim_dgp5_coverage/R/spatial_ml_competitors.R")
.carf <- "experiments/sim_dgp5_coverage/R/vendor/CAR-Forest.R"
if (file.exists(.carf)) source(.carf)
source("experiments/sim_dgp5_coverage/R/design_culture_competitors.R")  # fit_malasso_gen, fit_ppi_gen
source("experiments/sim_dgp5_coverage/R/dgp5_core.R")

dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
d0 <- if (opts$d0 < 0L) as.integer(ceiling(sqrt(log(opts$N)))) else opts$d0
methods <- strsplit(opts$methods, ",")[[1]]
cat(sprintf("PAIRED dgp5: N=%d A=%d reps=%d d0(SD)=%d methods=%s\n",
            opts$N, opts$A, opts$reps, d0, paste(methods, collapse=",")))

truth        <- make_truth5(alpha_0 = opts$alpha0, seed_truth = opts$seed_truth)
admin_struct <- make_admin_structure(A = opts$A, seed = 1L)
mfn <- dgp5_method_fn()

# Latent PIXEL coverage guardrail: fraction of clusters whose (1-alpha) credible
# interval on the LATENT p_i (posterior quantiles of p_post[,i]) covers the TRUE
# p_i^* surface value. This is the fine-scale calibration the depth-floor must
# NOT break (spec: pixel cov80 must stay in the conservative band, >= 0.90).
# Distinct from compute_cluster_cov (which scores the binomial posterior-
# predictive K/n against the noisy observed k/n -> the discreteness floor).
compute_pixel_latent_cov <- function(p_post, p_star, alpha = 0.20) {
  N <- ncol(p_post); cov_i <- numeric(N)
  for (i in seq_len(N)) {
    iv <- quantile(p_post[, i], probs = c(alpha/2, 1 - alpha/2), names = FALSE, na.rm = TRUE)
    cov_i[i] <- as.integer(p_star[i] >= iv[1] & p_star[i] <= iv[2])
  }
  mean(cov_i)
}

# Score one method's (S x N) p_post into a one-row data.frame.
score_method <- function(p_post, design, admin_id, tr, m, rep_id, el, d0_used) {
  ad80 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.20)
  ad95 <- compute_admin_cov(p_post, design$n_i, admin_id, tr$p_star, alpha = 0.05)
  cl80 <- tryCatch(compute_cluster_cov(p_post, design$n_i, tr$k_i, alpha = 0.20),
                   error = function(e) list(mean_cov = NA_real_))
  pix80 <- tryCatch(compute_pixel_latent_cov(p_post, tr$p_star, alpha = 0.20),
                    error = function(e) NA_real_)
  rmse <- sqrt(mean((ad80$p_a_post_mean - ad80$p_a_true)^2))
  bias <- mean(ad80$p_a_post_mean - ad80$p_a_true)
  psd  <- mean(ad80$p_a_post_sd)
  data.frame(
    method = m, rep = rep_id, d0 = d0_used,
    admin_cov80 = ad80$mean_cov, admin_cov95 = ad95$mean_cov,
    admin_rmse = rmse, admin_bias = bias, admin_post_sd = psd,
    admin_ratio = rmse / psd,
    cluster_cov80 = if (!is.null(cl80$mean_cov)) cl80$mean_cov else NA_real_,
    pixel_cov80 = pix80,
    elapsed_sec = el, stringsAsFactors = FALSE)
}

per_rep <- list()
for (r in seq_len(opts$reps)) {
  set.seed(opts$base_seed + r)
  design   <- draw_design5(opts$N)
  admin_id <- assign_admin(design$locs, admin_struct)
  tr       <- apply_truth5(design, truth)

  for (m in methods) {
    t0 <- Sys.time()
    p_post <- tryCatch({
      if (m == "BARTSIMP-PG-SD") {
        set.seed(opts$base_seed * 7L + r)   # fixed sampler seed, matched to legacy
        fit_bartsimp_pg_sd(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                           S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn,
                           m_trees = opts$m_trees, d0 = d0)
      } else if (m == "BARTSIMP-PG") {
        set.seed(opts$base_seed * 7L + r)   # SAME sampler seed as SD arm (paired)
        fit_bartsimp_pg(design$X, design$locs, design$n_i, tr$k_i, admin_id,
                        S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn,
                        m_trees = opts$m_trees)
      } else if (m %in% c("BARTSIMP-PG-soft", "BARTSIMP-BB-RE")) {
        mfn[[m]](design$X, design$locs, design$n_i, tr$k_i, admin_id,
                 S = opts$n_draws, n_iter = opts$n_iter, burn = opts$burn,
                 m_trees = opts$m_trees)
      } else {
        mfn[[m]](design$X, design$locs, design$n_i, tr$k_i, admin_id, S = opts$n_draws)
      }
    }, error = function(e) {
      cat(sprintf("  [%s rep %d ERROR] %s\n", m, r, conditionMessage(e))); NULL
    })
    el <- as.numeric(Sys.time() - t0, units = "secs")
    if (is.null(p_post)) next
    d0_used <- if (m == "BARTSIMP-PG-SD") d0 else 0L
    row <- score_method(p_post, design, admin_id, tr, m, r, el, d0_used)
    per_rep[[length(per_rep) + 1L]] <- row
  }
  # checkpoint after each rep
  pr <- do.call(rbind, per_rep)
  write.csv(pr, file.path(opts$out_dir, sprintf("dgp5_paired_%s_per_rep.csv", opts$tag)),
            row.names = FALSE)
  bp <- pr[pr$method == "BARTSIMP-PG"    & pr$rep == r, ]
  sp <- pr[pr$method == "BARTSIMP-PG-SD" & pr$rep == r, ]
  cat(sprintf("[rep %2d/%d] legacy cov80=%.2f bias=%+.4f | SD cov80=%.2f bias=%+.4f\n",
              r, opts$reps,
              if(nrow(bp)) bp$admin_cov80 else NA, if(nrow(bp)) bp$admin_bias else NA,
              if(nrow(sp)) sp$admin_cov80 else NA, if(nrow(sp)) sp$admin_bias else NA))
}

per_rep <- do.call(rbind, per_rep)
write.csv(per_rep, file.path(opts$out_dir, sprintf("dgp5_paired_%s_per_rep.csv", opts$tag)),
          row.names = FALSE)

# Aggregate: mean +/- MC SE per (method, metric), long format.
mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
metrics <- c("admin_cov80","admin_cov95","admin_rmse","admin_bias",
             "admin_post_sd","admin_ratio","cluster_cov80","pixel_cov80","elapsed_sec")
agg <- do.call(rbind, lapply(unique(per_rep$method), function(m) {
  d <- per_rep[per_rep$method == m, , drop = FALSE]
  do.call(rbind, lapply(metrics, function(mt) data.frame(
    method = m, metric = mt, value = mean(d[[mt]], na.rm = TRUE),
    mc_se = mc_se(d[[mt]]), n_rep = sum(is.finite(d[[mt]])),
    stringsAsFactors = FALSE)))
}))
agg$value <- round(agg$value, 5); agg$mc_se <- round(agg$mc_se, 5)
write.csv(agg, file.path(opts$out_dir, sprintf("dgp5_paired_%s_aggregated.csv", opts$tag)),
          row.names = FALSE)

# Paired within-rep deltas (SD - legacy), the cleanest signal.
wide <- reshape(per_rep[per_rep$method %in% c("BARTSIMP-PG","BARTSIMP-PG-SD"),
                        c("rep","method","admin_cov80","admin_bias","admin_post_sd","admin_rmse","admin_cov95")],
                idvar = "rep", timevar = "method", direction = "wide")
pair_ok <- complete.cases(wide)
wide <- wide[pair_ok, ]
dcov80 <- wide$`admin_cov80.BARTSIMP-PG-SD` - wide$`admin_cov80.BARTSIMP-PG`
dbias  <- abs(wide$`admin_bias.BARTSIMP-PG-SD`) - abs(wide$`admin_bias.BARTSIMP-PG`)
dpsd   <- wide$`admin_post_sd.BARTSIMP-PG-SD` - wide$`admin_post_sd.BARTSIMP-PG`
paired <- data.frame(
  quantity = c("delta_cov80","delta_abs_bias","delta_post_sd",
               "frac_reps_bias_reduced","frac_reps_cov80_up"),
  value = c(mean(dcov80), mean(dbias), mean(dpsd),
            mean(dbias < 0), mean(dcov80 > 0)),
  mc_se = c(mc_se(dcov80), mc_se(dbias), mc_se(dpsd), NA, NA),
  n_pair = sum(pair_ok))
paired$value <- round(paired$value, 5); paired$mc_se <- round(paired$mc_se, 5)
write.csv(paired, file.path(opts$out_dir, sprintf("dgp5_paired_%s_deltas.csv", opts$tag)),
          row.names = FALSE)

cat("\n==== AGGREGATED (admin cov80 / cov95 / ratio) ====\n")
show <- agg[agg$metric %in% c("admin_cov80","admin_cov95","admin_ratio","admin_bias","cluster_cov80"), ]
print(show[order(show$metric, show$method), ], row.names = FALSE)
cat("\n==== PAIRED DELTAS (SD - legacy) ====\n")
print(paired, row.names = FALSE)
