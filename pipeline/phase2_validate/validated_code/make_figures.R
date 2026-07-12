#!/usr/bin/env Rscript
# make_figures.R — validation figures from the phase-2 CSVs.
#   Fig 1: paired Admin-1 cov80 (legacy vs SD) per rep + means, primary run.
#   Fig 2: d0 sensitivity/ablation curve (admin cov80 & pixel cov80 vs d0).
#   Fig 3: scale + adversarial cov80 (legacy vs SD) across conditions.
# Reads validated_results/*.csv; writes validated_results/figures/*.png.
RES <- "/Users/alexziyujiang/Documents/GitHub/BARTSIMP/pipeline/phase2_validate/validated_results"
FIG <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) { p <- file.path(RES, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }

# ---- Fig 1: primary paired cov80 ----
pr <- rd("dgp5_paired_primary_per_rep.csv")
if (!is.null(pr)) {
  bp <- pr[pr$method == "BARTSIMP-PG", ]
  sp <- pr[pr$method == "BARTSIMP-PG-SD", ]
  m <- merge(bp[, c("rep","admin_cov80")], sp[, c("rep","admin_cov80")],
             by = "rep", suffixes = c(".legacy",".SD"))
  png(file.path(FIG, "fig1_primary_paired_cov80.png"), width = 1100, height = 650, res = 130)
  par(mar = c(4.5,4.5,3,1))
  yl <- range(c(m$admin_cov80.legacy, m$admin_cov80.SD, 0.80), na.rm = TRUE)
  plot(m$rep, m$admin_cov80.legacy, type = "b", pch = 1, col = "grey40", ylim = yl,
       xlab = "replicate", ylab = "Admin-1 cov80",
       main = "Primary paired A/B: Admin-1 cov80 (legacy d0=0 vs SD)")
  lines(m$rep, m$admin_cov80.SD, type = "b", pch = 19, col = "firebrick")
  abline(h = 0.80, lty = 2, col = "blue")
  legend("bottomleft", c("BARTSIMP-PG (d0=0)","BARTSIMP-PG-SD","nominal 0.80"),
         pch = c(1,19,NA), lty = c(1,1,2), col = c("grey40","firebrick","blue"), bty = "n")
  dev.off()
  cat("wrote fig1\n")
}

# ---- Fig 2: d0 sweep ----
ds <- rd("dgp5_d0sweep_aggregated.csv")
if (!is.null(ds)) {
  ds <- ds[order(ds$d0), ]
  png(file.path(FIG, "fig2_d0_sweep.png"), width = 1100, height = 650, res = 130)
  par(mar = c(4.5,4.5,3,4.5))
  plot(ds$d0, ds$admin_cov80, type = "b", pch = 19, col = "firebrick",
       ylim = range(c(ds$admin_cov80, ds$admin_cov80 + ds$admin_cov80_se,
                      ds$admin_cov80 - ds$admin_cov80_se, 0.80), na.rm = TRUE),
       xlab = "depth-floor d0  (0 = legacy ablation)", ylab = "Admin-1 cov80",
       main = "d0 sensitivity + ablation")
  arrows(ds$d0, ds$admin_cov80 - ds$admin_cov80_se, ds$d0, ds$admin_cov80 + ds$admin_cov80_se,
         angle = 90, code = 3, length = 0.04, col = "firebrick")
  abline(h = 0.80, lty = 2, col = "blue")
  par(new = TRUE)
  plot(ds$d0, ds$pixel_cov80, type = "b", pch = 17, col = "darkgreen", axes = FALSE,
       xlab = "", ylab = "", ylim = range(c(ds$pixel_cov80, 0.90), na.rm = TRUE))
  axis(4, col = "darkgreen", col.axis = "darkgreen")
  mtext("pixel cov80 (guardrail)", side = 4, line = 3, col = "darkgreen")
  abline(h = 0.90, lty = 3, col = "darkgreen")
  legend("right", c("admin cov80","pixel cov80","nominal 0.80","pixel floor 0.90"),
         pch = c(19,17,NA,NA), lty = c(1,1,2,3),
         col = c("firebrick","darkgreen","blue","darkgreen"), bty = "n", cex = 0.85)
  dev.off()
  cat("wrote fig2\n")
}

# ---- Fig 3: scale + adversarial ----
plot_cond <- function(agg, xvar, xlab, main, fname) {
  if (is.null(agg)) return(invisible())
  lg <- agg[agg$arm == "legacy", ]; sd <- agg[agg$arm == "SD", ]
  lg <- lg[order(lg[[xvar]]), ]; sd <- sd[order(sd[[xvar]]), ]
  png(file.path(FIG, fname), width = 1100, height = 650, res = 130)
  par(mar = c(4.5,4.5,3,1))
  yl <- range(c(lg$admin_cov80, sd$admin_cov80, 0.80), na.rm = TRUE)
  plot(lg[[xvar]], lg$admin_cov80, type = "b", pch = 1, col = "grey40", ylim = yl,
       xlab = xlab, ylab = "Admin-1 cov80", main = main)
  if ("admin_cov80_se" %in% names(lg))
    arrows(lg[[xvar]], lg$admin_cov80-lg$admin_cov80_se, lg[[xvar]], lg$admin_cov80+lg$admin_cov80_se,
           angle=90, code=3, length=0.04, col="grey40")
  lines(sd[[xvar]], sd$admin_cov80, type = "b", pch = 19, col = "firebrick")
  if ("admin_cov80_se" %in% names(sd))
    arrows(sd[[xvar]], sd$admin_cov80-sd$admin_cov80_se, sd[[xvar]], sd$admin_cov80+sd$admin_cov80_se,
           angle=90, code=3, length=0.04, col="firebrick")
  abline(h = 0.80, lty = 2, col = "blue")
  legend("bottomleft", c("BARTSIMP-PG (d0=0)","BARTSIMP-PG-SD","nominal 0.80"),
         pch = c(1,19,NA), lty = c(1,1,2), col = c("grey40","firebrick","blue"), bty = "n")
  dev.off()
  cat("wrote", fname, "\n")
}
plot_cond(rd("dgp5_scale_aggregated.csv"), "A",
          "aggregation scale A (admin units)", "Scale sweep: cov80 vs A",
          "fig3a_scale.png")
plot_cond(rd("dgp5_adversarial_aggregated.csv"), "sigma_gp",
          "spatial amplitude sigma_gp", "Adversarial: cov80 vs spatial amplitude",
          "fig3b_adversarial.png")

cat("figures done\n")
