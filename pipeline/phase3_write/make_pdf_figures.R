#!/usr/bin/env Rscript
# make_pdf_figures.R — publication PDF figures from the validated CSVs.
#   fig2_d0_sweep.pdf       : HEADLINE calibration evidence (admin cov80 + pixel guardrail + ratio)
#   fig1_primary_paired.pdf : primary paired Admin-1 cov80, legacy vs SD, per rep + means
#   fig3a_scale.pdf         : scale sweep cov80 vs A (legacy vs SD)
#   fig3b_adversarial.pdf   : adversarial recentering cov80 vs sigma_gp (toward 0.80 from both sides)
DATA <- "/Users/alexziyujiang/Documents/GitHub/BARTSIMP/pipeline/phase3_write/data"
FIG  <- "/Users/alexziyujiang/Documents/GitHub/BARTSIMP/pipeline/phase3_write/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path(DATA, f), stringsAsFactors = FALSE)

CL_LEG <- "grey35"; CL_SD <- "firebrick"; CL_NOM <- "royalblue3"; CL_PIX <- "darkgreen"

# ---- Fig 2: d0 sweep (HEADLINE) ----
ds <- rd("dgp5_d0sweep_aggregated.csv"); ds <- ds[order(ds$d0), ]
pdf(file.path(FIG, "fig2_d0_sweep.pdf"), width = 7.2, height = 4.6, pointsize = 11)
par(mar = c(4.4, 4.4, 2.4, 4.6), las = 1)
yl <- range(c(ds$admin_cov80 + ds$admin_cov80_se, ds$admin_cov80 - ds$admin_cov80_se, 0.79, 0.84))
plot(ds$d0, ds$admin_cov80, type = "b", pch = 19, col = CL_SD, lwd = 2, ylim = yl,
     xlab = expression("depth-floor "*d[0]*"   (0 = legacy BARTSIMP-PG ablation)"),
     ylab = "Admin-1 cov80", main = expression("Calibration vs depth-floor "*d[0]))
arrows(ds$d0, ds$admin_cov80 - ds$admin_cov80_se, ds$d0, ds$admin_cov80 + ds$admin_cov80_se,
       angle = 90, code = 3, length = 0.04, col = CL_SD)
abline(h = 0.80, lty = 2, col = CL_NOM, lwd = 1.5)
par(new = TRUE)
plot(ds$d0, ds$pixel_cov80, type = "b", pch = 17, col = CL_PIX, lwd = 1.6, axes = FALSE,
     xlab = "", ylab = "", ylim = c(0.84, 0.91))
axis(4, col = CL_PIX, col.axis = CL_PIX)
mtext("pixel cov80 (fine-scale guardrail)", side = 4, line = 3, col = CL_PIX, las = 0)
legend("topleft", c("Admin-1 cov80 (±1 MC SE)", "pixel cov80", "nominal 0.80"),
       pch = c(19, 17, NA), lty = c(1, 1, 2), lwd = c(2, 1.6, 1.5),
       col = c(CL_SD, CL_PIX, CL_NOM), bty = "n", cex = 0.85)
dev.off(); cat("wrote fig2_d0_sweep.pdf\n")

# ---- Fig 1: primary paired ----
pr <- rd("dgp5_paired_primary_per_rep.csv")
bp <- pr[pr$method == "BARTSIMP-PG", c("rep", "admin_cov80")]
sp <- pr[pr$method == "BARTSIMP-PG-SD", c("rep", "admin_cov80")]
m <- merge(bp, sp, by = "rep", suffixes = c(".legacy", ".SD"))
pdf(file.path(FIG, "fig1_primary_paired.pdf"), width = 7.2, height = 4.6, pointsize = 11)
par(mar = c(4.4, 4.4, 2.4, 1), las = 1)
yl <- range(c(m$admin_cov80.legacy, m$admin_cov80.SD, 0.80), na.rm = TRUE)
plot(m$rep, m$admin_cov80.legacy, type = "b", pch = 1, col = CL_LEG, ylim = yl,
     xlab = "replicate (matched-seed paired design)", ylab = "Admin-1 cov80",
     main = "Primary paired A/B (24 reps, N=500, A=20)")
lines(m$rep, m$admin_cov80.SD, type = "b", pch = 19, col = CL_SD)
abline(h = 0.80, lty = 2, col = CL_NOM, lwd = 1.5)
mL <- mean(m$admin_cov80.legacy); mS <- mean(m$admin_cov80.SD)
abline(h = mL, lty = 3, col = CL_LEG); abline(h = mS, lty = 3, col = CL_SD)
legend("bottomleft",
       c(sprintf("BARTSIMP-PG d0=0  (mean %.3f)", mL),
         sprintf("BARTSIMP-PG-SD    (mean %.3f)", mS), "nominal 0.80"),
       pch = c(1, 19, NA), lty = c(1, 1, 2), col = c(CL_LEG, CL_SD, CL_NOM), bty = "n", cex = 0.82)
dev.off(); cat("wrote fig1_primary_paired.pdf\n")

# ---- Fig 3a/3b: scale + adversarial ----
plot_cond <- function(agg, xvar, xlab, main, fname, xexpr = NULL) {
  lg <- agg[agg$arm == "legacy", ]; sd <- agg[agg$arm == "SD", ]
  lg <- lg[order(lg[[xvar]]), ]; sd <- sd[order(sd[[xvar]]), ]
  pdf(file.path(FIG, fname), width = 6.4, height = 4.6, pointsize = 11)
  par(mar = c(4.4, 4.4, 2.4, 1), las = 1)
  yl <- range(c(lg$admin_cov80 + lg$admin_cov80_se, lg$admin_cov80 - lg$admin_cov80_se,
                sd$admin_cov80 + sd$admin_cov80_se, sd$admin_cov80 - sd$admin_cov80_se, 0.80))
  plot(lg[[xvar]], lg$admin_cov80, type = "b", pch = 1, col = CL_LEG, ylim = yl,
       xlab = if (is.null(xexpr)) xlab else xexpr, ylab = "Admin-1 cov80", main = main)
  arrows(lg[[xvar]], lg$admin_cov80 - lg$admin_cov80_se, lg[[xvar]], lg$admin_cov80 + lg$admin_cov80_se,
         angle = 90, code = 3, length = 0.04, col = CL_LEG)
  lines(sd[[xvar]], sd$admin_cov80, type = "b", pch = 19, col = CL_SD)
  arrows(sd[[xvar]], sd$admin_cov80 - sd$admin_cov80_se, sd[[xvar]], sd$admin_cov80 + sd$admin_cov80_se,
         angle = 90, code = 3, length = 0.04, col = CL_SD)
  abline(h = 0.80, lty = 2, col = CL_NOM, lwd = 1.5)
  legend("bottomleft", c("BARTSIMP-PG (d0=0)", "BARTSIMP-PG-SD", "nominal 0.80"),
         pch = c(1, 19, NA), lty = c(1, 1, 2), col = c(CL_LEG, CL_SD, CL_NOM), bty = "n", cex = 0.82)
  dev.off(); cat("wrote", fname, "\n")
}
plot_cond(rd("dgp5_scale_aggregated.csv"), "A", "aggregation scale A (no. admin units)",
          "Scale sweep", "fig3a_scale.pdf")
plot_cond(rd("dgp5_adversarial_aggregated.csv"), "sigma_gp", "spatial amplitude sigma_gp",
          "Adversarial: two-sided recentering", "fig3b_adversarial.pdf",
          xexpr = expression("spatial amplitude "*sigma[gp]))
cat("ALL PDF FIGURES DONE\n")
