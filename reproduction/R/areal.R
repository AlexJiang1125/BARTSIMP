# =============================================================================
# Section 6 helpers: stratified cluster CV + population-weighted areal estimates
# =============================================================================
# The application in Section 6 does two things beyond the Section 5 machinery:
#
#   1. Cross-validation that SPLITS BY CLUSTER (enumeration area) while
#      preserving the stratified design (paper: 1584 clusters -> 1267 train /
#      317 test, stratified, repeated 10x). Splitting by cluster (not by child)
#      is essential: children in the same cluster share a location, so a
#      row-wise split would leak spatial information into the test set.
#
#   2. Aggregation of point-level predictions to administrative areas by an
#      UNDER-FIVE POPULATION-WEIGHTED average (paper Eq. for WHZ_{R_i}):
#
#         WHZ_R = sum_{g in R} WHZ(s_g) d5(s_g) / sum_{g in R} d5(s_g).
#
#      The area's credible interval is computed from the AREAL posterior (apply
#      the weighted average to every posterior draw, then take quantiles), not
#      by combining cell-level intervals.
#
# These are generic: they take a cluster/design table and a grid with an area
# label + a weight column, so they work for the synthetic DGP and, unchanged,
# for the real Kenya data once you have assembled the same columns.
# =============================================================================

# Stratified split of CLUSTERS into train/test, preserving strata proportions.
#
#   clusters   : data frame with a cluster id column and a stratum column.
#   prop_train : fraction of each stratum's clusters assigned to training.
#   id_col, stratum_col : column names.
#   seed       : RNG seed.
#
# Returns list(train = <cluster ids>, test = <cluster ids>).
stratified_cluster_split <- function(clusters, prop_train = 0.8,
                                     id_col = "cluster",
                                     stratum_col = "stratum",
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ids   <- clusters[[id_col]]
  strat <- as.character(clusters[[stratum_col]])
  train <- unlist(lapply(split(ids, strat), function(g) {
    n_tr <- max(1L, min(length(g) - 1L, round(prop_train * length(g))))
    sample(g, n_tr)
  }), use.names = FALSE)
  # guard: if a stratum has a single cluster it is forced into train
  list(train = sort(train), test = sort(setdiff(ids, train)))
}

# Subset a per-child observation frame to a set of cluster ids.
subset_by_cluster <- function(obs, cluster_ids, id_col = "cluster") {
  obs[obs[[id_col]] %in% cluster_ids, , drop = FALSE]
}

# Population-weighted areal estimates from a cell-level posterior.
#
#   draws  : ndpost x |G| matrix of posterior draws of the mean field at the
#            grid cells (e.g. fit_bartsimp(..., return_draws = TRUE)$draws).
#   grid   : the grid data frame aligned with the COLUMNS of `draws`.
#   level  : name of the area-label column in `grid` (e.g. "admin1"/"admin2").
#   weight : name of the (positive) weight column, e.g. "pop_u5".
#   alpha  : central (1 - alpha) credible interval.
#
# Returns a data frame with one row per area: area, point (posterior mean),
# lower, upper (posterior quantiles), n_cells.
aggregate_areal <- function(draws, grid, level = "admin1",
                            weight = "pop_u5", alpha = 0.05) {
  stopifnot(ncol(draws) == nrow(grid))
  areas <- sort(unique(grid[[level]]))
  qs <- c(alpha / 2, 1 - alpha / 2)
  rows <- lapply(areas, function(a) {
    cols <- which(grid[[level]] == a)
    w <- grid[[weight]][cols]
    wn <- w / sum(w)
    areal_draws <- as.numeric(draws[, cols, drop = FALSE] %*% wn)  # length ndpost
    qq <- stats::quantile(areal_draws, probs = qs, names = FALSE)
    data.frame(area = a, point = mean(areal_draws),
               lower = qq[1], upper = qq[2], n_cells = length(cols))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Direct (design-agnostic) areal estimate for comparison, computed from the
# observations only. To respect clustering we average CLUSTER means within each
# area and form a normal interval from the between-cluster spread. This is a
# simplified stand-in for the paper's survey-weighted direct estimate (which
# uses the DHS design weights); with real data, pass design-weighted cluster
# means and a design-based variance instead.
direct_areal_estimate <- function(obs, level = "admin1", alpha = 0.05,
                                   id_col = "cluster") {
  z <- stats::qnorm(1 - alpha / 2)
  areas <- sort(unique(obs[[level]]))
  rows <- lapply(areas, function(a) {
    sub <- obs[obs[[level]] == a, , drop = FALSE]
    cl_means <- tapply(sub$y, sub[[id_col]], mean)
    k <- length(cl_means)
    m <- mean(cl_means)
    se <- if (k > 1) stats::sd(cl_means) / sqrt(k) else NA_real_
    data.frame(area = a, point = m,
               lower = m - z * se, upper = m + z * se, n_clusters = k)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# True area-level WHZ (population-weighted average of the true surface), for
# validating areal predictions against the synthetic ground truth.
true_areal <- function(grid, level = "admin1", weight = "pop_u5",
                       truth = "f_true") {
  areas <- sort(unique(grid[[level]]))
  val <- vapply(areas, function(a) {
    cols <- grid[[level]] == a
    stats::weighted.mean(grid[[truth]][cols], grid[[weight]][cols])
  }, numeric(1))
  data.frame(area = areas, truth = val)
}
