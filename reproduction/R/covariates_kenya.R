# =============================================================================
# Section 6 data access: assembling the REAL Kenya-DHS-like inputs
# =============================================================================
# This script documents (and provides runnable helpers for) how to rebuild the
# actual Section 6 inputs of
#
#   Jiang, A. Z., & Wakefield, J. (2025). BARTSIMP: Flexible spatial covariate
#   modeling and prediction using Bayesian Additive Regression Trees.
#   Spatial and Spatio-temporal Epidemiology, 55, 100757.
#
# from their PUBLIC sources. Nothing here downloads restricted data for you: the
# DHS microdata require a (free) registration and per-project approval, so you
# must obtain them yourself. Once you have the files on disk, the helpers below
# extract the six covariates at the cluster/grid locations and assemble the exact
# column layout that dgp_section6.R produces synthetically. After that, the rest
# of the pipeline -- cv_section6.R (stratified cluster CV) and areal.R
# (population-weighted areal aggregation) -- runs UNCHANGED on the real data.
#
# WHAT THE PAPER USES
#   Outcome  : child wasting, weight-for-height Z-score (WHZ), 2014 Kenya DHS.
#   Design   : 1584 enumeration areas (clusters); strata = county x urban/rural
#              (47 counties -> ~92 strata); ~25 children per cluster.
#   Areas    : Admin1 = 47 counties; Admin2 = 290 constituencies.
#   Weights  : under-five population density (WorldPop) for areal aggregation.
#
# THE SIX GEOSPATIAL COVARIATES  (names match SECTION6_COVARIATES in dgp_section6.R)
# ---------------------------------------------------------------------------
#   pop_density  Population density            WorldPop  (ppp, 100 m / 1 km)
#   night_light  Night-time lights             NOAA VIIRS / DMSP (annual composite)
#   ndvi         Vegetation index (NDVI)       MODIS MOD13 (250 m-1 km, 16-day)
#   temp         Mean annual temperature       WorldClim v2 (bio1, ~1 km)
#   precip       Annual precipitation          WorldClim v2 (bio12, ~1 km)
#   access_city  Travel time to nearest city   JRC / Malaria Atlas accessibility (~1 km)
# ---------------------------------------------------------------------------
# All six are rasters. We extract each at (a) the DHS cluster coordinates for the
# training/observation table and (b) the centroids of a regular prediction grid
# for the surface. The admin labels come from a boundary polygon set (GADM or the
# IEBC county/constituency shapefiles), and pop_u5 is WorldPop's under-5 count.
#
# This file is SOURCE-SAFE: it defines functions and prints guidance, but runs no
# extraction at source() time. Call the helpers explicitly once your files are in
# place. Every helper checks for the packages/paths it needs and errors with an
# actionable message rather than assuming the data exist.
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Where to get each input (URLs are stable landing pages, not direct files).
# ---------------------------------------------------------------------------
KENYA_DATA_SOURCES <- list(
  dhs = list(
    what = "2014 Kenya DHS: Children's Recode (KR) + Geographic (GE/GPS).",
    url  = "https://dhsprogram.com/data/",
    note = paste("Register a free DHS account, request the Kenya 2014 Standard",
                 "DHS, and download the KR (child anthropometry) .DTA/.SAS/.SPSS",
                 "recode and the GE (cluster GPS) shapefile. Cluster coordinates",
                 "are DELIBERATELY DISPLACED (urban <=2 km, rural <=5 km, 1% of",
                 "rural up to 10 km) -- keep this in mind for the SPDE mesh.")
  ),
  worldpop = list(
    what = "Population density + under-5 population counts (2014).",
    url  = "https://hub.worldpop.org/",
    note = "Kenya 'ppp' (people per pixel) 2014, and the age/sex structured u5 layer."
  ),
  night_light = list(
    what = "Night-time lights annual composite (2014).",
    url  = "https://www.ngdc.noaa.gov/eog/",
    note = "VIIRS annual VNL (2014+) or DMSP-OLS stable-lights for earlier years."
  ),
  ndvi = list(
    what = "MODIS NDVI (MOD13A2/MYD13A2), 2014 annual mean.",
    url  = "https://lpdaac.usgs.gov/products/mod13a2v061/",
    note = "Average the 16-day composites over 2014 to a single annual-mean raster."
  ),
  worldclim = list(
    what = "Mean annual temperature (bio1) and annual precipitation (bio12).",
    url  = "https://www.worldclim.org/data/worldclim21.html",
    note = "WorldClim v2.1 bioclimatic variables at 30s (~1 km)."
  ),
  accessibility = list(
    what = "Travel time to the nearest city / accessibility (2015).",
    url  = "https://malariaatlas.org/project-resources/accessibility-to-cities/",
    note = "Weiss et al. (2018) global accessibility; JRC GHSL is an alternative."
  ),
  boundaries = list(
    what = "Admin1 (47 counties) and Admin2 (290 constituencies) polygons.",
    url  = "https://gadm.org/download_country.html",
    note = "GADM level 1/2 for Kenya, or the official IEBC county/constituency shapefiles."
  )
)

# Print the data-access checklist (call interactively).
kenya_data_sources <- function() {
  for (nm in names(KENYA_DATA_SOURCES)) {
    s <- KENYA_DATA_SOURCES[[nm]]
    cat(sprintf("\n[%s]\n  what: %s\n  url : %s\n  note: %s\n",
                toupper(nm), s$what, s$url, s$note))
  }
  invisible(KENYA_DATA_SOURCES)
}

# ---------------------------------------------------------------------------
# 1. Extract the six covariate rasters at a set of point locations.
#
#   points     : data frame with longitude/latitude columns (WGS84, EPSG:4326).
#   raster_paths: NAMED list mapping each of SECTION6_COVARIATES to a raster file
#                 (GeoTIFF, .tif). Names MUST be the six covariate names.
#   lon_col,lat_col : coordinate column names in `points`.
#
# Returns `points` with six appended numeric columns (one per covariate), in the
# SECTION6_COVARIATES order. Uses terra; bilinear for continuous surfaces.
# ---------------------------------------------------------------------------
extract_covariates_at <- function(points, raster_paths,
                                   lon_col = "lon", lat_col = "lat",
                                   method = "bilinear") {
  need <- c("pop_density", "night_light", "ndvi", "temp", "precip", "access_city")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required to extract raster covariates.")
  miss <- setdiff(need, names(raster_paths))
  if (length(miss))
    stop("raster_paths is missing entries for: ", paste(miss, collapse = ", "),
         ".\nProvide a GeoTIFF path for each of: ", paste(need, collapse = ", "))
  bad <- raster_paths[need][!file.exists(unlist(raster_paths[need]))]
  if (length(bad))
    stop("These raster files do not exist:\n  ",
         paste(unlist(bad), collapse = "\n  "))

  xy <- as.matrix(points[, c(lon_col, lat_col)])
  for (nm in need) {
    r <- terra::rast(raster_paths[[nm]])
    v <- terra::extract(r, xy, method = method)
    points[[nm]] <- as.numeric(v[, 1])
  }
  # NDVI is a ratio in [-1, 1]; many products store it scaled by 1e4. If values
  # look scaled (|max| >> 1), rescale so downstream transforms behave.
  if (max(abs(points$ndvi), na.rm = TRUE) > 5) points$ndvi <- points$ndvi / 1e4
  points
}

# ---------------------------------------------------------------------------
# 2. Attach admin1 / admin2 labels by spatial join to a boundary polygon set.
#
#   points   : data frame with lon/lat columns.
#   admin_sf : an sf polygon object (e.g. GADM), one row per area at the level.
#   level_col: the column in admin_sf holding the area name/id to copy over.
#   out_col  : name for the new label column in `points` ("admin1"/"admin2").
#
# Returns `points` with an added `out_col`. Points falling outside all polygons
# (e.g. displaced offshore) get NA -- inspect and snap or drop those.
# ---------------------------------------------------------------------------
attach_admin_labels <- function(points, admin_sf, level_col,
                                out_col = "admin1",
                                lon_col = "lon", lat_col = "lat") {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for the admin spatial join.")
  pts <- sf::st_as_sf(points, coords = c(lon_col, lat_col), crs = 4326,
                      remove = FALSE)
  admin_sf <- sf::st_transform(admin_sf, 4326)
  ix <- as.integer(sf::st_within(pts, admin_sf, sparse = TRUE) |>
                     lapply(function(z) if (length(z)) z[1] else NA) |>
                     unlist())
  points[[out_col]] <- admin_sf[[level_col]][ix]
  points
}

# ---------------------------------------------------------------------------
# 3. Build the per-child observation table from the DHS Children's Recode.
#
#   kr : the DHS Children's Recode as a data frame (haven::read_dta on the .DTA).
#   ge : the DHS Geographic dataset (sf/data.frame from the GE shapefile), with
#        DHSCLUST, LATNUM, LONGNUM, URBAN_RURA (and optionally ADM1NAME).
#
# WHZ.  DHS ships a pre-computed WHO weight-for-height z-score in HW72, stored
# x100, with 9996-9999 as special/missing flags. We use it directly (this is
# what the paper's outcome is); alternatively recompute from weight/height/age/
# sex with the `anthro` package (WHO 2006 standards) -- see the commented block.
#
# Returns one row per child with: cluster, admin labels are added later by the
# caller via attach_admin_labels() on the cluster table; here we output cluster,
# urban, s1=LONGNUM, s2=LATNUM, and y=WHZ. Covariates are joined per-cluster
# (all children in a cluster share the cluster location), so extract them once on
# the cluster table (step 4) and merge.
# ---------------------------------------------------------------------------
build_dhs_children <- function(kr, ge,
                               clust_col = "v001",
                               whz_col   = "hw72",
                               weight_col = "v005") {
  if (!all(c(clust_col, whz_col) %in% names(kr)))
    stop("KR recode is missing '", clust_col, "' or '", whz_col,
         "'. Check the DHS variable names in the recode manual.")
  whz <- as.numeric(kr[[whz_col]]) / 100
  # DHS flags: 9996 (height out of plausible range), 9997 (age in months out of
  # range), 9998 (flagged), 9999 (missing) -- all become NA after the /100.
  whz[whz >= 99.9] <- NA_real_

  ge_df <- if (inherits(ge, "sf")) sf::st_drop_geometry(ge) else ge
  need_ge <- c("DHSCLUST", "LATNUM", "LONGNUM", "URBAN_RURA")
  miss <- setdiff(need_ge, names(ge_df))
  if (length(miss))
    stop("GE dataset is missing: ", paste(miss, collapse = ", "))
  # DHS marks unlocated clusters with (0,0) -- drop them.
  ge_df <- ge_df[!(ge_df$LATNUM == 0 & ge_df$LONGNUM == 0), , drop = FALSE]

  obs <- data.frame(
    cluster = as.integer(kr[[clust_col]]),
    whz     = whz,
    dhs_wt  = if (weight_col %in% names(kr)) as.numeric(kr[[weight_col]]) / 1e6 else 1
  )
  obs <- obs[is.finite(obs$whz), , drop = FALSE]
  loc <- data.frame(
    cluster = as.integer(ge_df$DHSCLUST),
    s1      = as.numeric(ge_df$LONGNUM),
    s2      = as.numeric(ge_df$LATNUM),
    urban   = as.integer(toupper(substr(ge_df$URBAN_RURA, 1, 1)) == "U")
  )
  out <- merge(obs, loc, by = "cluster")
  names(out)[names(out) == "whz"] <- "y"
  out[order(out$cluster), , drop = FALSE]

  # ---- ALTERNATIVE: recompute WHZ from raw measurements ---------------------
  # if (requireNamespace("anthro", quietly = TRUE)) {
  #   z <- anthro::anthro_zscores(
  #     sex = kr$b4, age = kr$hw1, is_age_in_month = TRUE,
  #     weight = kr$hw2 / 10, lenhei = kr$hw3 / 10,
  #     measure = ifelse(kr$hw15 == 1, "l", "h"))
  #   whz <- z$zwfl  # weight-for-length/height z-score (WHO 2006)
  # }
}

# ---------------------------------------------------------------------------
# 4. Build the cluster design table (one row per cluster) and attach covariates.
#
#   obs      : output of build_dhs_children() (per-child, with cluster/s1/s2/urban).
#   raster_paths, admin1_sf, admin2_sf : as in steps 1-2.
#
# Produces the cluster table with: cluster, s1, s2, urban, admin1, admin2,
# stratum (= admin1 x urban/rural, mirroring the paper), and the six covariates.
# Merging this back onto `obs` by `cluster` yields the exact per-child column
# layout that simulate_section6()$obs has -- ready for cv_section6.R.
# ---------------------------------------------------------------------------
build_dhs_clusters <- function(obs, raster_paths, admin1_sf, admin2_sf,
                               admin1_col, admin2_col) {
  first <- !duplicated(obs$cluster)
  cl <- obs[first, c("cluster", "s1", "s2", "urban"), drop = FALSE]
  # covariates: extract at cluster coordinates (lon = s1, lat = s2)
  cl <- extract_covariates_at(cl, raster_paths, lon_col = "s1", lat_col = "s2")
  cl <- attach_admin_labels(cl, admin1_sf, admin1_col, "admin1", "s1", "s2")
  cl <- attach_admin_labels(cl, admin2_sf, admin2_col, "admin2", "s1", "s2")
  cl$stratum <- paste0("a1_", cl$admin1, ifelse(cl$urban == 1L, "_U", "_R"))
  rownames(cl) <- NULL
  cl
}

# ---------------------------------------------------------------------------
# 5. Build the prediction grid + under-five population weight.
#
#   admin1_sf : boundary polygons (defines the country extent + admin1 labels).
#   raster_paths : the six covariate rasters.
#   pop_u5_path  : WorldPop under-5 count raster (for the areal weight).
#   n_side       : grid resolution (n_side x n_side cells over the bbox); or pass
#                  `template` (a terra SpatRaster) to predict on its cells.
#
# Returns a grid data frame matching simulate_section6()$grid's columns
# (cell, s1, s2, admin1, admin2, pop_u5, <six covariates>), minus the synthetic
# truth columns (z_star/g_cov/f_true), which do not exist for real data.
# ---------------------------------------------------------------------------
build_prediction_grid <- function(admin1_sf, admin2_sf, raster_paths,
                                   pop_u5_path, admin1_col, admin2_col,
                                   n_side = 100L) {
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE))
    stop("Packages 'terra' and 'sf' are required to build the prediction grid.")
  bb <- sf::st_bbox(sf::st_transform(admin1_sf, 4326))
  s1 <- seq(bb["xmin"], bb["xmax"], length.out = n_side)
  s2 <- seq(bb["ymin"], bb["ymax"], length.out = n_side)
  grid <- expand.grid(s1 = s1, s2 = s2)
  # keep only cells inside the country
  inside <- sf::st_within(
    sf::st_as_sf(grid, coords = c("s1", "s2"), crs = 4326, remove = FALSE),
    sf::st_union(sf::st_transform(admin1_sf, 4326)), sparse = FALSE)[, 1]
  grid <- grid[inside, , drop = FALSE]
  grid$cell <- seq_len(nrow(grid))

  grid <- extract_covariates_at(grid, raster_paths, lon_col = "s1", lat_col = "s2")
  grid <- attach_admin_labels(grid, admin1_sf, admin1_col, "admin1", "s1", "s2")
  grid <- attach_admin_labels(grid, admin2_sf, admin2_col, "admin2", "s1", "s2")

  if (!file.exists(pop_u5_path))
    stop("pop_u5 raster not found: ", pop_u5_path)
  pu5 <- terra::rast(pop_u5_path)
  grid$pop_u5 <- as.numeric(terra::extract(pu5, as.matrix(grid[, c("s1", "s2")]))[, 1])
  grid$pop_u5[!is.finite(grid$pop_u5) | grid$pop_u5 < 0] <- 0

  grid[, c("cell", "s1", "s2", "admin1", "admin2", "pop_u5",
           "pop_density", "night_light", "ndvi", "temp", "precip", "access_city")]
}

# ---------------------------------------------------------------------------
# END-TO-END RECIPE (pseudocode; fill in your local paths)
# ---------------------------------------------------------------------------
# library(haven); library(sf)
# kr  <- haven::read_dta("KEKR72DT/KEKR72FL.DTA")          # children's recode
# ge  <- sf::st_read("KEGE71FL/KEGE71FL.shp")              # cluster GPS
# a1  <- sf::st_read("gadm41_KEN_1.shp"); a2 <- sf::st_read("gadm41_KEN_2.shp")
# rp  <- list(pop_density="worldpop_ppp_2014.tif",
#             night_light="viirs_2014.tif", ndvi="modis_ndvi_2014.tif",
#             temp="wc2.1_bio1.tif", precip="wc2.1_bio12.tif",
#             access_city="accessibility_2015.tif")
#
# obs      <- build_dhs_children(kr, ge)
# clusters <- build_dhs_clusters(obs, rp, a1, a2, admin1_col="NAME_1", admin2_col="NAME_2")
# obs      <- merge(obs, clusters[, c("cluster","admin1","admin2","stratum",
#                                     SECTION6_COVARIATES)], by="cluster")
# grid     <- build_prediction_grid(a1, a2, rp, "worldpop_u5_2014.tif",
#                                   admin1_col="NAME_1", admin2_col="NAME_2", n_side=100)
#
# # From here the SYNTHETIC and REAL paths converge -- reuse Section 6 code as-is:
# sim_like <- list(obs = obs, grid = grid, clusters = clusters)
# res      <- run_cv_section6(sim_like, reps = 1:10, cov_names = SECTION6_COVARIATES)
# fit      <- fit_bartsimp(obs, grid, cov_names = SECTION6_COVARIATES, return_draws = TRUE)
# areas    <- aggregate_areal(fit$draws, grid, level = "admin1", weight = "pop_u5")
# =============================================================================
