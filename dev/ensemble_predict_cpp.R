# =============================================================================
# ensemble_predict_cpp.R
#
# R-side loader for the C++ kernel in ensemble_predict_cpp.cpp.
# Compiles on first source(); falls back to the pure-R version from rbart.R
# if Rcpp isn't available or compilation fails.
#
# Use:
#   source("dev/rbart.R")                  # always — defines tree_predict, etc.
#   source("dev/ensemble_predict_cpp.R")   # optional — defines fast variants
#   # then call:
#   ensemble_predict(ens, X)              # routed via .ensemble_predict_impl
#
# The pure-R version remains available as `ensemble_predict_r()` for fallback
# testing.
# =============================================================================

# Stash the original pure-R implementation under a separate name so we can
# fall back to it deliberately (in tests, or if Rcpp is unavailable).
if (exists("ensemble_predict", inherits = FALSE) ||
    exists("ensemble_predict", mode = "function")) {
  ensemble_predict_r <- ensemble_predict
} else {
  warning("ensemble_predict() not found — load dev/rbart.R first.")
}

.bartsimp_cpp_state <- new.env(parent = emptyenv())
.bartsimp_cpp_state$loaded <- FALSE
.bartsimp_cpp_state$error  <- NULL

.try_load_cpp <- function() {
  if (.bartsimp_cpp_state$loaded) return(TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    .bartsimp_cpp_state$error <- "Rcpp not installed"
    return(FALSE)
  }
  # Locate the .cpp sibling of this file.
  this_file <- tryCatch(
    sys.frame(1)$ofile %||% normalizePath(sys.frames()[[1]]$ofile),
    error = function(e) NULL)
  cpp_path <- if (!is.null(this_file)) {
    file.path(dirname(normalizePath(this_file)), "ensemble_predict_cpp.cpp")
  } else {
    # Fallback heuristic: BARTSIMP_DEV env var, current dir, or dev/ subdir.
    bsdev <- Sys.getenv("BARTSIMP_DEV", unset = "")
    cands <- c(
      if (nzchar(bsdev)) file.path(bsdev, "ensemble_predict_cpp.cpp"),
      "ensemble_predict_cpp.cpp",
      "dev/ensemble_predict_cpp.cpp",
      "~/Documents/GitHub/BARTSIMP/dev/ensemble_predict_cpp.cpp")
    hit <- cands[file.exists(cands)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  if (length(cpp_path) == 0 || is.null(cpp_path) || is.na(cpp_path) ||
      !file.exists(cpp_path)) {
    .bartsimp_cpp_state$error <-
      paste0("ensemble_predict_cpp.cpp not found (looked near ",
             if (is.null(this_file)) "??" else this_file, ")")
    return(FALSE)
  }
  # R's default Makevars adds BLAS/LAPACK + gfortran linker flags even for
  # pure-C++ code. On Apple Silicon installs without /opt/gfortran this
  # breaks the build. Override with a temporary user-Makevars that strips
  # those flags — the kernel doesn't need them.
  tmp_makevars <- tempfile(pattern = "Makevars-")
  writeLines(c(
    "FLIBS=",
    "SHLIB_LIBADD="
  ), tmp_makevars)
  old_mv <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  Sys.setenv(R_MAKEVARS_USER = tmp_makevars)
  on.exit({
    if (is.na(old_mv)) Sys.unsetenv("R_MAKEVARS_USER")
    else Sys.setenv(R_MAKEVARS_USER = old_mv)
    try(file.remove(tmp_makevars), silent = TRUE)
  }, add = TRUE)
  res <- try(Rcpp::sourceCpp(cpp_path), silent = TRUE)
  if (inherits(res, "try-error")) {
    .bartsimp_cpp_state$error <- as.character(res)
    return(FALSE)
  }
  .bartsimp_cpp_state$loaded <- TRUE
  TRUE
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Public-facing routed version. Tries Rcpp; falls back to pure R.
ensemble_predict <- function(ensemble, X) {
  if (.try_load_cpp() && exists("ensemble_predict_cpp", mode = "function")) {
    as.numeric(ensemble_predict_cpp(ensemble, X))
  } else {
    ensemble_predict_r(ensemble, X)
  }
}

# Diagnostic helpers
ensemble_predict_backend <- function() {
  if (.try_load_cpp()) "cpp" else "r"
}
ensemble_predict_cpp_error <- function() .bartsimp_cpp_state$error
