#' Preprocess Test Data for BARTSIMP Spatial Prediction
#'
#' Validates and formats test covariates and spatial coordinates for use in
#' BARTSIMP prediction routines. The function ensures that the covariate matrix
#' and spatial coordinate matrix have consistent dimensions and contain no
#' missing values.
#'
#' This utility is primarily used internally by prediction functions to enforce
#' a consistent representation of test inputs.
#'
#' @param x_test Matrix or data frame of test covariates. Each row corresponds
#'   to a test observation.
#'
#' @param s_test Optional matrix or data frame containing spatial coordinates
#'   for the test observations. Must contain exactly two columns representing
#'   spatial coordinates (e.g., longitude and latitude).
#'
#' @details
#' The function performs the following checks:
#'
#' * Converts `x_test` and `s_test` to data frames.
#' * Verifies that `s_test` has exactly two columns.
#' * Ensures the number of rows in `x_test` and `s_test` match.
#' * Ensures neither input contains missing values.
#'
#' If `x_test` is empty (zero rows or zero columns), the function returns empty
#' objects for both covariates and spatial coordinates to allow downstream
#' prediction code to proceed safely.
#'
#' @return A list containing:
#'
#' \describe{
#' \item{x_test}{
#'   A data frame of processed test covariates.
#' }
#'
#' \item{s_test}{
#'   A data frame with columns `s1` and `s2` representing spatial coordinates.
#' }
#' }
#'
#' @examples
#' x_test <- data.frame(x1 = rnorm(5), x2 = runif(5))
#' s_test <- data.frame(lon = runif(5), lat = runif(5))
#'
#' BARTSIMP:::preprocess_test_data(x_test, s_test)
#'
#' @keywords internal
preprocess_test_data <- function(x_test, s_test = NULL) {
  x_test <- as.data.frame(x_test)

  # allow empty x.test
  if (nrow(x_test) == 0 || ncol(x_test) == 0) {
    return(list(
      x_test = matrix(0.0, 0, 0),
      s_test = data.frame(s1 = numeric(0), s2 = numeric(0))
    ))
  }

  if (is.null(s_test)) {
    stop("For non-empty x.test, test spatial coordinates must also be provided.")
  }

  s_test <- as.data.frame(s_test)
  if (ncol(s_test) != 2) {
    stop("s.test must have exactly 2 columns.")
  }
  colnames(s_test) <- c("s1", "s2")

  if (nrow(x_test) != nrow(s_test)) {
    stop("x.test and s.test must have matching numbers of rows.")
  }
  if (anyNA(x_test) || anyNA(s_test)) {
    stop("x.test and s.test must not contain missing values.")
  }

  list(
    x_test = x_test,
    s_test = s_test
  )
}
