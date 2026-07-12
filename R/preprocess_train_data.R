#' Preprocess Training Data for BARTSIMP Spatial Models
#'
#' Validates and formats training covariates, responses, and spatial
#' coordinates prior to model fitting. The function ensures consistent
#' dimensions, removes missing values, and sorts observations by spatial
#' location so that repeated locations appear in contiguous blocks.
#'
#' This ordering is required by downstream routines that assume training
#' observations from the same spatial location are grouped together.
#'
#' @param x_train Matrix or data frame of training covariates. Each row
#'   corresponds to a training observation.
#'
#' @param y_train Numeric vector of responses corresponding to the rows
#'   of `x_train`.
#'
#' @param s1 Numeric vector of the first spatial coordinate for each
#'   observation.
#'
#' @param s2 Numeric vector of the second spatial coordinate for each
#'   observation.
#'
#' @param size Optional integer vector specifying the number of
#'   observations at each unique spatial location. If provided, the
#'   function checks whether it matches the grouping implied by the
#'   spatial coordinates. If inconsistent, a warning is issued and the
#'   grouping is recomputed.
#'
#' @details
#' The function performs several preprocessing steps:
#'
#' 1. Converts `x_train` to a data frame and `y_train` to a numeric vector.
#' 2. Verifies that all inputs have matching lengths.
#' 3. Ensures no missing values are present.
#' 4. Sorts observations by `(s1, s2)` so that identical spatial locations
#'    are contiguous.
#' 5. Computes cluster sizes corresponding to repeated spatial coordinates.
#'
#' Cluster sizes are calculated using run-length encoding of the sorted
#' coordinate pairs. These sizes represent the number of observations at
#' each unique spatial location.
#'
#' @return A list containing:
#'
#' \describe{
#' \item{x_train}{
#'   The sorted training covariate data frame.
#' }
#'
#' \item{y_train}{
#'   The reordered response vector.
#' }
#'
#' \item{s1}{
#'   Sorted first spatial coordinate.
#' }
#'
#' \item{s2}{
#'   Sorted second spatial coordinate.
#' }
#'
#' \item{size}{
#'   Integer vector containing the number of observations at each unique
#'   spatial location.
#' }
#'
#' \item{train_order}{
#'   The permutation applied to the original data when sorting by spatial
#'   location.
#' }
#' }
#'
#' @examples
#' x_train <- data.frame(x1 = rnorm(5), x2 = runif(5))
#' y_train <- rnorm(5)
#' s1 <- c(1, 1, 2, 2, 2)
#' s2 <- c(3, 3, 4, 4, 4)
#'
#' BARTSIMP:::preprocess_train_data(x_train, y_train, s1, s2)
#'
#' @keywords internal
preprocess_train_data <- function(x_train, y_train, s1, s2, size = NULL) {
  x_train <- as.data.frame(x_train)
  y_train <- as.numeric(y_train)
  s_train <- data.frame(s1 = s1, s2 = s2)

  if (nrow(x_train) != length(y_train)) {
    stop("x.train and y.train must have matching numbers of rows.")
  }
  if (nrow(x_train) != nrow(s_train)) {
    stop("x.train, y.train, s1, and s2 must have matching numbers of rows.")
  }
  if (anyNA(x_train) || anyNA(y_train) || anyNA(s_train)) {
    stop("x.train, y.train, s1, and s2 must not contain missing values.")
  }

  # sort training rows by spatial location so repeated locations are contiguous
  ord <- order(s_train$s1, s_train$s2)
  x_train <- x_train[ord, , drop = FALSE]
  y_train <- y_train[ord]
  s_train <- s_train[ord, , drop = FALSE]

  # compute cluster sizes from sorted spatial coordinates
  loc_key <- paste(s_train$s1, s_train$s2, sep = "\r")
  size_computed <- rle(loc_key)$lengths

  # if user supplied size, check consistency; otherwise overwrite with computed size
  if (!is.null(size)) {
    size <- as.integer(size)
    if (sum(size) != nrow(x_train)) {
      warning("Supplied 'size' is inconsistent with training rows; recomputing size from (s1, s2).")
      size <- size_computed
    } else if (!identical(as.integer(size), as.integer(size_computed))) {
      warning("Supplied 'size' does not match grouped training coordinates; recomputing size from (s1, s2).")
      size <- size_computed
    }
  } else {
    size <- size_computed
  }

  list(
    x_train = x_train,
    y_train = y_train,
    s1 = s_train$s1,
    s2 = s_train$s2,
    size = size,
    train_order = ord
  )
}
