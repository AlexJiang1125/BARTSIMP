#' Spatially Corrected Prediction for BARTSIMP Models
#'
#' Computes posterior predictive means for test locations by combining
#' BART posterior tree predictions with a spatial residual model using
#' a Matérn Gaussian random field fitted via INLA/SPDE.
#'
#' The function assumes that posterior draws from a BART model have already
#' been obtained for both training and test locations. For each posterior draw,
#' the residuals between the observed outcomes and the BART predictions at the
#' training locations are modeled using a spatial Matérn field. The fitted
#' spatial field is then projected to the test locations and added to the
#' BART predictions to obtain spatially adjusted predictions.
#'
#' @param y_train Numeric vector of observed outcomes at the training locations.
#'
#' @param tree_draws_train Numeric matrix of posterior BART predictions at the
#'   training locations. Rows correspond to posterior draws and columns
#'   correspond to training observations.
#'
#' @param tree_draws_test Numeric matrix of posterior BART predictions at the
#'   test locations. Rows correspond to posterior draws and columns correspond
#'   to test observations.
#'
#' @param s_train Matrix or data frame of training spatial coordinates.
#'   Must contain exactly two columns corresponding to spatial coordinates
#'   (e.g., longitude and latitude).
#'
#' @param s_test Matrix or data frame of test spatial coordinates.
#'   Must contain exactly two columns corresponding to spatial coordinates.
#'
#' @param kappas Numeric vector of Matérn range parameters for each posterior
#'   draw. Must have length equal to the number of posterior draws.
#'
#' @param sigmams Numeric vector of marginal standard deviations of the spatial
#'   field for each posterior draw. Must have length equal to the number of
#'   posterior draws.
#'
#' @param sigma_obs Optional numeric vector specifying the observation noise
#'   standard deviation for each posterior draw. If supplied, the observation
#'   precision is fixed in the INLA model. If `NULL`, INLA estimates the noise
#'   level.
#'
#' @param mesh An `INLA` mesh object defining the SPDE approximation to the
#'   spatial Gaussian field. Typically created using
#'   `INLA::inla.mesh.2d()`.
#'
#' @param nu Shape parameter used when fixing the observation precision prior
#'   (only relevant when `sigma_obs` is provided).
#'
#' @param lambda Rate parameter used when fixing the observation precision prior
#'   (only relevant when `sigma_obs` is provided).
#'
#' @param verbose Logical; if `TRUE`, progress messages are printed during
#'   posterior processing.
#'
#' @details
#' The model assumes the following decomposition for the outcome:
#'
#' \deqn{
#' y(s) = f_{\text{BART}}(x) + w(s) + \epsilon
#' }
#'
#' where
#'
#' * \eqn{f_{\text{BART}}(x)} is the prediction from the BART model,
#' * \eqn{w(s)} is a spatial Gaussian random field with a Matérn covariance
#'   structure,
#' * \eqn{\epsilon} is observational noise.
#'
#' For each posterior draw:
#'
#' 1. Residuals are computed:
#' \deqn{r = y_{\text{train}} - f_{\text{BART}}}
#'
#' 2. A spatial Matérn field is fitted to the residuals using the SPDE
#'    representation implemented in `INLA`.
#'
#' 3. The fitted spatial field is projected to the test locations using the
#'    mesh projection matrix.
#'
#' 4. The spatial prediction is added to the BART prediction to obtain the final
#'    posterior predictive mean.
#'
#' This approach allows the BART model to capture nonlinear feature effects
#' while the Gaussian field accounts for residual spatial dependence.
#'
#' @return A list with the following components:
#'
#' \describe{
#' \item{spatial_mean_test}{
#'   Matrix of posterior spatial predictions at the test locations.
#'   Rows correspond to posterior draws and columns correspond to test points.
#' }
#'
#' \item{prediction_mean_test}{
#'   Matrix of posterior predictive means at the test locations after combining
#'   BART predictions and spatial corrections.
#' }
#'
#' \item{s_train}{
#'   Processed training coordinates used internally by the function.
#' }
#'
#' \item{s_test}{
#'   Processed test coordinates used internally by the function.
#' }
#' }
#'
#' @examples
#' \dontrun{
#' result <- bartsimp_predict(
#'   y_train = y,
#'   tree_draws_train = bart_train_draws,
#'   tree_draws_test = bart_test_draws,
#'   s_train = coords_train,
#'   s_test = coords_test,
#'   kappas = kappas,
#'   sigmams = sigmams,
#'   mesh = mesh
#' )
#'
#' pred_mean <- colMeans(result$prediction_mean_test)
#' }
#'
#' @references
#' Lindgren, F., Rue, H., & Lindström, J. (2011).
#' An explicit link between Gaussian fields and Gaussian Markov random fields:
#' the stochastic partial differential equation approach.
#' *Journal of the Royal Statistical Society: Series B*, 73(4), 423–498.
#'
#' @export
bartsimp_predict <- function(y_train,
                             tree_draws_train,
                             tree_draws_test,
                             s_train,
                             s_test,
                             kappas,
                             sigmams,
                             sigma_obs = NULL,
                             mesh,
                             nu = 1,
                             lambda = sqrt(8) / 2.5,
                             verbose = FALSE) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required.")
  }

  # -----------------------------
  # basic coercion
  # -----------------------------
  tree_draws_train <- as.matrix(tree_draws_train)
  tree_draws_test  <- as.matrix(tree_draws_test)

  s_train <- as.data.frame(s_train)
  s_test  <- as.data.frame(s_test)

  if (ncol(s_train) != 2) {
    stop("s_train must have exactly 2 columns.")
  }
  if (ncol(s_test) != 2) {
    stop("s_test must have exactly 2 columns.")
  }

  colnames(s_train) <- c("s1", "s2")
  colnames(s_test)  <- c("s1", "s2")

  # -----------------------------
  # preprocess training data
  # sort by location so rows are grouped consistently
  # -----------------------------
  train_dat <- preprocess_train_data(
    x_train = data.frame(dummy = seq_along(y_train)),
    y_train = y_train,
    s1 = s_train$s1,
    s2 = s_train$s2,
    size = NULL
  )

  # recover sorted training order and apply it to tree draws
  ord_train <- order(s_train$s1, s_train$s2)

  y_train <- train_dat$y_train
  s_train <- data.frame(s1 = train_dat$s1, s2 = train_dat$s2)

  tree_draws_train <- tree_draws_train[, ord_train, drop = FALSE]

  # -----------------------------
  # preprocess test data
  # preserve row order
  # -----------------------------
  test_dat <- preprocess_test_data(
    x_test = data.frame(dummy = seq_len(nrow(s_test))),
    s_test = s_test
  )

  s_test <- data.frame(
    s1 = test_dat$s_test$s1,
    s2 = test_dat$s_test$s2
  )

  # -----------------------------
  # dimensions / consistency checks
  # -----------------------------
  ndpost <- nrow(tree_draws_train)
  n_train <- ncol(tree_draws_train)
  n_test <- nrow(s_test)

  if (length(y_train) != n_train) {
    stop("length(y_train) must equal ncol(tree_draws_train).")
  }
  if (nrow(tree_draws_test) != ndpost) {
    stop("nrow(tree_draws_test) must equal nrow(tree_draws_train).")
  }
  if (ncol(tree_draws_test) != n_test) {
    stop("ncol(tree_draws_test) must equal nrow(s_test).")
  }
  if (length(kappas) != ndpost) {
    stop("length(kappas) must equal nrow(tree_draws_train).")
  }
  if (length(sigmams) != ndpost) {
    stop("length(sigmams) must equal nrow(tree_draws_train).")
  }
  if (!is.null(sigma_obs) && length(sigma_obs) != ndpost) {
    stop("If provided, length(sigma_obs) must equal nrow(tree_draws_train).")
  }
  if (is.null(mesh)) {
    stop("mesh must be supplied.")
  }

  spatial_mean_test <- matrix(NA_real_, nrow = ndpost, ncol = n_test)

  # Projection matrices do not depend on draw
  A_train <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(s_train))
  A_test  <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(s_test))

  for (m in seq_len(ndpost)) {
    if (verbose && (m %% 10 == 0 || m == 1)) {
      message("Processing draw ", m, " / ", ndpost)
    }

    resid_m <- y_train - tree_draws_train[m, ]

    spde_m <- INLA::inla.spde2.pcmatern(
      mesh = mesh,
      prior.range = c(sqrt(8) / kappas[m], NA),
      prior.sigma = c(sigmams[m], NA)
    )

    stack_est <- INLA::inla.stack(
      data = list(y = resid_m),
      A = list(A_train),
      effects = list(
        i = seq_len(spde_m$n.spde)
      ),
      tag = "est"
    )

    stack_pred <- INLA::inla.stack(
      data = list(y = NA),
      A = list(A_test),
      effects = list(
        i = seq_len(spde_m$n.spde)
      ),
      tag = "pred"
    )

    stack_full <- INLA::inla.stack(stack_est, stack_pred)

    if (is.null(sigma_obs)) {
      result_m <- INLA::inla(
        y ~ -1 + f(i, model = spde_m),
        data = INLA::inla.stack.data(stack_full),
        control.predictor = list(
          A = INLA::inla.stack.A(stack_full),
          compute = TRUE
        ),
        control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
        verbose = FALSE
      )
    } else {
      hyper_fixed <- list(
        prec = list(
          prior = "loggamma",
          fixed = TRUE,
          initial = -2 * log(sigma_obs[m]),
          param = c(nu / 2, lambda * nu / 2)
        )
      )

      result_m <- INLA::inla(
        y ~ -1 + f(i, model = spde_m),
        data = INLA::inla.stack.data(stack_full),
        family = "gaussian",
        control.family = list(hyper = hyper_fixed),
        control.predictor = list(
          A = INLA::inla.stack.A(stack_full),
          compute = TRUE
        ),
        control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
        verbose = FALSE
      )
    }

    idx_pred <- INLA::inla.stack.index(stack_full, tag = "pred")$data
    spatial_mean_test[m, ] <- result_m$summary.fitted.values[idx_pred, "mean"]
  }

  list(
    spatial_mean_test = spatial_mean_test,
    prediction_mean_test = spatial_mean_test + tree_draws_test,
    s_train = s_train,
    s_test = s_test
  )
}
