## Installation and requirements

This package depends on the following R packages:

- **INLA** (for spatial modeling and mesh construction)
- **Rcpp** and **RcppArmadillo** (for compiled code support)

Because **INLA** is not available on CRAN, install it first:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

Then install **BARTSIMP** from GitHub:

```r
install.packages("remotes")
remotes::install_github("AlexJiang1125/BARTSIMP")
```


## Example

Below is a minimal example using the included toy dataset.

```r
library(dplyr)
library(BARTSIMP)

############################################
# Load example dataset
############################################
# The package includes a small simulated dataset
# with predictors, response, and spatial coordinates.
data("toy_data")

############################################
# Prepare model inputs
############################################

# Predictor matrix used by BART
x_all <- toy_data %>%
  select(x1, x2, x3, x4, x5) %>%
  data.frame()

# Response variable
y_all <- toy_data$y

# Spatial coordinates associated with each observation
s_all <- toy_data %>%
  select(s1, s2)

############################################
# Train / test split
############################################
# First 50 observations are used for training
# Remaining observations are used for testing
x_train <- x_all[1:50, ]
y_train <- y_all[1:50]
s_train <- s_all[1:50, ]

x_test <- x_all[51:75, ]
y_test <- y_all[51:75]
s_test <- s_all[51:75, ]

############################################
# Fit the BARTSIMP model
############################################
# The model first fits a BART regression to capture
# nonlinear relationships between predictors and response.
# A spatial random field will later be used to model
# residual spatial dependence.

fit <- bartsimp(
  x.train = x_train,
  y.train = y_train,
  x.test = x_test,
  s1 = s_train$s1,
  s2 = s_train$s2,
  s1.test = s_test$s1,
  s2.test = s_test$s2,
  ntree = 10,   # number of trees in the BART ensemble
  ndpost = 10,  # number of posterior samples to retain
  nwarmup = 10, # number of warmup / burn-in iterations
  seed = 1
)

############################################
# Spatial residual correction
############################################
# BART captures nonlinear covariate effects but does not
# explicitly model spatial correlation. We therefore model
# the residuals with a Matérn spatial Gaussian field and
# project the spatial effect to the test locations.

mod_predict <- bartsimp_predict(
  y_train = y_train,
  tree_draws_train = fit$yhat.train,
  tree_draws_test = fit$yhat.test,
  s_train = s_train,
  s_test = s_test,
  kappas = fit$kappas,
  sigmams = fit$sigmams,
  sigma_obs = fit$sigma,
  mesh = fit$mesh
)

############################################
# Posterior predictive means
############################################
# Matrix of predictions:
# rows  = posterior draws
# cols  = test observations
mod_predict$prediction_mean_test

# Posterior mean prediction
colMeans(mod_predict$prediction_mean_test)
# Root mean squared prediction error
sqrt(mean((colMeans(mod_predict$prediction_mean_test) - y_test)^2))
```
