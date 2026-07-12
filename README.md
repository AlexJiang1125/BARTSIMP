# BARTSIMP

**BART for Spatial INLA Modeling and Prediction.**

BARTSIMP models a continuous outcome as a flexible covariate mean plus a spatial
field: a Bayesian Additive Regression Trees (BART) ensemble captures non-linear
covariate effects, and a Matérn Gaussian random field (via the SPDE / INLA
representation) captures residual spatial dependence,

$$ y_{ik} = \mathrm{BART}(x_i) + z(s_i) + \varepsilon_{ik}. $$

It is designed for geostatistical problems where covariate effects are non-linear
and predictions (and areal aggregates) need calibrated uncertainty. The package
structure and interface are inspired by the **BART** R package.

Method and applications are described in:

> Jiang, A. Z., & Wakefield, J. (2025). *BARTSIMP: Flexible spatial covariate
> modeling and prediction using Bayesian Additive Regression Trees.* Spatial and
> Spatio-temporal Epidemiology, **55**, 100757.
> <https://doi.org/10.1016/j.sste.2025.100757>

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

INLA ships compiled code, so it must match your R build. If `library(INLA)`
crashes, you are likely mixing R installations (e.g. a Homebrew/conda R with an
INLA built for CRAN R) --- use one consistent R. `library(BARTSIMP)` attaches
INLA automatically (it is in `Depends`).

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

### Faster prediction (Woodbury)

The spatial residual correction above originally refit INLA once per posterior
draw. With the Matérn hyperparameters fixed at each draw, that fit is a
Gaussian-Markov-random-field posterior mean available in closed form, so it can
be computed with a single sparse Cholesky. Pass `method = "woodbury"` to use it:

```r
mod_predict <- bartsimp_predict(
  y_train = y_train,
  tree_draws_train = fit$yhat.train,
  tree_draws_test = fit$yhat.test,
  s_train = s_train, s_test = s_test,
  kappas = fit$kappas, sigmams = fit$sigmams,
  sigma_obs = fit$sigma, mesh = fit$mesh,
  method = "woodbury"          # closed-form GMRF mean; needs sigma_obs
)
```

It is numerically identical to the default `method = "inla"` (agreement ~1e-8)
but orders of magnitude faster, and the gap grows with the number of draws, the
mesh size, and the number of prediction points.

## Reproducing the paper (Sections 5 & 6)

The [`reproduction/`](reproduction/) folder contains self-contained code for the
paper's simulation study (Section 5) and application (Section 6), with two knitr
vignettes you can run without any restricted downloads:

- **Section 5** (`reproduction/vignette_section5.Rmd`) --- the gridded simulation
  comparing BARTSIMP against BART / SPDE / SPDE0 on RMSE, AIL, ACR, and AIS, plus
  a side-by-side check of the Woodbury predictor.
- **Section 6** (`reproduction/vignette_section6.Rmd`) --- the Kenya-DHS-style
  application: stratified cluster cross-validation, population-weighted areal
  aggregation, and partial dependence, on a synthetic dataset that mirrors the
  real design.

The DHS microdata are not redistributable, so the Section 6 vignette runs on
synthetic data by default; `reproduction/R/covariates_kenya.R` documents how to
obtain the real 2014 Kenya DHS outcome and the six public geospatial covariates
(WorldPop, NOAA, MODIS, WorldClim, JRC) and assemble them into the same layout.
See [`reproduction/README.md`](reproduction/README.md) for the full guide.

## Citation

If you use BARTSIMP, please cite:

```bibtex
@article{jiang2025bartsimp,
  title   = {{BARTSIMP}: Flexible spatial covariate modeling and prediction
             using {B}ayesian Additive Regression Trees},
  author  = {Jiang, Alex Ziyu and Wakefield, Jon},
  journal = {Spatial and Spatio-temporal Epidemiology},
  volume  = {55},
  pages   = {100757},
  year    = {2025},
  doi     = {10.1016/j.sste.2025.100757}
}
```

## License

MIT (see [LICENSE](LICENSE)).
