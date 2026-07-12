# BARTSIMP paper reproduction materials

Reproducible code for the simulation study (**Section 5**) and the application
(**Section 6**) of

> Jiang, A. Z., & Wakefield, J. (2025). *BARTSIMP: Flexible spatial covariate
> modeling and prediction using Bayesian Additive Regression Trees.* Spatial and
> Spatio-temporal Epidemiology, **55**, 100757.
> <https://doi.org/10.1016/j.sste.2025.100757>

Everything here runs on **self-contained synthetic data** that mirrors the
structure of the paper's inputs, so you can exercise the full workflow without
any restricted downloads. Section 6 additionally documents how to assemble the
*real* Kenya DHS inputs (see `R/covariates_kenya.R`); once assembled in the same
column layout, the same code runs unchanged.

## Requirements

* **R** with a working **INLA** install (`install.packages("INLA", repos =
  c(INLA = "https://inla.r-inla-download.org/R/stable"))`). INLA ships compiled
  code, so it must match your R build's ABI; if `library(INLA)` segfaults, you
  are almost certainly mixing R builds (e.g. a Homebrew/conda R against an INLA
  built for CRAN R). Use one consistent R installation.
* **dbarts** (the BART competitor and the partial-dependence demo).
* The **BARTSIMP** package itself, installed from this repository
  (`R CMD INSTALL .` from the repo root, or `devtools::install()`).
* For the *real* Section 6 data only: **terra**, **sf**, **haven**.

`library(BARTSIMP)` attaches INLA automatically (it is in `Depends`).

## What's here

Top level:

| File                     | Purpose                                                    |
|--------------------------|------------------------------------------------------------|
| `vignette_section5.Rmd`  | Section 5 walkthrough: one scenario, a sweep, and the fast Woodbury predictor. |
| `run_section5.R`         | Headless Section 5 driver (writes CSVs to `results/`).     |
| `vignette_section6.Rmd`  | Section 6 walkthrough: cluster CV, areal aggregation, partial dependence. |

`R/` building blocks (sourcing `R/run_scenario.R` pulls in the Section 5 set):

| File                     | Purpose                                                    |
|--------------------------|------------------------------------------------------------|
| `dgp_section5.R`         | Section 5 data-generating process (Matérn GRF + covariate surfaces). |
| `metrics.R`              | RMSE, AIL, ACR, AIS (interval scoring).                    |
| `competitors.R`          | BART (`dbarts`), SPDE, SPDE0 (`INLA`) fitters.             |
| `bartsimp_fit.R`         | BARTSIMP two-stage fitter (sampler + Woodbury predictor).  |
| `run_scenario.R`         | One-scenario runner; sources the four files above.         |
| `dgp_section6.R`         | Section 6 synthetic Kenya-DHS-like DGP.                    |
| `areal.R`                | Stratified cluster split + population-weighted areal estimates. |
| `cv_section6.R`          | Section 6 stratified cross-validation runner.              |
| `covariates_kenya.R`     | How to assemble the **real** DHS + raster inputs (data-access + extraction). |

(Scripts prefixed with `_` are developer smoke/validation checks, not part of the
public workflow.)

## Section 5 quick start

Knit the vignette, or run the driver:

```r
rmarkdown::render("vignette_section5.Rmd")
```

```sh
Rscript run_section5.R            # quick demo config
Rscript run_section5.R --paper    # full paper-scale sweep (long)
```

## Section 6 quick start

```r
rmarkdown::render("vignette_section6.Rmd")
```

or interactively:

```r
library(BARTSIMP)
.repro_R <- "R"
source(file.path(.repro_R, "run_scenario.R"))
source(file.path(.repro_R, "dgp_section6.R"))
source(file.path(.repro_R, "areal.R"))
source(file.path(.repro_R, "cv_section6.R"))

sim <- simulate_section6(seed = 1)
cv  <- run_cv_section6(sim, reps = 1:2)          # cluster CV, four methods
fit <- fit_bartsimp(sim$obs, sim$grid, cov_names = SECTION6_COVARIATES,
                    return_draws = TRUE)          # areal needs cell-level draws
areas <- aggregate_areal(fit$draws, sim$grid, level = "admin1", weight = "pop_u5")
```

To run on the real 2014 Kenya DHS data, follow `R/covariates_kenya.R`
(`kenya_data_sources()` prints the download checklist) to build `obs`, `grid`,
and `clusters` with the same columns, then reuse the calls above.

## The fast Woodbury predictor

Stage 2 of BARTSIMP (recovering the residual Matérn field per posterior draw)
originally refit INLA per draw. With the hyperparameters fixed at each draw, that
fit is a Gaussian-Markov-random-field posterior mean available in closed form, so
`BARTSIMP::bartsimp_predict(..., method = "woodbury")` computes it with one sparse
Cholesky --- numerically identical to `method = "inla"` but far cheaper. The
reproduction fitter `fit_bartsimp()` uses `"woodbury"` by default. See the
"fast spatial predictor" section of the Section 5 vignette for a side-by-side
equivalence-and-timing check.
