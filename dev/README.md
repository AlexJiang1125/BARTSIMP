# `dev/` — GLM extension prototypes

Scratch directory for the binary/GLM extension of BARTSIMP. Code here is
not part of the package (excluded via `.Rbuildignore`).

## Stage 1 roadmap (explicit-z PG Gibbs sampler)

The goal of Stage 1 is a pure-Gibbs sampler that maintains the spatial GP
explicitly, using Pólya–Gamma augmentation. Built up in three steps so
each piece can be validated before the next is added:

- `01_pg_matern_intercept.R` — PG + Matérn GP, intercept only.
  Validates the augmentation + spatial conditional math without any
  covariate or SPDE complications. Recovery target: `(beta_0, sigma_m, rho)`.

- `02_pg_matern_linear.R` — adds a linear covariate term.
  Adds a Gibbs step for `beta` and checks linear coefficient recovery.

- `03_pg_matern_bart.R` — replaces the linear predictor with a BART
  sum-of-trees. This is the Stage 1 target.

Once Stage 1 is solid we move to Stage 2 (integrate `z` out via INLA),
which is roughly "BARTSIMP-Gaussian applied to PG pseudo-data."
