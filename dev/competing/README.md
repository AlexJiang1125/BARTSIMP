# Competing methods (paper §4.2)

Wrappers expose a single common interface:

```r
fit_<method>(
  X_train, locs_train, n_i, k_i,
  X_pred,  locs_pred,
  n_iter = 1500, burn = 500,
  verbose = TRUE
)
```

Returns a list with:

| field            | shape           | meaning                                       |
|------------------|-----------------|-----------------------------------------------|
| `eta_pred_mean`  | length n_pred   | posterior mean of η at each pred point        |
| `eta_pred_q025`  | length n_pred   | 2.5%  quantile of η at each pred point        |
| `eta_pred_q975`  | length n_pred   | 97.5% quantile of η at each pred point        |
| `p_pred_mean`    | length n_pred   | posterior mean of p = sigma(η)                |
| `p_pred_q025`    | length n_pred   | 2.5%  quantile of p                           |
| `p_pred_q975`    | length n_pred   | 97.5% quantile of p                           |
| `eta_pred_draws` | n_keep × n_pred | per-iter draws (NULL for INLA-based methods)  |
| `elapsed_sec`    | scalar          | wall time                                     |
| `method`         | string          | method name                                   |

## Methods

- `bart_logistic.R`     — `BART::lbart`, no spatial component
- `spde_binomial.R`     — INLA `family="binomial"`, linear x1 + x2 + SPDE Matérn
- `spde0_binomial.R`    — INLA `family="binomial"`, intercept-only + SPDE Matérn
- `bartsimp_pg_wrap.R`  — same interface around `run_sampler_bartsimp_pg`

## V3 metrics computed downstream (one helper in `metrics.R`)

For each method given `(eta_pred_*, p_pred_*, eta_true, p_true, k_i, n_i)`:

- RMSE(η), RMSE(p)
- Brier score (over per-trial outcomes if observed at pred locs)
- log-loss (likewise)
- AIL(η), ACR(η), AIS(η) at α = 0.05
- ACR(p) on the grid (monotonic transform check)
