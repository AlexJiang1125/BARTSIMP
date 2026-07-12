# Modeling Briefing (paper-modeler-owned)

## Evaluations carried into the manuscript

All four reuse the Phase-2 validated CSVs (`validated_results/` → copied to
`data/`); the validation report's verdict (CONDITIONAL GO) is honored. Per the
budget-discipline instruction, the validated CSVs are the empirical core (a
complete, coherent story: primary + sweep + scale + adversarial).

| Evaluation | Harness | Design | Reps | CSV |
|------------|---------|--------|------|-----|
| Primary paired A/B | `sim_dgp5_coverage` | N=500, A=20, matched seeds, light MCMC (m_trees=20, n_iter=600, burn=200) | 24 | `dgp5_paired_primary_{aggregated,per_rep,deltas}.csv` |
| Sensitivity / ablation (d₀ sweep) | `sim_dgp5_coverage` | d₀∈{0,1,2,3,4}, d₀=0 = legacy ablation | 16 | `dgp5_d0sweep_{aggregated,per_rep}.csv` |
| Scale sweep | `sim_dgp5_coverage` | A∈{5,20,40}, paired | 12 | `dgp5_scale_{aggregated,per_rep}.csv` |
| Adversarial (spatial amplitude) | `sim_dgp5_coverage` | σ_gp∈{0.5,1.0,1.5}, paired | 12 | `dgp5_adversarial_{aggregated,per_rep}.csv` |

Known-ground-truth: yes — `sim_dgp5_coverage` scores against the known
5-covariate nonlinear non-additive truth + Matérn field; the admin-1
area-average prevalence is the scored estimand with the realized-cluster
(in-sample) weighting matching `compute_admin_cov`.

MC standard errors: present for every Monte-Carlo metric (SE = sd/√B), in the
aggregated CSVs (`mc_se` / `_se` columns).

## Comparator coverage confirmation (`comparator_inventory.md`, gate COMPLETE)

| FULL ID | Primary table | Scale/roster context |
|---------|---------------|----------------------|
| BARTSIMP-PG-SD (proposed) | YES (full rep) | YES |
| C1 BARTSIMP-PG (d₀=0, arm-to-beat) | YES (full rep) | YES |
| C2 GLM | YES | YES |
| C3 Vanilla-BART | YES | YES |
| C4 GLMM-iCAR | YES | YES |
| C5 GLMM-SPDE | YES | YES |
| C6 disaggregation | roster-context (reduced reps) | YES |
| C7 RF-GLS | roster-context (≤5 reps, ~480 s/rep) | context only |
| C8 SPARforest | roster-context (≤5 reps) | context only |

The primary paired loop carries proposed + C1–C5 at full rep; C6–C8 are
compute-throttled roster context (documented in `comparator_inventory.md`, with
the d₀-gate arm BARTSIMP-PG d₀=0 at full rep). The Results tables carry all
FULL IDs; reduced-rep comparators are flagged with the caveat (anomaly A6).

## Production figures (PDF, in `figures/`)

- `fig2_d0_sweep.pdf` — HEADLINE calibration evidence (admin cov80 + ratio +
  pixel guardrail vs d₀).
- `fig1_primary_paired.pdf` — primary paired Admin-1 cov80, legacy vs SD.
- `fig3a_scale.pdf` — scale sweep cov80 vs A.
- `fig3b_adversarial.pdf` — adversarial two-sided recentering vs σ_gp.

Regenerated as PDF from the validated CSVs by `make_pdf_figures.R` (Framework R).

## Production scale-up (best-effort, non-blocking)

A 60-rep paired production run was launched DETACHED (`run_paired_dgp5.R --reps
60 --tag production_scaleup`, base_seed 90000 matched to validation). At
~5 min/rep it does not complete within budget; the manuscript uses the validated
24-rep numbers and states the planned production N (≥60–100) honestly as future
work to resolve the primary contrast. Partial output → `scaleups/`.

## Formula↔code audit: see `formula_code_audit.md` (zero unresolved MISMATCH).
## Anomaly log: see `anomaly_log.md` (no Type-1; A1/A2 are honesty-critical).
