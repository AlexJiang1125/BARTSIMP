# Iteration Log — Phase 2 Validate: BARTSIMP-PG-SD (semi-dense depth-floor)

Append-only chronological log. Method: BARTSIMP-PG-SD = BARTSIMP-PG (BART +
SPDE–Matérn GP + Pólya–Gamma binomial-logit) + hard coarse-layer depth-floor d₀
on the trees. Reduces to legacy BARTSIMP-PG at d₀=0. Gate comparator: BARTSIMP-PG
(d₀=0), paired rep-for-rep.

## Iteration 0 — SETUP
- Read `pipeline/phase1_think/methodology_specification.md`. Extracted: core method
  (depth-floor as coverage-recentering knob), parameter d₀ ∈ {0,1,2,3,4} with
  default ⌈√(log N)⌉ ≈ 2–3, gate comparator BARTSIMP-PG (d₀=0), success criterion
  (admin cov80 → 0.80 without breaking pixel calibration; min advantage ≥+0.02 AND
  paired |bias| reduction ≥2 MC SE; robustness ≥70% of stress grid).
- Surveyed code: depth-floor implementable as three localized tree-move edits in
  `dev/rbart.R` (SD1 ensemble_init = complete depth-d₀ trees; SD2 reject prunes with
  NOG depth ≤ d₀; SD3 deeper moves unchanged). dgp5 + cluster-harness comparators
  reused (GLM, Vanilla-BART, GLMM-iCAR, GLMM-SPDE; slow spatial-ML at low reps).

## Iteration 1 — IMPLEMENT depth-floor + SMOKE TEST
- Implemented `fit_bartsimp_pg_sd` wrapping the d₀ knob; legacy = d₀=0.
- Smoke test (2 reps) → `dgp5_paired_smoke_per_rep.csv`. Valid output, no NaN, ranges
  plausible, SD output distinct from legacy, ~40 s/rep. PASS.

## Iteration 2 — PRIMARY paired A/B (N=500, A=20, 24 reps)
- Ran full roster + paired legacy (d₀=0) vs SD (d₀=3). Log: `primary_run.log`.
- Outputs: `dgp5_paired_primary_aggregated.csv`, `_per_rep.csv`, `_deltas.csv`.
- Result: admin cov80 0.762 → 0.777; paired Δcov80 = +0.0145 ± 0.0102 (≈1.4 MC SE,
  one-sided p≈0.08 — directional, not significant). cov95 0.944→0.942, ratio
  1.036→1.034, cluster cov80 ~1.000, pixel cov80 0.869→0.866 (all preserved).
  Δ|bias| = −0.00066 ± 0.00033 (2 MC SE); frac_reps_bias_reduced 0.625,
  frac_reps_cov80_up 0.417. JUDGE: right sign, modest magnitude, sub-2-MC-SE →
  proceed to stress tests for corroboration (not a kill; mechanism intact).

## Iteration 3 — SENSITIVITY: d₀ sweep (d₀ ∈ {0,1,2,3,4}, 16 reps)
- Log: `d0_sweep.log`. Output: `dgp5_d0sweep_aggregated.csv`, `_per_rep.csv`.
- cov80 by d₀ = {0:0.793, 1:0.815, 2:0.824, 3:0.815, 4:0.827}: broad plateau,
  +0.02–0.03 over legacy, clears the +0.02 bar. pixel cov80 flat ~0.86–0.87, cov95
  ~0.93–0.94 throughout. Cleanest positive evidence; knob robust, not knife-edge.
- Doubles as ABLATION: d₀=0 row = legacy (floor removed).

## Iteration 4 — SCALE + ADVERSARIAL
- SCALE A ∈ {5,20,40} paired (12 reps). Log: `scale.log`. Output:
  `dgp5_scale_aggregated.csv`, `_per_rep.csv`. SD ≥ legacy at every A
  (5:0.833→0.850, 20:0.782→0.816, 40:0.797→0.803); ratio in-band; pixel unchanged.
- ADVERSARIAL σ_gp ∈ {0.5,1,1.5} paired (12 reps). Log: `adversarial.log`. Output:
  `dgp5_adversarial_aggregated.csv`, `_per_rep.csv`. Floor recenters toward 0.80
  from both sides: σ=0.5 0.782→0.816 (up), σ=1.0 0.820→0.812, σ=1.5 0.821→0.804
  (gently down). Recentering/calibration, not monotone inflation.

## Iteration 5 — FINALIZE
- Verified every reported number against the CSVs (primary aggregated/deltas, d₀
  sweep, scale, adversarial) and the `primary_run.log` PAIRED DELTAS block.
- Confirmed figures generated (`make_figures.log`: fig1–fig3b present, ~70–107 KB).
- Verified `comparator_inventory.md` matches what ran (proposed + C1–C5 FULL in
  primary table; C6–C8 low-rep context; E1–E4 deferred) — gate COMPLETE, no edit
  needed.
- Verdict: **CONDITIONAL GO**. Wrote `validation_report.md`. All four stress types
  present and CSV-backed (SENSITIVITY=d₀ sweep, ABLATION=d₀=0 arm, SCALE,
  ADVERSARIAL); none INCOMPLETE.
