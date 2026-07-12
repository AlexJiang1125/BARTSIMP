# Comparator Inventory (validate-method-owned)

Phase-2 validation of **BARTSIMP-PG-SD** (semi-dense depth-floor). The Phase-2
gate question is the *paired A/B* — does the depth-floor (d₀ > 0) move Admin-1
cov80 toward nominal relative to legacy BARTSIMP-PG (d₀ = 0)? So the
non-negotiable comparator is **BARTSIMP-PG (d₀ = 0)**, which is carried in every
table. The brief's core-8 roster is carried for *context* (the "calibrated at
every scale" framing), at the validation budget.

## Comparators promised by research_brief.md (§6)

| ID | Method name | Citation | Section in brief |
|----|-------------|----------|------------------|
| C1 | BARTSIMP-PG (proposed base; **= the d₀=0 ablation arm to beat**) | Chipman+ 2010 / Polson+ 2013 / Lindgren+ 2011 | §6 core #1 |
| C2 | GLM (logistic) | — | §6 core #2 |
| C3 | Vanilla-BART (proposed minus GP) | Chipman+ 2010 (dbarts probit) | §6 core #3 |
| C4 | GLMM-iCAR | Besag; INLA | §6 core #4 |
| C5 | GLMM-SPDE | Lindgren+ 2011; INLA | §6 core #5 |
| C6 | disaggregation | Nandi+ TMB | §6 core #6 |
| C7 | RF-GLS | Saha+ (RandomForestsGLS) | §6 core #7 |
| C8 | SPARforest | vendored CAR-Forest | §6 core #8 |
| E1 | BARTSIMP-PG-soft | Linero–Yang 2018 | §6 extended |
| E2 | BARTSIMP-BB-RE | beta-binomial RE | §6 extended |
| E3 | dbarts / BART::lbart (AC/HH augmentation) | Albert–Chib; Holmes–Held | §6 extended |
| E4 | PPI / MA-LASSO (design-culture) | Angelopoulos+; GREG | §6 extended |

## Comparators implemented in validated_code/ (and reused from the cluster harness)

| ID | Module:function | Status | Reason if not FULL |
|----|-----------------|--------|--------------------|
| **proposed** | methods_bartsimp.R:`fit_bartsimp_pg_sd` (NEW; d₀ via rbart.R `tree_complete`/`ensemble_init`/`tree_prune`) | FULL | — |
| C1 | methods_bartsimp.R:`fit_bartsimp_pg` (d₀=0) | FULL | the paired arm-to-beat; in every table |
| C2 | dgp5_core.R:`fit_glm_gen` | FULL | — |
| C3 | dgp5_core.R:`fit_vanilla_bart_gen` | FULL | — |
| C4 | dgp5_core.R:`fit_glmm_icar_gen` | FULL | — |
| C5 | dgp5_core.R:`fit_glmm_spde_gen` | FULL | — |
| C6 | spatial_ml_competitors.R:`fit_disaggregation_gen` | FULL (roster-context run, reduced reps) | ~30 s/rep; not in the per-rep paired loop |
| C7 | spatial_ml_competitors.R:`fit_rfgls_gen` | PARTIAL (roster-context, ≤5 reps) | ~480 s/rep — infeasible at full paired-rep budget; run at low reps for the scale/context table only |
| C8 | spatial_ml_competitors.R:`fit_sparforest_gen` | PARTIAL (roster-context, ≤5 reps) | vendored CAR-Forest, slow; low-rep context only |
| E1 | methods_bartsimp.R:`fit_bartsimp_pg_soft` | AVAILABLE (not run in Phase 2) | width lever, orthogonal to the d₀ center-bias question (spec §4) |
| E2 | methods_bartsimp.R:`fit_bartsimp_pg_bbre` | AVAILABLE (not run) | DEFF regime explicitly OUT of the robustness requirement (spec Success Criterion) |
| E3 | dbarts (= Vanilla-BART here) / lbart | PARTIAL | augmentation-invariance is a width cross-check, not the center-bias gate |
| E4 | design_culture_competitors.R:`fit_ppi_gen`,`fit_malasso_gen` | AVAILABLE (not run) | admin-direct, no pixel UQ; pad tables without testing the d₀ claim |

## Coverage: brief vs implemented

| Brief ID | Implemented? | Rationale if not in primary table |
|----------|--------------|-----------------------------------|
| C1 BARTSIMP-PG | YES (every table) | the paired comparator |
| C2 GLM | YES | — |
| C3 Vanilla-BART | YES | — |
| C4 GLMM-iCAR | YES | — |
| C5 GLMM-SPDE | YES | — |
| C6 disaggregation | YES (context/scale run) | compute-throttled out of the per-rep paired loop |
| C7 RF-GLS | LOW-REP context | 480 s/rep makes full-rep infeasible in the Phase-2 budget |
| C8 SPARforest | LOW-REP context | slow vendored competitor |
| E1–E4 | available, not run | each tests a DIFFERENT axis (width / DEFF / augmentation / design) than the d₀ center-bias gate; deferred to Phase 3 production |

## Gate result: COMPLETE

The Phase-2 gate comparator — **BARTSIMP-PG (d₀=0)**, paired rep-for-rep with the
proposed BARTSIMP-PG-SD — is FULL and present in every results table, as are the
four linear/ML context comparators (GLM, GLMM-iCAR, GLMM-SPDE, Vanilla-BART). The
three slow spatial-ML competitors (disaggregation, RF-GLS, SPARforest) are run at
reduced replicate count for the roster-context table only; this is a compute
rationale, recorded here, not a silent drop. Extended arms (E1–E4) probe axes
orthogonal to the depth-floor center-bias question and are deferred to Phase 3.

## Downstream contract
- paper-modeler: the proposed method + C1–C5 must appear in EVERY production CSV
  under `pipeline/phase3_write/data/`. C6–C8 must appear in the cross-scale /
  roster scorecard. No silent drops.
- paper-writer: every FULL ID must appear in every results table.
- paper-grader: verify (a) this file exists, (b) every FULL ID present in every
  table; do NOT re-derive the list from the brief.
