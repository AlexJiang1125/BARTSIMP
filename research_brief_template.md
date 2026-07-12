# Research Brief — BARTSIMP-PG: a Bayesian MBG framework for small area estimation

> **Framing.** This is *not* a sequel to the Gaussian SSTE paper. The target
> contribution is a **general Bayesian model-based-geostatistics (MBG) framework
> for small area estimation (SAE)** under complex survey designs, with **calibrated
> joint uncertainty** as the headline property. BART supplies the flexible mean, a
> structured SPDE–Matérn GP supplies the spatial field, Pólya–Gamma augmentation
> handles the binary/binomial likelihood, and the framework delivers coherent UQ at
> pixel, cluster, and area (Admin 1/2) scales simultaneously.

---

## 1. Problem Statement

- **The challenge:** Produce **calibrated** small-area prevalence estimates (and
  full pixel-level surfaces) for a binary health outcome from complex-survey data,
  where the covariate–response relationship is nonlinear/non-additive and there is
  residual spatial structure. The deliverable is area-level (policy), but modeling
  is point/cluster-level — so the method must carry honest uncertainty *through
  aggregation* to any partition.
- **Why existing methods fall short (the gap each competitor exposes):**
  - **Design-based / SAE** (direct/Hájek, model-assisted GREG, PPI): correct where
    their iid-residual variance model holds, but give **no pixel/cluster UQ** and
    understate SE under design effects (overdispersion).
  - **Linear model-based** (GLM, GLMM-iCAR, GLMM-SPDE, disaggregation): the linear
    mean cannot represent curvature/interactions → **mean-bias under-coverage** at
    area level (calibration ratio 1.2–1.4).
  - **Flexible-but-non-spatial ML** (vanilla BART, dbarts): good pointwise, but
    **collapse at the area functional** (admin cov80 ≈ 0.28) — no spatial term to
    aggregate.
  - **Spatial-ML** (RF-GLS, SPARforest): hit nominal coverage only by **laundered
    width** (interval wider than their error) or a constant bias that makes
    coverage **scale-dependent**, and emit only *marginal* per-cluster intervals —
    no coherent **joint** area-mean law.
  - The framework is the **only method calibrated across aggregation scales** with
    genuine joint field draws — that is the contribution.
- **Abstract restatement:** Bayesian estimation of a bounded linear functional
  (area average) of a latent surface `η = f(x) + z(s)` — an unknown nonlinear mean
  plus a GP — under a non-conjugate (binomial-logit) likelihood and an informative
  sampling design, with the requirement that the *posterior* recenter and rewidth
  to the efficient variance of that functional (a forward functional-BvM problem).

## 2. Data Assets

- **Real application:** 2014 **Kenya DHS** child wasting (under-5, WHZ < −2);
  Nigeria DHS is a secondary deployment. Stratified two-stage cluster sampling:
  ~1,584 EAs across 92 strata, 25 households/EA → 20,977 child measurements;
  production fit uses ~1,380 clusters.
- **Outcome:** binary wasting indicator, aggregated per cluster as `k_i` successes
  of `n_i` trials (binomial). (Gaussian sibling models WHZ continuously.)
- **Inputs:** per-cluster covariate matrix `X` (≈6 DHS geospatial rasters: pop
  density, night-light, vegetation index, temperature, precipitation, access to
  city; up to `p≈49` in the high-dim theory regime), cluster coordinates `locs`,
  `n_i`, `k_i`.
- **Simulation DGPs (the core empirical engine):**
  - **`sim_dgp5_coverage`** — known **5-covariate nonlinear, non-additive** truth +
    Matérn field; admin-1 area-average prevalence is the scored estimand.
  - **`sim_mbg_surface`** — pixel-level latent **probability surface** `p*(s)` on a
    G²=1600 grid (no binomial floor), for clean surface coverage / interpolation /
    extrapolation diagnostics.
  - Plus `sim_07` (thm verification), `sim_design_coverage`, `sim_sparsity`.
- **Known data quirks:** complex design (weights, stratification, clustering);
  within-cluster overdispersion / design effect (Beta-Binomial DEFF) that
  model-based SD does **not** track; urban–rural omitted as covariate; sparse
  Admin areas; **BART can absorb the cluster RE** (seen on real Nigeria data and in
  the BB-RE arm), collapsing the overdispersion term.

## 3. Domain Context

- **Terminology:** MBG (model-based geostatistics); SAE; WHZ / wasting; EA/cluster;
  Admin 0/1/2; SPDE/Matérn/GRF; direct vs. model-based estimates; calibration ratio
  `RMSE/post_sd`; "laundered" vs "earned" coverage; binomial information floor; DEFF.
- **Constraints:** policy deliverable is area-level, modeling is point-level → must
  aggregate with honest joint UQ; the informative design must be respected.
- **Field conventions:** maps of posterior medians + CI widths; comparison vs.
  design-based direct estimates; coverage reported with accuracy; INLA/SPDE is the
  LMIC default.
- **Practitioner expectations:** honest (near-nominal) intervals across scales even
  at some point-accuracy cost; reproducible maps; feasible on standard HPC.

## 4. Target Venue

- **Primary targets:** **JASA (Applications & Case Studies)** and **Annals of
  Applied Statistics (AOAS)** — methods-with-substantive-application venues that
  reward a principled framework + serious data application + (ideally) supporting
  theory. Frame as a methods/framework paper, not an extension note.
- **Document class:** currently elsarticle (`theory/paper.tex`, `NOTES.md`); would
  move to the target journal's class. natbib author-year citations.
- **Quality bar:** high on **all three axes** — (a) methodological novelty (joint
  BART+GP+PG MBG for SAE), (b) empirical rigor (10-method benchmark, multi-scale
  calibration, stress tests), (c) **theory** (BvM/coverage results; see §11).
  Reproducibility expected (SLURM arrays, seeded, byte-identical anchors).

## 5. Success Criteria

- **Primary:** admin-level (and pixel-level) **coverage near nominal** (cov80≈0.80,
  cov95≈0.95) with **calibration ratio `RMSE/post_sd ≈ 1` across aggregation scales**
  (A ∈ {5,20,40}) — the property no competitor holds at every scale.
- **Secondary:** admin/pixel RMSE, bias, post_sd, CRPS, Brier, log-loss; interval
  width/AIL/AIS; "earned vs laundered" diagnostic (ratio read alongside absolute
  RMSE + bias + a genuine joint-draw law).
- **Robustness / stress:** Beta-Binomial DEFF (ρ), spatial-amplitude crank,
  aggregation-scale sweep, weak-vs-strong covariate signal, mesh/hyperprior
  sensitivity, GP-limit checks, soft-vs-hard trees, PG vs Albert–Chib vs Holmes–Held
  augmentation invariance.
- **Compute budget:** production fit (`m_trees=100, n_iter=5000, burn=2000,
  n_cl≈1380, n_mesh≈730`) ~10–15 min with `chol_perm=TRUE`; sims run as cluster
  SLURM arrays under a 2 h/task wall (light config `m_trees=20, n_iter=600`).
- **Interpretability:** partial-dependence plots, spatial field maps, area summaries.

## 6. Comparator Methods

> The pipeline enforces that **every listed method appears in every results table**.
> The **core roster (8)** is the spine of every scorecard; the **extended arms** are
> run on the *same per-rep designs* (matched seeds) and appear where relevant.
> Wrappers in `dev/competing/` (`fit_<method>(...)`) + `experiments/.../R/`.

**Core roster (every table):**

| # | Method | Role / what it isolates | Implementation |
|---|--------|-------------------------|----------------|
| 1 | **BARTSIMP-PG** (proposed) | BART mean + SPDE-Matérn GP + PG | `dev/bartsimp_pg.R` + `dev/rbart.R` |
| 2 | GLM (logistic) | linear, non-spatial floor | `glm` |
| 3 | Vanilla-BART | proposed **minus** GP (spatial ablation) | `BART`/PG, no field |
| 4 | GLMM-iCAR | areal random-effect baseline | INLA iCAR |
| 5 | GLMM-SPDE | INLA workhorse (linear + Matérn) | INLA SPDE binomial |
| 6 | disaggregation | native binomial+areal+joint TMB UQ | `disaggregation` |
| 7 | RF-GLS | spatial random forest (marginal UQ) | `RandomForestsGLS` |
| 8 | SPARforest | spatial-ML, per-cluster Gaussian draws | `R/spatial_ml_competitors.R` |

**Extended arms (run on matched designs; cited where they bite):**
- **BARTSIMP-PG-soft** (Linero–Yang soft trees) — Pareto-frontier sharpening of #1.
- **BARTSIMP-BB-RE** — Beta-Binomial random-effect variant (DEFF probe).
- **dbarts**, **BART::lbart** — Albert–Chib(probit) / Holmes–Held(logit)
  augmentation cross-checks (surface study).
- **Design-culture arm:** **PPI** (prediction-powered inference), **MA-LASSO**
  (model-assisted GREG) — area-mean-direct, no pixel UQ.
- *Deferred (no joint law):* spatial conformal `scp`, NN-GLS `geospaNN`.

## 7. Adjacent Fields

- **Semiparametric / debiased ML:** RoBART (Breunig–Liu–Yu 2025), one-step /
  augmented-IPW efficient variance, Riesz representers — source of the bias-
  correction idea and the missing-data-vs-deterministic-representer subtlety.
- **Bayesian nonparametric UQ:** Castillo–Ročková 2021 (UQ for Bayesian CART),
  KvdVvZ 2011 (GP functional BvM), Yang–Bhattacharya–Pati 2017 — the BvM toolkit.
- **Survey statistics / SAE:** model-assisted estimation, PPI, finite-pop vs
  superpopulation variance (Little–Rubin), design effects.
- **Disease mapping / spatial epi:** binomial CAR/BYM, geostatistical prevalence.
- **Latent-variable augmentation:** PG beyond logistic (negbin, multinomial),
  Albert–Chib probit, Holmes–Held logit — conjugacy-restoring tricks.

## 8. Evaluation Design

- **Replicates / limits:** ~40–50 paired reps per scenario (matched seeds: truth
  drawn first → byte-identical truth+data across arms), cluster SLURM 1–40 arrays,
  2 h/task wall; light MCMC for the grid, a production-config probe to split
  structural-δ from config-δ.
- **Evaluation types:** admin-1 area-mean coverage (`sim_dgp5_coverage`), pixel
  surface coverage (`sim_mbg_surface`), theorem-verification (`sim_07`), real-data
  Kenya/Nigeria application + out-of-sample CV (`app_oos_cv`), pixel maps.
- **Seed convention:** seed = base + rep; sampler bit-reproducible at
  `chol_perm=FALSE`; runs echo config; **Framework R only** + `env -u
  MallocNanoZone` (miniforge `Rscript` segfaults Rcpp/INLA/BART — recurring hazard).
- **Stress scenarios:** aggregation-scale sweep A∈{5,20,40}; Beta-Binomial DEFF
  ρ∈{0,0.1,0.2}; spatial-amplitude σ∈{0.5,1,1.5}; tree-hyperparameter ablation
  (m_trees,k,α,β); 1D hole+tails interpolation/extrapolation sandbox; augmentation
  invariance (PG/AC/HH); frequentist width-calibration (`post_sd/freq_sd` vs N,n_i).

## 9. Existing Gaps (the lingering issues to tackle)

- **The coverage asymmetry across scales (the central gap).** The framework
  **over-covers at the pixel/cluster scale** (cluster: binomial discreteness floor,
  forced cov80≈0.99–1.0, Thm 9.X; pixel surface: deep-tree partition instability +
  boundary GP widening, cov80≈0.98) **but under-covers at the aggregation
  (admin-1) scale** (cov80≈0.77, just below nominal). This is the asymmetry to fix:
  honest-to-conservative at fine scales, but the area-mean deliverable comes in
  *below* nominal. The admin under-shoot is a **center bias**, not a width problem
  (cov95 near-nominal at ≈0.94; calibration ratio RMSE/post_sd ≈ 1.08 ⇒ standardized
  bias δ≈0.4 SD) — the predicted **BART non-Donsker** slowly-vanishing mean bias that
  survives aggregation while the over-coverage cushion washes out. **So the task is
  to reduce that aggregation-scale bias** (recenter the admin posterior) without
  destroying the calibrated fine-scale behavior.
- **Bias-reduction candidates (to evaluate):**
  1. **PG-RoBART center-bias correction** (BLY25 Algorithm 1, binomial-logit port) —
     targets exactly the non-Donsker mean bias; the Gaussian sibling moves CP
     0.607→0.93. *Caveat:* the BLY route is **holed as theory** for the
     deterministic admin representer (Mistake 13), but remains a usable empirical
     correction.
  2. **Castillo–Ročková depth-floor prior knob** — a *prior modification* (keep
     coarse tree layers dense; hard depth-floor d₀≈2–3) rather than debiasing;
     testable in `sim_dgp5_coverage`, discriminates non-Donsker-bias vs
     ensemble-shrinkage mechanisms.
  3. **Soft trees** — orthogonal lever: sharpens width (−7% post_sd, Pareto) but
     does **not** move center bias; stacks with (1)/(2).
- **DEFF robustness** is an acknowledged weak regime (BB-RE buys no admin-level DEFF
  robustness — BART absorbs the RE).

## 10. Scope & Non-Goals

- **Do NOT** frame as a Gaussian-paper sequel — it is a standalone MBG-for-SAE
  framework.
- **Do NOT** claim universal robustness — be explicit that admin-level
  observation-overdispersion (DEFF) is a regime where the method is *not* the most
  robust; its strength is principled **joint** UQ under its modeling assumptions.
- **Do NOT** re-attempt the **BLY25/RoBART theory route** for the admin-1 BvM — it
  is fatally holed at O1 (confirmed by two independent audits); pursue the forward
  (Castillo–Ročková ⊕ KvdVvZ) route instead. (RoBART as an *empirical* correction is
  still fine.)
- **Do NOT** rewrite the sampler in compiled code for speed (deliberate pure-R
  design; perf work tracked separately, must stay bit-verifiable).
- **Limitations to acknowledge, not fix:** weak covariate signal → worse point
  accuracy (a finding); wider intervals than competitors (price of honest joint
  coverage); ~100× vanilla-BART compute; BART's flat extrapolation bias in the
  far field (the GP's job, not a tree knob); cubic-rate theory condition fails at
  deployed dimension (`p≈49`).

## 11. Theory Arc

> The user wants the theory strengthened beyond the single existing result. Docs
> live in `…/BARTSIMP_glm/theory/` — read **`HANDOFF.md`** first.

- **Theorem 1 (Matérn-only conservative/nominal coverage)** — essentially rigorous
  (`draft_thm_continuous.tex`). **Known bug:** mis-regimed (cites KvdVvZ Thm 5.3
  conservative regime while sitting at the Thm 5.4 nominal-BvM regime; the
  `β⋆ > ν/2 + 3/4` vs `+1/4` arithmetic flag is the symptom). Three NEEDS-DETAIL
  items (boundary curvature, n→N_a rate, width constant C_max).
- **Theorem 2 (combined BART+Matérn admin-1 BvM)** — the open prize.
  - **v3 / BLY25-RoBART route: DEAD** (O1 fatal, Mistake 13) — do not revive.
  - **Forward route (viable, not yet a theorem):** KvdVvZ Thm 5.4 (GP, exact
    conjugacy) ⊕ Castillo–Ročková 2021 Thm 4/9 (forward tree-prior functional BvM,
    semi-dense prior), combined by **conditioning on trees** + semi-dense
    marginalization. Open pieces P1 (tree⊕GP product domain), P2 (representer decay,
    `s_a > p/2` — a **dimension wall** at p=49), P3 (BARTSIMP violates semi-dense:
    ensemble out of scope, soft prior permits stumps, data-adaptive ≠ dyadic splits).
  - **Front-runner for a clean provable result: single-tree Gaussian BARTSIMP**
    (`draft_thm2_single_tree.md`) — one Bayesian-CART tree + Gaussian likelihood
    (keeps BART, drops binomial). Audited "MAJOR REVISION, no fatal holes"; P2
    confirmed; key correction = posterior calibrates the **in-sample** efficient
    variance `V_smp = σ²/π_a` (not the superpopulation `V_0`).
- **Goal for this thread:** land at least one rigorous *combined* (or clean
  single-tree) BvM that explains the multi-scale coverage behavior, and connect it
  to the bias-reduction levers in §9. PG/binomial conservative coverage remains
  conjectural; the Gaussian sibling is the anchor.
- **Discipline:** the `nonparametric-bayes-proof-discipline` skill (14 catalogued
  mistakes) auto-triggers — use its 11-item self-audit on any new proof.

## 12. Reference Papers

- **In-repo theory dir:** `…/BARTSIMP_glm/theory/` — `HANDOFF.md`,
  `draft_thm2_single_tree.md`, `forward_bvm_exploration.md`, `draft_thm_continuous.tex`,
  `draft_thm2_combined.tex` (holed), `limit_theory_bartsimp_pg.tex`, audit reports.
- **Empirical dirs:** `…/bartsimp_cluster/experiments/{sim_dgp5_coverage,
  sim_mbg_surface,sim_07,app_oos_cv}/` (each has `README.md` + `DESIGN.md`).
- **Manuscripts:** `theory/paper.tex` (binary draft), `theory/NOTES.md`,
  `bartsimp_for_ss/main.tex` (Gaussian parent), `main.bib`.
- **Must-read papers:**
  - Chipman, George & McCulloch (2010) — BART.
  - Polson, Scott & Windle (2013) — Pólya–Gamma augmentation.
  - Rue et al. (2009) / Lindgren et al. (2011) — INLA / SPDE-Matérn.
  - Gómez-Rubio & Rue (2018) — INLA-within-MCMC.
  - **Castillo–Ročková (2021, arXiv:1910.07635)** — UQ for Bayesian CART (forward
    route linchpin).
  - **Knapik–van der Vaart–van Zanten (2011)** — GP functional BvM (Thm 5.3/5.4).
  - **Breunig–Liu–Yu (2025, arXiv:2509.24634)** — RoBART (empirical correction; holed
    as theory here).
  - Jeong–Ročková (2023, arXiv:2008.06620) — adaptive minimax additive BART rates.
  - Petrillo (2024, arXiv:2410.20289) — BART→GP limit kernel.
  - Rao & Molina (2015) — SAE reference.
