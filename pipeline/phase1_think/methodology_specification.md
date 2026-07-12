# Methodology Specification — BARTSIMP-PG with Semi-Dense Coarse-Layer Prior

> Phase 1 (methodology-architect) output. Target venue: JASA Applications &
> Case Studies / AOAS. This is a **framework paper** with a calibrated-joint-UQ
> headline, a bias-reduction extension to the existing BARTSIMP-PG sampler, and
> a strengthened theory arc. It is *not* a Gaussian-paper sequel.

## Proposed Method: BARTSIMP-PG-SD (semi-dense BARTSIMP-PG for calibrated SAE)

### One-Sentence Summary

A Bayesian model-based-geostatistics framework for small-area estimation that
couples a BART mean, a structured SPDE–Matérn Gaussian-process spatial field,
and Pólya–Gamma augmentation for the binomial-logit likelihood, and adds a
**semi-dense coarse-layer tree prior (a hard depth-floor d₀)** that recenters
the aggregation-scale (Admin-1) posterior to remove the BART non-Donsker mean
bias — delivering credible intervals that are calibrated *jointly* at pixel,
cluster, and area scales simultaneously.

### Core Idea

BARTSIMP-PG already produces honest-to-conservative coverage at the fine
(pixel/cluster) scale but **under-covers at the Admin-1 aggregation scale**
(cov80 ≈ 0.77 vs 0.80; cov95 ≈ 0.94). The diagnosis is decisive: because cov95
is near-nominal while cov80 lags, and the calibration ratio RMSE/post_sd ≈ 1.08
implies a standardized center bias δ = √(1.08²−1) ≈ 0.41 SD, the admin
under-shoot is a **center bias, not a width problem**. The bias is the
predicted *BART non-Donsker* effect: the flexible tree mean carries a slowly
vanishing low-frequency bias that the credible interval cannot see, and that
*survives aggregation* precisely because area averaging loads on the coarse
(low-resolution) component of the mean — the same component a pruning tree
prior under-represents. Castillo–Ročková (2021) prove that tree-prior
functional under-coverage is caused by the prior *pruning the coarse layers*,
and that the remedy is a **prior modification (keep the coarse layers dense)**,
not estimator-side debiasing. We therefore force the coarsest d₀ tree layers to
be always present (a deployable hard depth-floor, d₀ ≈ 2–3 = √(log N)), which
recenters the admin posterior toward nominal **without widening the fine-scale
bands** (the floor touches only the coarse partition; deep-tree
partition-instability and boundary-GP widening — the fine-scale conservatism —
are untouched). This is a *Bayesian-coherent* fix that, unlike the
RoBART/BLY25 debiasing correction, is a single prior knob with no new estimator
and is the empirical instantiation of the viable forward-route theory
(Theorem 2 target).

### Mathematical Specification

**Data.** Clusters i = 1,…,N. Cluster i has a covariate vector
x_i ∈ X = [0,1]^p (p ≈ 5 in simulation, up to ≈49 in the high-dim regime),
a spatial location s_i ∈ D ⊂ ℝ², trials n_i, and successes
k_i ∈ {0,…,n_i}. Outcome model:

  k_i | p_i ~ Binomial(n_i, p_i),  logit(p_i) = η_i = α₀ + f(x_i) + z(s_i) [+ ξ_i].

- α₀: scalar intercept, absorbed as the GLM-style fixed offset
  binary_offset = logit(p̄), p̄ = Σk_i / Σn_i (so trees model deviations).
- f: X → ℝ: BART ensemble mean, f(x) = Σ_{t=1}^{m} g(x; T_t, M_t), with m trees.
  Each tree (T_t, M_t) has topology T_t and leaf parameters β_{t,ℓ} | T_t ~
  N(0, τ²) iid, τ² = (3 / (k_BART √m))². The tree-topology prior is the Chipman
  branching prior P(split at depth d) = α_tree (1+d)^{−β_tree}
  (defaults α_tree = 0.95, β_tree = 2), **modified by the semi-dense
  depth-floor below**.
- z: D → ℝ: a mean-zero Gaussian process with a Matérn-ν covariance (ν = 1
  fixed), represented through the SPDE/GMRF approximation of Lindgren et al.
  (2011): z(s) = Σ_j A(s)_j u_j with u ~ N(0, Q_ψ^{−1}), Q_ψ =
  inla.spde.precision(spde; θ = (log ρ, log σ_m)), A the sparse
  mesh-projection matrix. Hyperparameters ψ = (ρ, σ_m) (spatial range,
  marginal SD) carry a PC prior.
- ξ_i (optional): cluster random effect ξ_i ~ N(0, σ²_ξ) for overdispersion;
  default off.

**Pólya–Gamma augmentation (Polson–Scott–Windle 2013).** Introduce
ω_i ~ PG(n_i, η_i). Then conditional on ω = (ω_i),

  k_i − n_i/2 = ω_i η_i + (noise),  equivalently the pseudo-response
  ȳ_i = (k_i − n_i/2)/ω_i  satisfies  ȳ_i | ω_i, η_i ~ N(η_i, 1/ω_i),

rendering the model **conditionally Gaussian** in (f, u, ξ) with
heteroscedastic precision Ω = diag(ω_i). Optional DEFF-tempering replaces
(n_i, k_i) by (w_i n_i, w_i(k_i − n_i/2)) with w_i = 1/[1+(n_i−1)ρ_DEFF] to
deflate effective trials by the design effect.

**Semi-dense coarse-layer prior (THE EXTENSION).** Define a depth-floor
d₀ = d₀(N) ≈ ⌈√(log N)⌉ ∈ {2,3} for deployed N ∈ [200, 10⁴]. Modify each
tree's topology prior so that:

  (SD1) every tree is **initialized complete** to depth d₀ (a full binary tree
        of 2^{d₀} leaves), not as a stump;
  (SD2) **pruning is forbidden at any node of depth ≤ d₀**: in the
        birth/death (prune) move, any proposal whose collapsed internal node
        sits at depth ≤ d₀ is rejected (probability-0 under the modified prior);
  (SD3) for depth > d₀ the Chipman branching prior and all
        grow/prune/change/swap moves are **unchanged** (sparsity and adaptivity
        operate only below the floor);
  (SD4) leaf parameters remain random, β_{t,ℓ} | T_t ~ N(0, τ²), at every depth.

Formally, writing f in an L²(X)-normalized multiresolution (Haar/CART) basis
{ψ_{ℓk}} with level ℓ and cell index k, the modification sets the inclusion
indicators γ_{ℓk} = 1 for all k at all levels ℓ ≤ d₀ (complete coarse layers),
matching the Castillo–Ročková "semi-dense" hypothesis (their j₀(n) ≍ √(log n)).

**Estimand (the deliverable).** For an Admin-1 region a (a *spatial* partition
cell), the area-average prevalence functional

  χ^a = Σ_{i ∈ a} w_i p_i,  w_i = n_i / Σ_{j∈a} n_j  (the in-sample / realized-
        cluster weighting used by `compute_admin_cov`),

with latent-scale companion χ^a_η = Σ_{i∈a} w_i η_i = α₀ + ḡ_a + z̄_a,
ḡ_a = Σ w_i f(x_i), z̄_a = Σ w_i z(s_i). The posterior of χ^a is obtained by
**joint** draws of (f, z) across all clusters in a (so the area-mean law carries
the cross-cluster covariance — the property no marginal-interval competitor
has).

**Posterior summary.** Central (1−γ) credible interval I^a_N(γ) from the
empirical quantiles of the post-burn χ^a draws. Pixel/cluster intervals are the
analogous per-point quantiles of η (or p) draws.

**Assumptions (stated).**
(A1) Conditional independence of k_i given η_i (binomial). (A2) Additive
latent structure η = α₀ + f(x) + z(s): no f×z interaction (an interaction
induces a Θ(1) contraction-bias floor — Mistake 14). (A3) Matérn ν fixed
(= 1); SPDE/GMRF approximation accurate at the chosen mesh. (A4) Informative
sampling design respected through n_i weighting; within-cluster overdispersion
(DEFF) is *not* fully modeled by default (acknowledged weak regime). (A5) For
the theory companion only: Gaussian likelihood and the dimension condition
s_a > p/2 (see Theoretical Properties).

### Algorithm

Gibbs/MH sampler (pure R; the existing `run_sampler_bartsimp_pg` with the
depth-floor flag added to `dev/rbart.R`). Per iteration:

```
Inputs: X (N×p), locs (N×2), n_i, k_i; m_trees, k_BART, alpha_tree, beta_tree,
        d0 (depth-floor), PC-prior + MH step sizes; mesh/SPDE objects.

Init:   binary_offset = logit( sum(k_i)/sum(n_i) )
        ensemble = m complete trees of depth d0   # (SD1) — replaces stump init
        u = 0, z = 0, sigma_m, rho = init; f = 0

For it = 1..n_iter:
  1. eta_i      = binary_offset + f_i + z_i [+ xi_i]
  2. omega_i    ~ PG(n_i, eta_i);   ybar_i = (k_i - n_i/2)/omega_i      # PG step
  3. Q_psi      = SPDE precision(sigma_m, rho)
  4. BART sweep on residual r_i = ybar_i - binary_offset - z_i [- xi_i],
       precision Omega:                                                  # mean
       for each tree t: birth/death/change/swap MH using the marginal
       leaf-integrated likelihood, with cutpoint grid;
       *** (SD2) reject any prune whose NOG node depth <= d0 ***
       *** (SD3) depth>d0 moves unchanged ***
       draw leaf params beta ~ N(., tau^2) (SD4); f_i = sum_t g_t(x_i)
  5. GP step (Woodbury): yhat_i = ybar_i - binary_offset - f_i [- xi_i];
       MH update of (sigma_m, rho) under PC prior; draw u | . (Gaussian),
       z = A u   (optionally restricted-spatial-regression projected)
  6. [optional] xi_i, sigma2_xi Gibbs (overdispersion);  DART var-probs
  7. store draws; per-iter grid/area predictions f_pred, z_pred
Post: aggregate area functional chi^a = sum_{i in a} w_i plogis(eta_i) over
      post-burn draws -> credible intervals at pixel/cluster/area scales.
```

The depth-floor changes only steps Init and 4 (three localized edits to
`rbart.R`: complete-tree init, prune guard, death-move guard). Everything else
— PG, Woodbury GP, hyperparameter MH, prediction — is byte-identical to the
current sampler when d₀ = 0.

### Parameters and Tuning

| Parameter | Role | Default | Valid range / How to set |
|-----------|------|---------|--------------------------|
| **d₀** (depth-floor) | semi-dense coarse-layer floor — the bias-reduction knob | ⌈√(log N)⌉ (≈2–3) | {0,1,2,3,4}; 0 = legacy BARTSIMP-PG. Set to √(log N); sensitivity-sweep {2,3} on dgp5. Larger d₀ overcompletes coarse partition (variance cost). |
| m_trees | ensemble size | 100 (prod), 20 (light) | 20–200; inert above ≈25 for tree size (surface ablation). |
| k_BART | leaf-prior amplitude; τ = 3/(k√m) | 2 | 1–5; amplitude/width lever, depth-decoupled. |
| α_tree | Chipman branch base prob | 0.95 | (0,1); mild depth modulator. |
| β_tree | Chipman branch decay (depth lever) | 2 | 1–5; structural depth/width lever (below the floor). |
| ν (Matérn) | spatial smoothness | 1 (fixed) | fixed; ν≤1/2 favored by Theorem 1 regime. |
| ρ₀, σ₀ (PC) | spatial range/SD PC-prior anchors | 0.05, 0.2 | anchored at two-stage diagnostic so prior doesn't fight data. |
| ρ_DEFF | DEFF-tempering intra-cluster corr | 0 | ≥0; >0 inflates post_sd ~√DEFF (DEFF robustness arm). |
| bd_mode | tree base distribution | "conditional" | {"conditional","marginal","soft"}; "soft" = SoftBART width lever. |
| soft_bw | SoftBART gate bandwidth | 0.10 | >0; →0 recovers hard CART. |
| enable_re | cluster overdispersion RE | FALSE | {T,F}; BART tends to absorb it (weak DEFF gain). |
| n_iter / burn | MCMC length | 5000 / 2000 (prod) | light 600/200 for sim grid. |
| chol_perm | fill-reducing Cholesky | FALSE | TRUE ≈28× faster, distributionally identical. |

### Theoretical Properties (expected)

Each property marked PROVEN / CONJECTURED / UNKNOWN. The honest theory state
(per `theory/HANDOFF.md`) is preserved; the BLY25/RoBART route is **dead** and
not invoked.

1. **(Theorem 1) Matérn-only (spatial-only) functional coverage at Admin-1.**
   For the GP component alone (linear/known mean), the area-mean posterior
   attains nominal-to-conservative coverage under the KvdVvZ (2011) functional-
   BvM machinery with spatial representer index q = 1/4, p = 0, prior index
   ν/2. — **PROVEN (modulo two fixes, essentially rigorous).** Known
   correction: the draft is *mis-regimed* (cites KvdVvZ Thm 5.3 conservative
   while q = 1/4 > p = 0 is the Thm 5.4 nominal-BvM regime; the β⋆ > ν/2 + 3/4
   vs + 1/4 arithmetic flag is the symptom). After the fix it is a *nominal*
   coverage result; three NEEDS-DETAIL items remain (boundary curvature,
   n→N_a rate, width constant C_max). Gaussian likelihood.

2. **(Theorem 2, single-tree Gaussian — the front-runner deliverable)
   Nominal forward BvM for the Admin-1 functional of a single semi-dense
   Bayesian-CART tree + Matérn GP, Gaussian likelihood.** The marginal
   posterior of χ^a satisfies, in P_{η⋆}-probability,
   β_R( L(√N(χ^a − χ̂^a) | Y), N(0, V_smp^a) ) → 0 with the **in-sample**
   efficient variance V_smp^a = σ²/π_a, so the plain credible interval is
   **calibrated for the in-sample/finite-population estimand with NO debiasing**.
   Assembled by conditioning on the tree (η|T is a single GP on the product
   domain ⇒ KvdVvZ Prop 3.2/Thm 5.4 apply conditionally) and marginalizing via
   Castillo–Ročková (2021) Thm 9 semi-dense machinery. — **CONJECTURED
   (audited "MAJOR REVISION, no fatal holes").** The new piece P2
   (representer decay) is **CONFIRMED**: the C–R summability Σ_ℓ ℓ c_ℓ < ∞ holds
   iff the within-region covariate density is Hölder-smoother than p/2, i.e.
   **s_a > p/2**. Remaining live math: E-ii (in-sample↔population representer
   reconciliation through C–R Step 9), E-i (tree⊕GP product-domain extension),
   E-iii (dyadic vs CGM axis-split). The semi-dense depth-floor (SD1–SD4) is the
   exact prior modification this theorem requires — so the **method is built to
   meet the theorem's hypothesis**.

3. **(Mechanism prediction) The semi-dense depth-floor recenters the Admin-1
   posterior toward nominal.** Castillo–Ročková: tree-prior functional
   under-coverage is caused by the prior pruning coarse layers; forcing them
   dense restores finite-dimensional-distribution Gaussian convergence and
   nominal coverage. Empirically: depth-floored BARTSIMP-PG admin-1 cov80 →
   0.80 relative to vanilla 0.77, without widening pixel bands. —
   **CONJECTURED.** C–R prove only semi-dense ⇒ BvM (no converse), and a
   competing **ensemble-CLT + coarse-scale shrinkage** story (τ² ∝ m^{−1/2}
   damps the coarse mean) also fits the cov80 ≈ 0.78 dip; the A/B experiment
   discriminates them.

4. **(Cluster floor, Theorem 9.X) Forced over-coverage at the cluster scale.**
   The average cluster (1−α) coverage converges to a floor strictly above
   nominal, forced by binomial discreteness at small n_i (analytic_floor_80 ≈
   0.92–0.93). — **PROVEN.** This is why the fine scale must not be "fixed"; the
   floor is a feature, and the depth-floor leaves it intact.

5. **(Dimension wall) The clean single-tree theorem requires s_a > p/2.**
   p = 49 ⇒ s_a > 24.5 — no single-tree sparsity rescue; deployability needs
   BART's low-effective-dimension/additive-sieve structure beyond a single
   tree. — **PROVEN-as-limitation (the bound is rigorous; the rescue is
   UNKNOWN).** Mirrors the binomial p<2 cubic-rate wall.

6. **(PG/binomial transfer) Conservative/nominal coverage for the deployed
   binomial-logit model.** — **CONJECTURED / partly UNKNOWN.** The Gaussian
   sibling (properties 1–2) is the anchor; the continuous→binomial transfer is
   the standing open gap. PG vs Albert–Chib vs Holmes–Held augmentation
   invariance of the *point estimate* is empirically PROVEN (surface study);
   only band width is augmentation-dependent (a k knob).

### Comparators for Evaluation

Core-8 (every results table; the spine):

1. **BARTSIMP-PG-SD (proposed)** — BART + SPDE-Matérn GP + PG + semi-dense
   depth-floor.
2. **GLM (logistic)** — linear, non-spatial floor; isolates the value of
   flexibility + spatial structure.
3. **Vanilla-BART** — proposed minus the GP (spatial ablation); isolates the
   GP's area-level work (admin cov80 collapses to ≈0.28 without it).
4. **GLMM-iCAR** — areal random-effect baseline (INLA iCAR); the disease-mapping
   default.
5. **GLMM-SPDE** — INLA workhorse (linear mean + Matérn); the LMIC standard and
   the strongest linear-mean spatial comparator.
6. **disaggregation** — native binomial + areal + joint TMB UQ; the fair native
   analog with a genuine joint law.
7. **RF-GLS** — spatial random forest, marginal per-cluster UQ; the
   constant-bias / laundered-width camp.
8. **SPARforest** — spatial-ML with independent per-cluster Gaussian draws; the
   "near-nominal-by-coincidence, scale-dependent" cautionary case.

Extended arms (matched per-rep designs; cited where they bite): **BARTSIMP-PG
(hard, d₀=0)** — the direct ablation that isolates the depth-floor;
**BARTSIMP-PG-soft** — SoftBART width lever (stacks with the floor);
**BARTSIMP-BB-RE** — Beta-Binomial RE (DEFF probe); **dbarts**, **BART::lbart**
— Albert–Chib(probit)/Holmes–Held(logit) augmentation cross-checks; **PPI**,
**MA-LASSO** — design-culture area-mean-direct (no pixel UQ). Diagnostic
overlay: the **calibration ratio RMSE/post_sd** ("earned vs laundered"), read
across aggregation scales, with the **constant-bias camp** (SPARforest, RF-GLS)
vs **averageable-misfit camp** (GLMM-SPDE, disaggregation) framing.

Deferred (not carried): spatial conformal `scp`, NN-GLS `geospaNN` — marginal /
Gaussian-only, no joint area-mean law.

### Expected Strengths

1. **Calibrated joint UQ at *every* aggregation scale — the discriminating
   property.** BARTSIMP-PG already holds ratio ≈ 1 at A ∈ {40,20,5} (1.01/1.08/
   0.96) with bias ≈ 0.001, where the constant-bias camp (SPARforest, RF-GLS)
   sweeps its ratio *up* through 1 and the averageable-misfit camp (GLMM-SPDE,
   disaggregation) drifts *down* to 1 — each calibrated at only one crossing
   point. The depth-floor pushes the one soft spot (admin cov80 0.77 → ~0.80)
   onto the line without disturbing the others. Mechanism: joint field draws +
   a recentred coarse mean.
2. **Bias-reduction is a single Bayesian-coherent prior knob, not a bolt-on
   estimator.** The depth-floor is a prior modification (≈0 extra compute, three
   localized `rbart.R` edits, fully bit-verifiable) that is *the same object*
   the forward-route theorem requires — so the empirical fix and the theory
   target are mechanistically identical (a rare alignment). Contrast RoBART,
   which needs a new estimator and whose theory route is holed for this
   estimand.
3. **Honest fine-scale conservatism is preserved by construction.** The floor
   touches only the coarse partition; deep-tree partition instability and
   boundary-GP widening (the pixel/cluster over-coverage) are untouched, so the
   fix cannot trade fine-scale calibration for area-scale calibration — the
   explicit design constraint.
4. **Mechanism-discriminating design.** The depth-floor A/B distinguishes the
   non-Donsker-pruning story from the ensemble-CLT-shrinkage story — a
   scientific result either way, and the kind of mechanistic clarity a top
   venue rewards.

### Expected Weaknesses

1. **Admin-level observation overdispersion (DEFF) remains a non-robust
   regime.** BART absorbs the cluster RE (seen on Nigeria and in the BB-RE arm),
   so the model-based posterior SD does not track within-cluster DEFF; at
   ρ = 0.2 BARTSIMP variants are the worst-covered of the roster (ratio ≈ 1.84).
   The depth-floor does *not* help here (it is a mean-recentering knob, not a
   variance term). Honest scope statement; the DEFF-tempered likelihood
   (Route 3) is a partial, not full, mitigation.
2. **Wider intervals and ≈100× vanilla-BART compute.** The price of honest
   joint coverage; SoftBART recovers ≈7% width but not the cost.
3. **The clean theorem is single-tree, Gaussian-likelihood, low-p.** The
   deployed object (ensemble, binomial, p≈49) is covered only by conjecture; the
   s_a > p/2 dimension wall and the continuous→binomial transfer are open. The
   paper must state this prominently rather than over-claim.
4. **Weak covariate signal → worse point accuracy** (a finding, not a bug); and
   the floor adds a small variance cost when d₀ overcompletes the coarse
   partition relative to the true coarse structure.

### Minimum Viable Evaluation

**Experiment.** Run `sim_dgp5_coverage` (known 5-covariate nonlinear,
non-additive truth + Matérn field; admin-1 area-average prevalence the scored
estimand) as a **paired A/B** of **BARTSIMP-PG-SD (d₀ = ⌈√(log N)⌉)** vs
**BARTSIMP-PG (d₀ = 0)**, matched seeds (truth drawn first → byte-identical
truth+data per rep), N = 500, A = 20, ≥40 reps, light MCMC config
(m_trees = 20, n_iter = 600). Carry the full core-8 roster on the same designs.
Add the aggregation-scale sweep A ∈ {40, 20, 5} (the curve, not one point) and
the `sim_mbg_surface` pixel-coverage guardrail.

**Metrics.** Primary: admin-1 cov80 and the calibration ratio RMSE/post_sd at
each A. Guardrail: pixel cov80 (must stay in the conservative band, not
collapse). Secondary: admin RMSE, bias, post_sd, cov95; paired within-rep
deltas (Δcov80, Δbias, Δpost_sd).

**Evidence that constitutes success.** The paired depth-floor arm moves admin-1
cov80 from ≈0.77 toward ≈0.80 (reduced |bias|, Δbias < 0 in a majority of reps)
**while** pixel cov80 stays ≥ 0.90 (fine-scale calibration intact) and the
calibration ratio stays ≈ 1 across A ∈ {40,20,5}. A null result (floor doesn't
move cov80) is itself informative — it adjudicates toward the ensemble-CLT-
shrinkage mechanism and redirects the lever to τ² ∝ m^{−1/2}.

**Run note.** Framework R only: `env -u MallocNanoZone
/Library/Frameworks/R.framework/Resources/bin/Rscript` (bare `Rscript` resolves
to miniforge and segfaults Rcpp/INLA/BART).

### Success Criterion

- **Primary metric.** Admin-1 cov80 closer to nominal (0.80) than vanilla
  BARTSIMP-PG, *without breaking pixel calibration* (pixel cov80 must remain in
  the conservative band, i.e. ≥ 0.90 and ≤ 1.0), measured on `sim_dgp5_coverage`
  (admin) and `sim_mbg_surface` (pixel) at matched seeds.
- **Minimum meaningful advantage.** Admin-1 cov80 improvement of ≥ +0.02
  (e.g. 0.767 → ≥ 0.787, recovering ≥ ~⅔ of the 0.03 gap to nominal) **and** a
  paired admin |bias| reduction significant at ≥ 2 Monte-Carlo SE, with the
  calibration ratio held in [0.95, 1.10] at every A ∈ {40,20,5}.
- **Comparison target.** The arm to beat is **BARTSIMP-PG (hard, d₀ = 0)** — the
  same model minus the floor, paired rep-for-rep (this isolates the
  contribution). The framework must also remain the only roster member with
  ratio ≈ 1 at *every* aggregation scale (vs SPARforest/RF-GLS constant-bias and
  GLMM-SPDE/disaggregation averageable-misfit camps).
- **Robustness requirement.** The advantage (admin cov80 up, pixel calibration
  intact, ratio ≈ 1 across A) must hold in ≥ 70% of scenarios across the stress
  grid where it is applicable: aggregation sweep A ∈ {5,20,40}, spatial
  amplitude σ ∈ {0.5,1,1.5}, weak-vs-strong covariate signal, and the
  production-vs-light MCMC config probe (to split structural-δ from config-δ).
  It is **not** required to hold under the Beta-Binomial DEFF stress
  (ρ ∈ {0.1,0.2}) — that is the acknowledged non-robust regime, reported
  honestly, not a failure of the depth-floor.

### Novelty Statement

1. **Has this combination been published?** No. BART+GP+PG MBG for SAE with
   coherent joint multi-scale UQ is itself unpublished (the framework's first
   claim). The **semi-dense coarse-layer prior as a coverage-recentering knob
   for the aggregation-scale functional** has not been deployed in a spatial
   SAE / MBG model: Castillo–Ročková (2021) introduce semi-dense priors for a
   single regression-tree functional BvM (no spatial GP, no survey aggregation,
   no binomial likelihood), and no MBG paper imports it.
2. **Specific novel element beyond domain transfer.** (a) The *diagnosis* that
   Admin-1 under-coverage is a non-Donsker center bias loading on the coarse
   tree component (the cov95-vs-cov80 + ratio decomposition is a new diagnostic
   read); (b) the *mechanistic identification* of the semi-dense depth-floor —
   rather than RoBART debiasing — as the Bayesian-coherent remedy, with an A/B
   that discriminates pruning-bias from ensemble-shrinkage; (c) the
   **forward-route theorem target** (single-tree Gaussian BARTSIMP admin-1 BvM,
   in-sample V_smp = σ²/π_a, with the confirmed s_a > p/2 representer
   computation) whose hypothesis is *exactly the prior the method deploys* —
   tying a new theoretical result to the deployed knob.
3. **Reply to "this is just BART+GP applied to DHS data."** The contribution is
   not the application but the *calibrated-joint-UQ-across-scales* property and
   the prior-side recentering mechanism that delivers it, backed by a new
   functional-BvM theorem. A reviewer who says "just BART+GP" must explain why
   no BART+GP method on the board (Vanilla-BART, dbarts) holds the area
   functional, and why off-the-shelf spatial-ML (RF-GLS, SPARforest) only hits
   nominal at a single crossing scale. The novelty is the **scale-invariant
   calibration** and its prior-theoretic basis, not the ingredients.
   *Honest note:* the binomial-logit conservative-coverage theory remains
   conjectural; the rigorous theorem is the Gaussian single-tree sibling. We
   claim the framework + mechanism + Gaussian theorem, and mark the binomial
   transfer open.
