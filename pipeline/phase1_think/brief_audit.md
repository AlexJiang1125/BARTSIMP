# Brief Audit — BARTSIMP-PG MBG-for-SAE framework

Phase 0.5 contract. Every idea, example, and argument in `research_brief.md` is
inventoried below with a DEVELOP / COMPLEMENT / ACKNOWLEDGE / REJECT label.
Every DEVELOP and COMPLEMENT item must appear in `methodology_specification.md`;
every REJECT carries a substantive (non-"out of scope") reason.

## Core framing & problem (brief §1, §3)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-1 | General Bayesian MBG framework for SAE under complex survey designs (not a Gaussian-paper sequel) | DEVELOP | The paper's identity; spec is framed as a framework. |
| IDEA-2 | Calibrated **joint** uncertainty across pixel/cluster/area scales as the headline property | DEVELOP | Primary contribution and success metric. |
| IDEA-3 | BART mean + SPDE–Matérn GP + Pólya–Gamma binomial likelihood as the model class | DEVELOP | The base model class; the extension is built on it. |
| IDEA-4 | Estimand = bounded linear functional (area average) of latent surface η = f(x)+z(s) under non-conjugate likelihood + informative design ("forward functional-BvM problem") | DEVELOP | Defines estimand and the theory target. |
| IDEA-5 | Carry honest UQ *through aggregation* to any partition | DEVELOP | Drives the aggregation-scale-sweep evaluation. |

## The central gap and bias-reduction levers (brief §1, §9)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-6 | **Coverage asymmetry across scales** — over-cover pixel/cluster, under-cover admin-1 | DEVELOP | The problem the extension solves. |
| IDEA-7 | Admin under-shoot is a **center bias** (cov95≈0.94 near-nominal, ratio 1.08 ⇒ δ≈0.4 SD), not a width problem | DEVELOP | Diagnosis that selects the mechanism (recenter, not rewidth). |
| IDEA-8 | **Castillo–Ročková depth-floor prior knob** (keep coarse layers dense; d₀≈2–3) | DEVELOP | **Selected primary mechanism** — prior-side, Bayesian-coherent, theory-linked. |
| IDEA-9 | **PG-RoBART center-bias correction** (BLY25 Algorithm 1, binomial-logit port) | COMPLEMENT | Carried as a secondary/empirical comparator arm and ablation, NOT the headline (its theory route is holed for the deterministic admin representer — Mistake 13). |
| IDEA-10 | **Soft / SoftBART trees** (Linero–Yang) | COMPLEMENT | Orthogonal width/efficiency lever; stacks with the depth-floor but does not move center bias (verified: dgp5 bias +0.0014→+0.0013). |
| IDEA-11 | Discriminate non-Donsker-bias vs ensemble-shrinkage mechanisms via the depth-floor A/B | DEVELOP | Built into the evaluation as a mechanism-discrimination experiment. |
| IDEA-12 | BART non-Donsker slowly-vanishing mean bias survives aggregation while fine-scale over-coverage washes out | DEVELOP | The mechanistic story the spec rests on. |

## Theory arc (brief §11)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-13 | Theorem 1 (Matérn-only conservative/nominal coverage), essentially rigorous | DEVELOP | Stated as PROVEN-modulo-fixes; mis-regime fix noted. |
| IDEA-14 | Theorem 1 mis-regime bug (cites KvdVvZ 5.3 conservative while at 5.4 nominal; β⋆>ν/2+3/4 vs +1/4) | COMPLEMENT | Logged as a known correction in Theoretical Properties. |
| IDEA-15 | **Theorem 2 (combined BART+Matérn admin-1 BvM)** as the open prize | DEVELOP | The theory target. |
| IDEA-16 | Forward route: KvdVvZ Thm 5.4 ⊕ Castillo–Ročková 2021 Thm 4/9 via condition-on-trees + semi-dense marginalization | DEVELOP | The viable route, marked CONJECTURED with open pieces P1/P2/P3. |
| IDEA-17 | **Single-tree Gaussian BARTSIMP** clean provable result (front-runner; audited MAJOR REVISION, no fatal holes) | DEVELOP | The deliverable theorem; in-sample V_smp=σ²/π_a calibration. |
| IDEA-18 | BLY25/RoBART theory route for admin-1 BvM | REJECT | **Fatally holed at O1 (Mistake 13), confirmed by two independent audits**: BLY's debiasing relies on a pointwise conditional-zero-mean that the random missing-data indicator R has but the deterministic admin representer 1_a/π_a lacks; leftover Θ(ε_N) bias ⇒ √N(b−b̂)→∞ in the non-Donsker regime. Do not revive as theory (kept only as an empirical correction, IDEA-9). |
| IDEA-19 | In-sample (V_smp=σ²/π_a) vs superpopulation (V_0) estimand distinction | DEVELOP | Sharpens what the posterior calibrates; one-line real-data caveat. |
| IDEA-20 | Dimension wall s_a>p/2 (p=49 ⇒ s_a>24.5); cubic-rate condition fails at deployed dim | ACKNOWLEDGE | Stated as an honest scope limit of the clean theorem; not fixed (needs a sparsity/additive-sieve argument beyond a single tree). |
| IDEA-21 | PG/binomial conservative coverage remains conjectural; Gaussian sibling is the anchor | ACKNOWLEDGE | Honest theory-status statement. |
| IDEA-22 | nonparametric-bayes-proof-discipline 11-item self-audit | ACKNOWLEDGE | Process guard for any new proof; not a paper artifact. |

## Comparators (brief §6)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-23 | Core-8 roster: BARTSIMP-PG, GLM, Vanilla-BART, GLMM-iCAR, GLMM-SPDE, disaggregation, RF-GLS, SPARforest | DEVELOP | Carried verbatim as the evaluation spine. |
| IDEA-24 | Extended arms: BARTSIMP-PG-soft, BARTSIMP-BB-RE, dbarts, BART::lbart, PPI, MA-LASSO | COMPLEMENT | Run on matched designs; cited where they bite. |
| IDEA-25 | Deferred: spatial conformal `scp`, NN-GLS `geospaNN` | REJECT | Emit only *marginal* per-cluster intervals / Gaussian-only ⇒ no honest joint area-mean law, which is the contribution's discriminating axis; including them adds rows without testing the joint-UQ claim. |
| IDEA-26 | Calibration-ratio "earned vs laundered" diagnostic | DEVELOP | The adjudication instrument; primary secondary-metric. |
| IDEA-27 | Constant-bias camp (SPARforest, RF-GLS) vs averageable-misfit camp (GLMM-SPDE, disaggregation) framing | DEVELOP | The cross-aggregation-scale narrative that proves "calibrated at every scale". |

## Data & evaluation (brief §2, §5, §8)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-28 | Kenya DHS 2014 child wasting application (Nigeria secondary) | COMPLEMENT | Real-data deployment; the empirical evaluation is sim-led. |
| IDEA-29 | `sim_dgp5_coverage` (5-cov nonlinear non-additive truth + Matérn) admin-1 coverage | DEVELOP | The Minimum Viable Evaluation harness. |
| IDEA-30 | `sim_mbg_surface` pixel surface coverage (G²=1600, no binomial floor) | DEVELOP | Fine-scale-calibration guardrail (must not break). |
| IDEA-31 | Calibration-ratio RMSE/post_sd ≈ 1 across A∈{5,20,40} as primary success metric | DEVELOP | Primary success criterion. |
| IDEA-32 | Stress suite: BB-DEFF ρ, spatial-amplitude, aggregation sweep, weak/strong signal, mesh/hyperprior, GP-limit, soft-vs-hard, PG/AC/HH augmentation invariance | COMPLEMENT | Robustness battery; the extension must survive it. |
| IDEA-33 | Matched-seed paired reps (truth drawn first → byte-identical across arms) | DEVELOP | Evaluation design constraint baked into the MVE. |
| IDEA-34 | Compute budget (production ~10–15 min; light config for sims) | ACKNOWLEDGE | Feasibility constraint; depth-floor adds ~0 cost (prior knob). |

## Scope, non-goals, limitations (brief §10)

| ID | Item | Label | Note |
|----|------|-------|------|
| IDEA-35 | Do NOT claim universal robustness; admin-level DEFF is a non-robust regime (BART absorbs the RE) | ACKNOWLEDGE | Stated as an Expected Weakness. |
| IDEA-36 | Do NOT rewrite the sampler in compiled code (pure-R, bit-verifiable) | ACKNOWLEDGE | Implementation constraint honored (depth-floor is a small pure-R prior edit). |
| IDEA-37 | Wider intervals than competitors (price of honest joint coverage); ~100× vanilla-BART compute | ACKNOWLEDGE | Honest Expected Weakness. |
| IDEA-38 | BART flat extrapolation bias in the far field is the GP's job, not a tree knob | ACKNOWLEDGE | Bounds what the depth-floor can do; far-field stays GP-handled. |
| IDEA-39 | BART can absorb the cluster RE / overdispersion (seen on Nigeria + BB-RE arm) | COMPLEMENT | Explains the DEFF weakness; the DEFF-tempered likelihood (Route 3, already in code) is the partial mitigation noted. |

## Reject summary (substantive reasons)

- **IDEA-18 (BLY25/RoBART theory route for admin-1 BvM)** — fatally holed at O1
  (Mistake 13), confirmed twice; the deterministic admin representer lacks the
  pointwise conditional-zero-mean BLY's debiasing requires, so √N(b−b̂)→∞. A
  genuine mathematical dead end, not a scoping choice. Retained only as an
  empirical correction (IDEA-9).
- **IDEA-25 (spatial conformal `scp`, NN-GLS `geospaNN`)** — produce only
  marginal per-cluster intervals (or are Gaussian-only), so they cannot
  instantiate the *joint* area-mean law that is the paper's discriminating
  contribution; they would pad tables without engaging the claim under test.
