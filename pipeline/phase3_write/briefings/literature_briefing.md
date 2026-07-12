# Literature Briefing (Phase 3, writing/positioning) — BARTSIMP-PG-SD

> Writing-focused synthesis for the Paper Writer. Builds INCREMENTALLY on
> `pipeline/phase1_think/briefings/literature_briefing.md` (the curated
> in-program set P1–P10, already digested in the theory drafts and brief §12).
> No fresh deep readers were dispatched: the reference set is fixed and curated,
> and re-reading risks citation drift against the already-anchored sources.
> This briefing adds the positioning layer (comparative table, gap analysis,
> quality benchmarks, related-work prose, action items).

## Executive Summary

Model-based geostatistics (MBG) for small-area estimation (SAE) is dominated by
linear-mean spatial models — GLMM-iCAR (Besag) and the INLA SPDE–Matérn
workhorse (Lindgren–Rue–Lindström 2011; Rue et al. 2009) — and by
disaggregation regression (binomial areal + TMB joint UQ). These calibrate the
area functional only at the coarse crossing scale because a linear mean cannot
represent curvature/interactions. Flexible-mean alternatives (BART, Chipman–
George–McCulloch 2010; spatial random forests RF-GLS, Saha–Basu–Datta 2023;
SPARforest) either lack a spatial term that aggregates (vanilla BART collapses
at the area functional) or hit nominal coverage only by laundered width / a
constant bias that makes coverage scale-dependent. Our framework couples a BART
mean, an SPDE–Matérn GP, and Pólya–Gamma augmentation (Polson–Scott–Windle
2013) for the binomial-logit likelihood, delivering coherent JOINT uncertainty
at pixel/cluster/area scales — the property no competitor holds at every scale.
The single methodological novelty beyond this combination is a semi-dense
coarse-layer tree prior (a hard depth-floor d₀≈√(log N), after Castillo–Ročková
2021) that recenters the one soft spot — Admin-1 cov80 — toward nominal without
disturbing fine-scale calibration, and whose hypothesis is exactly that of the
forward-route functional-BvM theorem we target (Knapik–van der Vaart–van Zanten
2011 ⊕ Castillo–Ročková 2021). The contribution is the scale-invariant
calibration property and its prior-theoretic basis, not the ingredients.

## Comparative Table

| Aspect | GLMM-SPDE (INLA) | disaggregation | RF-GLS / SPARforest | Vanilla-BART | Castillo–Ročková 2021 | Ours (BARTSIMP-PG-SD) |
|--------|------------------|----------------|---------------------|--------------|------------------------|------------------------|
| Mean model | linear | linear | forest (non-spatial-aware mean) | BART | single CART tree | BART ensemble |
| Spatial term | Matérn SPDE | areal/CAR | GLS spatial corr. | none | none | Matérn SPDE |
| Likelihood | binomial-logit | binomial | Gaussian (transformed) | probit/logit | Gaussian | binomial-logit (PG) |
| Joint area-mean law | yes | yes (TMB) | no (marginal only) | degenerate (no field) | n/a (no aggregation) | **yes (joint field draws)** |
| Calibrated at every scale | no (coarse only) | no (coarse only) | no (fine only / scale-dep.) | no (collapses) | n/a | **yes (the claim)** |
| Coverage recentering knob | — | — | — | — | semi-dense prior (theory) | **deployed depth-floor d₀** |
| Theory | BvM (linear) | — | — | contraction only | forward functional BvM (single tree) | Thm 1 (GP) proven; Thm 2 (combined) conjectured |

## Gap Analysis (PRESENTATION gaps — method is validated)

- **Gap: real-data deployment figure.** Reference MBG papers (INLA/SPDE, dis-
  aggregation) show posterior-median maps + CI-width maps on DHS-style data.
  Priority: important. We frame Kenya DHS 2014 child wasting as the motivating
  application; the empirical core is simulation-led (the honest scope of this
  submission). A real-data map is described as the deployment target and a Web
  Appendix placeholder; full real-data fit is future production work.
- **Gap: two-camp scale-crossing figure.** The "constant-bias camp
  (RF-GLS/SPARforest, calibrated at the fine crossing) vs averageable-misfit camp
  (GLMM-SPDE/disaggregation, calibrated at the coarse crossing)" narrative is the
  rhetorical spine; it should be made visible via the calibration-ratio column
  read across A. Priority: critical for positioning.
- **Gap: theorem count vs an AOAS/JASA theory paper.** Reference theory papers
  (KvdVvZ; Castillo–Ročková) state 2–4 theorems. We carry Theorem 1 (proven
  modulo re-regime), the single-tree Gaussian Theorem 2 (conjectured, no fatal
  holes), the cluster-floor result, and the dimension-wall proposition. Priority:
  important — honesty labels (PROVEN/CONJECTURED) are themselves a credibility move.

## Quality Benchmarks (to match the target venue)

- Sections: Abstract, Introduction, Related Work, Methods, Theoretical
  Properties, Results (multi-experiment), Discussion (5 topics), Limitations,
  Web Appendix. (~8–10 numbered sections.)
- Tables: one per evaluation type (primary paired, d₀ sweep, scale, adversarial),
  each carrying the full comparator roster; a calibration-ratio diagnostic column.
- Figures: ≥2 publication PDF (d₀ sweep calibration curve; primary paired;
  scale; adversarial recentering).
- Theory: ≥2 stated results with explicit hypotheses and honest status labels;
  the s_a > p/2 representer computation; the in-sample vs superpopulation estimand
  distinction.
- Reproducibility: seeded paired designs (truth drawn first → byte-identical
  across arms), Framework-R run note, validated-code provenance.

## Related Work Draft (seed prose for the Writer)

Small-area estimation under complex survey designs has two cultures. The
*design-based* culture (direct/Hájek estimators, model-assisted GREG, and
recently prediction-powered inference) is correct where its iid-residual variance
model holds but provides no pixel- or cluster-level uncertainty and understates
standard errors under design effects. The *model-based geostatistics* culture
(Diggle–Tanton–Moyeed; Rue et al. 2009; Lindgren et al. 2011) places a latent
Gaussian field over space and aggregates it to areas; the INLA SPDE–Matérn model
is the LMIC standard and disaggregation regression (TMB) its native binomial,
joint-UQ analog. Both rest on a *linear* mean, so where the covariate–response
relationship is nonlinear or non-additive the area functional inherits a
mean-misfit bias that surfaces as scale-dependent under-coverage. Flexible-mean
machine learning enters from two directions. Bayesian additive regression trees
(Chipman et al. 2010) give a calibrated-pointwise nonparametric mean but, absent
a spatial term, collapse at the area functional. Spatial random forests (RF-GLS;
SPARforest) restore a spatial mean but emit only marginal per-cluster intervals
and reach nominal coverage by a constant bias or laundered width that holds at
one aggregation scale and fails at others. None delivers a coherent *joint*
area-mean law calibrated across scales. On the theory side, Knapik–van der
Vaart–van Zanten (2011) give a functional Bernstein–von Mises theorem for
Gaussian-process priors, and Castillo–Ročková (2021) prove that tree-prior
functional under-coverage is a coarse-layer-pruning phenomenon remediable by a
*prior modification* (a semi-dense, dense-coarse-layers prior) rather than by
estimator-side debiasing. Our framework imports the latter as a deployable
depth-floor and assembles the former two results, by conditioning on trees, into
the forward route for a combined functional BvM — distinct from the
RoBART/Breunig–Liu–Yu (2025) debiasing route, which we show does not transfer to
the deterministic area representer. SoftBART (Linero–Yang 2018) supplies an
orthogonal width lever we note but do not feature.

## Flags for Authors

- **Novelty boundary (Castillo–Ročková 2021).** Their semi-dense result is the
  closest prior art for the recentering mechanism. State explicitly that theirs
  is single-tree, Gaussian, no spatial GP, no survey aggregation, no binomial —
  so our deployment (ensemble + GP + PG + area functional) is not a special case
  of their theorem; our combined statement is conjectured, not inherited.
- **Citation hygiene (BLY25/RoBART).** Cite Breunig–Liu–Yu (2025) ONLY as an
  empirical debiasing comparator, NEVER as a theory route for the admin
  functional. The deterministic representer lacks the pointwise conditional
  zero-mean their debiasing requires (Mistake 13); reintroducing it would import
  a known error.
- **Theorem 1 re-regime.** The Matérn-only result must cite KvdVvZ Thm 5.4
  (nominal BvM, q=1/4 > p=0) not Thm 5.3 (conservative, requires q<p); the
  β⋆ > ν/2 + 3/4 vs + 1/4 flag is the symptom.

## Action Items for Writer

1. Emphasize the **scale-invariant calibration** property and the **two-camp**
   framing as the contribution; treat BART/GP/PG as known ingredients.
2. Position the depth-floor as a **single Bayesian-coherent prior knob** that is
   *the same object* the forward-route theorem requires (theory/empirics alignment).
3. Carry honest status labels on every theorem; state the binomial-transfer gap,
   the s_a > p/2 dimension wall, and the DEFF non-robust regime up front.
4. Use the calibration-ratio RMSE/post_sd ("earned vs laundered") as the
   adjudication instrument in every results table.
5. Report the CONDITIONAL-GO modest-effect honestly: d₀ sweep is the headline;
   primary paired Δcov80 < 2 MC SE at 24 reps is a stated limitation.
