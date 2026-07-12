# Literature Briefing — BARTSIMP-PG-SD (Phase 1 → Phase 3 digest)

> Synthesis digest produced by the methodology-architect phase for the
> downstream writing/positioning phase (Phase 3). This is a **revision-mode**
> briefing: the relevant literature is already digested in the in-repo theory
> drafts and the research brief's reference list, so no fresh paper-reader
> subagents were dispatched (that would risk citation drift against the
> already-curated set). The tables below consolidate what each source provides,
> what the method still needs, the positioning seeds, and the coverage gaps the
> writing phase must close.

## Papers Read (curated set — sources already in-program)

| # | Citation (short) | Role in this paper | Where digested |
|---|------------------|--------------------|----------------|
| P1 | **Castillo & Ročková (2021)**, UQ for Bayesian CART | **Load-bearing.** Semi-dense coarse-layer prior → functional BvM; the depth-floor is its deployment | `forward_bvm_exploration.md`, brief §11 |
| P2 | **Kuelbs–van der Vaart–van Zanten (2011)**, functional BvM for GP | Thm 1 (Matérn-only coverage); conditioning-on-trees route for Thm 2 | `draft_thm2_single_tree.md`, HANDOFF |
| P3 | **Polson, Scott & Windle (2013)**, Pólya–Gamma | The binomial-logit augmentation; conditional-Gaussian structure | `dev/bartsimp_pg.R`, brief §3 |
| P4 | **Lindgren, Rue & Lindström (2011)**, SPDE/GMRF Matérn | The spatial field representation; PC priors | brief §3 |
| P5 | **Chipman, George & McCulloch (2010)**, BART | The mean ensemble + branching prior that the floor modifies | base model class |
| P6 | **BLY25 (Bayes-debiased one-step / RoBART)** | **Dead theory route** (Mistake 13); empirical comparator arm only | HANDOFF, brief §9 |
| P7 | **Linero & Yang (2018)**, SoftBART | Orthogonal width lever (complement, not headline) | soft-vs-hard sim, brief §9 |
| P8 | **Disaggregation / TMB areal binomial** (Nandi et al.) | The fair native joint-UQ comparator | brief §6 |
| P9 | **RF-GLS, SPARforest** (spatial ML) | Constant-bias / laundered-width camp comparators | brief §6 |
| P10 | **GLMM-iCAR, GLMM-SPDE (INLA)** | Disease-mapping / LMIC-standard linear-mean baselines | brief §6 |

## Per-Paper Takeaways

- **P1 Castillo–Ročková (2021).** The single most important external result.
  Tree-prior functional under-coverage is a *coarse-layer pruning* phenomenon;
  the fix is a prior modification (semi-dense, j₀ ≍ √(log n)), not debiasing.
  Provides the mechanism, the theorem template (Thm 9 marginalization), and the
  exact form of the depth-floor. **Limit:** their result is single-tree,
  Gaussian, no spatial GP, no survey aggregation, no binomial — all of which the
  paper must extend or mark open.
- **P2 KvdVvZ (2011).** Functional BvM for GP-priors with a representer-decay
  condition. Provides Thm 1 directly and, via condition-on-trees, the GP half of
  Thm 2. **Watch:** regime selection (Thm 5.3 conservative vs 5.4 nominal) — the
  current draft is mis-regimed and must be corrected to the nominal statement.
- **P3 PG (2013).** Renders the model conditionally Gaussian; the surface study
  confirms point-estimate invariance across PG/AC/HH (width-only difference).
  Provides the likelihood-side machinery and the Gaussian-sibling bridge.
- **P4 SPDE (2011).** Standard; provides the sparse-precision Matérn field and
  PC priors. No open issues.
- **P5 BART (2010).** Provides the mean ensemble and the Chipman branching prior
  P(split|d) = α(1+d)^{−β} that the depth-floor edits. The non-Donsker
  slowly-vanishing mean bias is the BART-specific phenomenon the paper diagnoses.
- **P6 BLY25/RoBART.** *Do not cite as a theory route for the admin functional.*
  Its debiasing requires a pointwise conditional-zero-mean the deterministic
  admin representer lacks (Mistake 13). Cite only as the empirical correction
  arm we benchmark the prior-side fix against.
- **P7 SoftBART.** Narrows bands (≈ −7% post_sd) but does not move center bias.
  Position as an orthogonal, stackable width lever — not a competitor to the
  floor.
- **P8 disaggregation.** The honest native analog with a genuine joint law; the
  "averageable-misfit camp" that is calibrated only at the coarse crossing scale.
- **P9 RF-GLS / SPARforest.** The "constant-bias / laundered-width camp,"
  calibrated only at the fine crossing scale; the cautionary cases for
  scale-dependent coverage.
- **P10 GLMM-iCAR / GLMM-SPDE.** The field-standard linear-mean baselines; the
  comparison that shows what the flexible BART mean buys.

## Provides / Needs Matrix (reference)

| Need (for the contribution) | Provided by | Status |
|-----------------------------|-------------|--------|
| Coarse-layer recentering mechanism | P1 (semi-dense prior) | PROVEN (their setting); deployed as depth-floor |
| GP-only Admin-1 coverage | P2 (functional BvM) | PROVEN modulo re-regime fix |
| Tree ⊕ GP product-domain BvM | P1 ⊕ P2 (condition-on-trees) | CONJECTURED (Thm 2, no fatal holes) |
| Representer-decay condition | P1 Thm 9 + own P2 computation | CONFIRMED: s_a > p/2 |
| Binomial-logit deployment | P3 (PG) | PROVEN (sampler); coverage transfer CONJECTURED |
| Spatial field | P4 (SPDE) | PROVEN/standard |
| Flexible mean + non-Donsker diagnosis | P5 (BART) | mean PROVEN; bias diagnosis = own contribution |
| Bias-reduction theory route | NOT P6 (dead); P1 forward route | viable route = P1⊕P2; P6 rejected |
| Width control (secondary) | P7 (SoftBART) | complement |
| Joint-law comparators | P8, P9, P10 | comparator roster |

## Positioning Seeds (for the Introduction / Related Work)

1. **"Calibrated joint UQ at every aggregation scale"** is the headline claim;
   frame the paper as a *framework*, not a Gaussian-paper sequel. The novelty is
   the scale-invariant calibration property + its prior-theoretic basis, not the
   ingredients (BART, GP, PG are all known).
2. **The two-camp framing** (constant-bias: SPARforest/RF-GLS, calibrated only at
   the fine crossing; averageable-misfit: GLMM-SPDE/disaggregation, calibrated
   only at the coarse crossing) is the rhetorical spine that proves "calibrated
   at *every* scale" — no competitor holds ratio ≈ 1 across all A.
3. **Diagnosis-then-mechanism** narrative: the cov95-vs-cov80 + ratio
   decomposition identifies a non-Donsker *center* bias loading on the coarse
   tree component; the semi-dense floor is the matched, Bayesian-coherent remedy.
   This is a new diagnostic read, not just a method.
4. **Theory/empirics alignment** as a selling point: the deployed knob (semi-dense
   floor) is *exactly* the hypothesis of the forward-route theorem — a rare
   coincidence worth foregrounding.
5. **Honest scope** as a credibility move: state the binomial-transfer gap, the
   s_a > p/2 dimension wall, and the DEFF non-robust regime up front; the rigorous
   theorem is the Gaussian single-tree sibling.

## Coverage Gaps for Phase 3 (what the writing phase must close)

- **Re-regime Theorem 1** to the nominal-BvM (KvdVvZ Thm 5.4) statement and close
  three NEEDS-DETAIL items (boundary curvature, n→N_a rate, width constant C_max).
- **Finish Theorem 2 live math:** E-i (tree ⊕ GP product-domain extension),
  E-ii (in-sample V_smp = σ²/π_a ↔ population reconciliation via C–R Step 9),
  E-iii (dyadic vs CGM axis-split). Mark CONJECTURED until done.
- **State the binomial transfer** as an explicit open gap; do not let the
  Gaussian-sibling theorem's strength leak into an over-claim about the deployed
  binomial model.
- **Empirical confirmation pending:** the depth-floor A/B has not yet been run;
  the spec's MVE defines it. Phase 3 prose must condition claims on that result
  (and prepare the null-result mechanism-adjudication framing).
- **Citation hygiene:** ensure BLY25/RoBART is cited *only* as an empirical
  comparator and never as a theory route for the admin functional, to avoid
  reintroducing the Mistake-13 error into the manuscript.
