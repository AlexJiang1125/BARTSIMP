# Methodology Rationale — BARTSIMP-PG-SD

> Phase 1 (methodology-architect) companion to `methodology_specification.md`.
> Records *why* each design decision was made: the literature basis, the
> combination reasoning, the alternatives considered and rejected, the impact
> justification, the supporting stress-test evidence, and the open questions.
> This is a **revision-mode** methodology: the BARTSIMP-PG sampler already
> exists and is validated; the contribution is a single prior-side extension
> (the semi-dense depth-floor) plus a strengthened theory arc, framed as a
> calibrated-joint-UQ MBG-for-SAE framework paper.

## 1. The decision that organizes everything: recenter, do not rewiden

The whole design follows from one diagnostic read of the existing scorecard
(`sim_dgp5_coverage`): BARTSIMP-PG under-covers at Admin-1 (cov80 ≈ 0.767 vs
0.80) but is near-nominal at cov95 (≈ 0.941) and over-covers at the fine
(pixel/cluster) scale. The calibration ratio RMSE/post_sd ≈ 1.08 with bias
≈ 0.0017.

- If the deficit were a **width** problem, cov95 would lag too and the ratio
  would be ≫ 1 driven by post_sd being too small everywhere. Instead cov95 is
  fine and only the *tighter* interval (cov80) misses.
- The standardized center bias implied is δ = √(1.08² − 1) ≈ 0.41 SD — a real
  off-center shift of the area-mean posterior, small enough that the wide 95%
  band still catches the truth but large enough to eat into the 80% band.

**Decision:** the lever must move the posterior *center*, not its *width*. This
immediately (a) selects prior-side recentering mechanisms over width knobs
(SoftBART, k_BART, DEFF-tempering all become secondary), and (b) imposes the
hard design constraint that whatever we do must not widen the already-conservative
fine scale — i.e. the lever must act *only on the coarse partition that area
averaging loads on*.

## 2. Literature basis for the chosen mechanism

**Castillo–Ročková (2021), "Uncertainty quantification for Bayesian CART."**
This is the load-bearing citation. Their result: a single Bayesian regression
tree's posterior for a *linear functional* under-covers precisely because the
Galton–Watson tree prior **prunes the coarse (low-resolution) layers**, so the
coarse Haar/CART coefficients — the ones a smooth functional loads on — are
shrunk toward zero and the credible interval is centered off-truth. Their
remedy is a **"semi-dense" prior**: force the inclusion indicators γ_{ℓk} = 1
for all cells at all levels ℓ ≤ j₀(n) with j₀(n) ≍ √(log n), keeping the
coarse layers complete while leaving deep layers adaptive. Under semi-dense
priors they prove a functional Bernstein–von Mises (BvM) theorem.

This maps onto our problem object-for-object: the BART non-Donsker mean bias is
a slowly-vanishing low-frequency bias; area averaging is a smooth linear
functional that loads on exactly the coarse component; the Chipman branching
prior P(split | depth d) = α(1+d)^{−β} is a pruning prior that under-represents
coarse layers. The depth-floor d₀ ≈ √(log N) is the literal deployment of C–R's
semi-dense hypothesis with j₀ ↔ d₀.

**Kuelbs–van der Vaart–van Zanten (2011) functional BvM machinery** supplies
the GP-side coverage result (Theorem 1) and the conditioning-on-trees route for
Theorem 2: conditional on a tree topology T, η | T is a single Gaussian process
on the product domain X × D, so KvdVvZ Prop 3.2 / Thm 5.4 apply conditionally,
and C–R Thm 9 carries the marginalization over T.

**Polson–Scott–Windle (2013) Pólya–Gamma** and **Lindgren–Rue–Lindström (2011)
SPDE/GMRF** are the inherited, already-validated pieces of the base model class;
they are not the contribution but are restated for completeness.

## 3. Why combine these specific pieces (and not others)

The combination is BART mean + SPDE-Matérn GP + PG binomial likelihood + a
semi-dense depth-floor. Reasoning:

- **BART + GP is non-negotiable for the headline.** The spatial ablation
  (Vanilla-BART, GP removed) collapses Admin-1 cov80 to ≈ 0.28: the GP is what
  carries cross-cluster covariance into the area-mean law. The framework's
  discriminating property (a genuine *joint* area-mean posterior) requires both
  a flexible covariate mean *and* a structured spatial field drawn jointly.
- **PG is the right augmentation** because it renders the model conditionally
  Gaussian, which (i) lets the GP step reuse the exact Woodbury/Cholesky machinery
  and (ii) makes the Gaussian-likelihood theory sibling a faithful anchor for the
  binomial deployment. The surface study confirmed PG / Albert–Chib / Holmes–Held
  give an *identical point estimate* (only band width differs), so the choice of
  augmentation is a width knob, not an identification choice — PG is kept for its
  Gaussian-conditional convenience.
- **The depth-floor is the only addition**, deliberately minimal: it is the
  single object that is simultaneously (a) the empirically-indicated recentering
  lever, (b) the exact hypothesis of the viable forward-route theorem, and (c) a
  pure-R prior edit honoring the "do not rewrite the sampler in compiled code"
  non-goal.

## 4. Alternatives considered and why rejected/demoted

| Candidate | Disposition | Reason |
|-----------|-------------|--------|
| **PG-RoBART / BLY25 debiasing** (binomial-logit port of BLY25 Algorithm 1) | **Demoted to empirical comparator arm only** | Its *theory* route to the Admin-1 BvM is **fatally holed (Mistake 13/O1)**: BLY's one-step debiasing needs a pointwise conditional-zero-mean that the random missing-data indicator R satisfies but the **deterministic** admin representer 1_a/π_a does not, leaving a Θ(ε_N) bias ⇒ √N(b−b̂)→∞ in the non-Donsker regime. Confirmed by two independent audits. It also adds a *new estimator* (not Bayesian-coherent) and is not bit-verifiable against the base sampler. Kept only as an empirical correction to benchmark against. |
| **SoftBART / soft trees** (Linero–Yang) | **Complement (width lever)** | Orthogonal: it narrows bands (≈ −7% post_sd in paired sims) but the paired soft-vs-hard run showed **no center-bias change** (dgp5 bias +0.0014 → +0.0013). It addresses width, which is not the problem. It *stacks* with the floor (can be run on top) but cannot be the headline. |
| **k_BART amplitude / β_tree decay tuning** | Rejected as primary | These are global width/depth knobs that move the whole partition; they cannot recenter the coarse component without disturbing the fine scale, violating the design constraint. |
| **DEFF-tempered likelihood (Route 3)** | Retained as a robustness arm, not the fix | It inflates post_sd ~ √DEFF to address the overdispersion regime, a *variance* remedy for a *different* (DEFF) weakness — not the center-bias problem. |
| **Cluster random effect ξ_i** | Default off | BART tends to absorb the RE (seen on Nigeria + BB-RE arm), so it buys little and costs identifiability; available but off by default. |
| **scp (spatial conformal), geospaNN (NN-GLS)** | Not carried | Emit only *marginal* per-cluster intervals (or are Gaussian-only); they cannot instantiate the *joint* area-mean law that is the contribution's discriminating axis, so they would pad tables without testing the claim. |

## 5. Impact justification — why the chosen mechanism is worth it

- **It is the right size of intervention.** The Admin-1 gap is only 0.03 in
  cov80; this does not call for a new estimator or a re-architected sampler, it
  calls for a targeted prior nudge on the coarse partition. The depth-floor is
  exactly that: three localized edits to `rbart.R` (complete-tree init; forbid
  prune at depth ≤ d₀; guard the death move), ≈ 0 extra compute, and d₀ = 0
  recovers the legacy sampler byte-for-byte. This is the minimal change that can
  move the diagnosed quantity.
- **Theory/empirics alignment is rare and valuable.** The same object that the
  experiment toggles (the semi-dense floor) *is* the hypothesis the
  forward-route theorem needs. Most methods papers bolt a heuristic fix onto a
  model and prove a theorem about an idealized cousin; here the deployed knob and
  the proven hypothesis coincide. That is a genuine selling point for a top venue.
- **The framework story survives regardless of the floor's empirical verdict.**
  The headline is calibrated-joint-UQ-across-scales, which BARTSIMP-PG already
  largely delivers (ratio 1.01 / 1.08 / 0.96 across A = 40/20/5). The floor
  pushes the one soft spot onto the line; even a *null* floor result is publishable
  as a mechanism adjudication (it would redirect the lever to the ensemble-CLT-
  shrinkage story, τ² ∝ m^{−1/2}).

## 6. Supporting stress-test / empirical evidence (from the existing program)

These are the numbers that justify the design choices (sources:
`sim_dgp5_coverage/README.md`, `sim_mbg_surface/README.md`, `theory/HANDOFF.md`):

- **Headline scorecard (dgp5, admin-1):** BARTSIMP-PG cov80 0.767, cov95 0.941,
  RMSE 0.0310, post_sd 0.0287, bias 0.0017, ratio 1.08. → the center-bias
  diagnosis.
- **Cross-A sweep:** ratio 1.01 (A=40) / 1.08 (A=20) / 0.96 (A=5) → the
  scale-invariant-calibration claim, with A=20 as the soft spot.
- **Spatial ablation:** removing the GP collapses admin cov80 to ≈ 0.28 → the GP
  is essential to the joint area law.
- **Soft-vs-hard paired:** soft −7% post_sd, bias unchanged (+0.0014 → +0.0013)
  → SoftBART is a width lever, not a recentering lever.
- **DEFF stress (ρ=0.2):** BARTSIMP variants worst-covered (ratio ≈ 1.84) → the
  acknowledged non-robust regime; excluded from the robustness requirement.
- **Pixel surface:** cov80 0.976 (over-covers); freq width-calibration ratio → 1
  only as N → ∞ → the fine-scale conservatism is a feature to *preserve*.
- **Augmentation invariance:** PG / AC / HH identical point estimate, width-only
  difference → augmentation choice is a width knob.
- **Theory state:** Thm 1 essentially rigorous but mis-regimed (cite KvdVvZ 5.4
  nominal, not 5.3 conservative; β⋆ > ν/2 + 1/4); single-tree Gaussian Thm 2
  audited "MAJOR REVISION, no fatal holes" with P2 (s_a > p/2) confirmed; BLY
  route dead.

## 7. Open questions carried forward (for Phase 3 / writing and Phase theory)

1. **Mechanism adjudication.** Does the floor actually move cov80, or does the
   competing ensemble-CLT + coarse-shrinkage story (τ² ∝ m^{−1/2}) own the dip?
   The paired A/B decides; both outcomes are reportable.
2. **d₀ calibration.** Is ⌈√(log N)⌉ the right floor, or does the optimal d₀
   depend on the coarse structure of the truth? Sensitivity sweep {2,3}.
3. **Theorem 2 remaining math.** E-i (tree ⊕ GP product-domain extension),
   E-ii (in-sample ↔ population representer reconciliation via C–R Step 9),
   E-iii (dyadic vs CGM axis-split). None fatal, all live.
4. **Theorem 1 fix.** Re-regime to the nominal-BvM (Thm 5.4) statement and close
   the three NEEDS-DETAIL items (boundary curvature, n→N_a rate, width constant).
5. **Dimension wall.** s_a > p/2 (p=49 ⇒ s_a > 24.5) has no single-tree rescue;
   does BART's additive-sieve / low-effective-dimension structure beyond one tree
   provide one? UNKNOWN.
6. **Binomial transfer.** The conservative/nominal coverage proof for the
   deployed binomial-logit model is the standing open gap; the Gaussian sibling
   is the anchor and the paper must mark this open rather than over-claim.
