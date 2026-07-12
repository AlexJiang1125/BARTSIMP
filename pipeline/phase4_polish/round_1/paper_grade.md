# Paper Assessment

### Manuscript: Calibrated Joint Uncertainty Across Scales for Small-Area Estimation: A Bayesian Model-Based Geostatistics Framework with a Semi-Dense Coarse-Layer Tree Prior (BARTSIMP-PG-SD)
### Date: 2026-06-26
### Round: Phase-4 POLISH, round 1 (revised manuscript)

---

## Round-0 fixable-issue resolution check

| # | Round-0 issue (severity) | Status | Evidence |
|---|--------------------------|--------|----------|
| 1 | Introduction δ-arithmetic self-contradiction (HIGH) | **RESOLVED** | Intro now states two regimes, each labeled: validation ratio ≈1.034 ⇒ δ≈0.26, production-surface ratio ≈1.08 ⇒ δ≈0.4 (lines 173–178); §estimand mirrors it (lines 419–421). Arithmetic verifies: √(1.034²−1)=0.265, √(1.08²−1)=0.40. |
| 2 | abstract/body d₀-plateau rounding inconsistent | **RESOLVED** | Abstract "0.79 → 0.815–0.827" (line 92), §d0sweep and tab_d0sweep both "0.815–0.827"; matches CSV rows (0.815, 0.824, 0.815, 0.827). |
| 3 | Overfull \hbox | **RESOLVED** | round_1/manuscript.log: `grep -c "Overfull \hbox"` = 0. |
| 4 | hyperref Unicode warnings | **RESOLVED** | 0 warnings; the sole "Unicode" string in the log is the `puenc.def` package banner, not a warning. |
| 5 | overloaded `p` notation in Theorem 1 | **RESOLVED** | Thm 1 now writes the ill-posedness index as `p_obs=0` and parenthesizes "(the KvdVvZ ill-posedness index, distinct from the covariate dimension p of Section 3.1)" (lines 476–477). |

All five round-0 fixable issues are resolved. PDF compiles clean to 17 pp.

---

### 1. Correctness: 4/5
Every load-bearing equation is covered by the canonical `formula_code_audit.md` (F1–F16, **zero unresolved MISMATCH**), and the manuscript LaTeX matches the reconciled column on spot-check (F12 area functional χ^a = Σ w_i p_i, w_i=n_i/Σn_j; F14 δ=√(ratio²−1)). All headline numbers reproduce against the validated CSVs: primary Δcov80=+0.0145±0.0102 (deltas CSV 0.01447/0.01023), d₀ plateau 0.815–0.827, scale deltas +0.017/+0.034/+0.006, adversarial two-sided recentering — all exact. The δ-arithmetic that was internally contradictory in round 0 is now correct and regime-labeled. Not 5: the headline theory (Thm 2) is conjectured, and Thm 1 is "proven modulo" three stated gaps, so not every boundary case is closed.

### 2. Completeness: 4/5
Theory (two theorems with honest status labels + three propositions), a ten-method benchmark, a primary paired contrast, and three corroborating sweeps (d₀ sensitivity, aggregation-scale, adversarial amplitude) with a built-in ablation (d₀=0 arm). The comparator inventory's FULL IDs (proposed + C1–C5) appear in the primary table; C6–C8 are disclosed as reduced-rep roster context (anomaly A6), not silently dropped. Held at 4, not 5: no real-data DHS case study (a Phase-1/2 gap for an Applications venue — the Kenya/Nigeria deployment is future work), and the two `\begin{definition}`-class estimands are exercised but the production-surface harness is only referenced, not run here.

### 3. Rigor: 5/5
Proper UQ throughout: Monte-Carlo SE on every metric, paired matched-seed design isolating the single component, pre-registered +0.02 bar honestly reported as missed in the primary arm. Heterogeneity (scale sweep), sensitivity (d₀ sweep), and adversarial stress (amplitude sweep) are all present, plus formal theoretical guarantees with explicit PROVEN/CONJECTURED labels. Every Type 1–3 anomaly (A1, A2) is referenced by ID in Discussion/Limitations; A3–A5 (Type 4) also cited. `anomaly_log.md` confirms no Type 1 anomalies. The honest disclosure of sub-2-MC-SE evidence is itself a rigor strength.

### 4. Clarity: 4/5
Well-structured, consistent notation (the round-0 `p` overload is fixed), strong two-cultures motivation, and a precise contribution statement ("calibrated joint UQ across scales" headline + "one prior knob" novelty). The honest-framing paragraph and per-anomaly Discussion read cleanly. Not 5: the prose is dense and qualification-heavy — the Introduction and Theory sections carry long hedged sentences (e.g., the combined-BvM "status, route, and what is confirmed" paragraph) that a tired reviewer must parse twice; elegant economy is not yet there.

### 5. Novelty: 4/5
A genuine cross-disciplinary import: the Castillo–Ročková semi-dense coarse-layer prior (a tree-theory device) is operationalized as a deployable depth-floor d₀≈√(log N) inside a BART+SPDE-Matérn+Pólya-Gamma MBG sampler, and the same prior object is the exact hypothesis the target functional-BvM theorem requires — a non-obvious alignment of deployed method and theory. The Remark distinguishing this prior-side route from the RoBART debiasing route (which provably fails for a deterministic area representer) is a real technical insight. Not 5: the individual ingredients are all known, and the new mathematical element (the s_a>p/2 representer-summability condition) is confirmed only in the single-tree Gaussian sibling.

### 6. Impact: 3/5
The target (LMIC disease mapping / DHS SAE) is a real, decision-relevant domain, and "one fitted model, calibrated UQ at pixel/cluster/area simultaneously" is a workflow-compatible value proposition. But the demonstrated benefit is modest (primary Δcov80=+0.0145, sub-2-MC-SE) and confined to simulation; there is no real-data deployment, the method costs ~2 orders of magnitude more compute than vanilla BART, and the design-effect (DEFF) regime is acknowledged non-robust. A practitioner's decision would change only marginally on present evidence. Honest about all of this, which prevents a lower score, but the demonstrated magnitude caps it at 3.

### 7. Performance: 4/5
On the benchmark the framework is the unique roster member with calibration ratio ≈1 (1.034) at the area scale while competitors are calibrated at a single crossing (GLM 9.47, Vanilla-BART 3.27, GLMM-iCAR 1.36, GLMM-SPDE 1.20), with the joint draws that deliver calibration at every scale — a clear, well-quantified, UQ-backed advantage that holds across the d₀, scale, and amplitude sweeps. Not 5: the *within-family* depth-floor improvement (the paper's specific contribution) is statistically modest in the primary contrast (below 2 MC SE), so the performance win is "clear" rather than "dominant/practically significant."

---

### Raw Score: 28/35
### Weighted Score: 41/50
  Execution (x1.0): Correctness 4 + Completeness 4 + Rigor 5 + Clarity 4 = 17/20
  Research  (x2.0): (Novelty 4 + Impact 3 + Performance 4) x 2.0 = 22/30
  Weighted total = 17 + 22 = 41/50
### Grade: B (>=34)

---

### Code Audit Summary

Per the skill's single-source-of-truth protocol, the canonical audits were VERIFIED (not re-derived):

| Category | Result |
|----------|--------|
| Formula Correctness (`formula_code_audit.md`) | PASS — exists, complete, 16 rows all MATCH, zero unresolved MISMATCH, no NEW_IN_PAPER. Two spot-checks (F12, F14) agree with manuscript LaTeX. |
| Parameter Consistency | PASS — d₀=⌈√log N⌉=3 at N=500 (F11) consistent across §results and tables; production config (m=100, 5000 iter, 2000 burn-in) consistent §sampler/Web-Appendix; replicate counts in table captions (24/16/12/12) match CSV n_rep columns. |
| Comparator coverage (`comparator_inventory.md`) | PASS — every FULL ID (proposed, C1–C5) present in tab_primary; C6–C8 disclosed as reduced-rep context (A6), no silent drop. |
| Anomaly coverage (`anomaly_log.md`) | PASS — no Type 1 anomalies; A1 and A2 (the honesty-critical items) both referenced by ID in Limitations/Discussion; A3–A5 also cited. |

No code-audit-driven score deductions triggered.

### Top 3 Strengths
1. Calibrated-joint-UQ-across-scales is a crisp, defensible headline backed by the calibration-ratio column: the framework is the *unique* benchmark member near ratio 1 at the area scale (Table 1), and the claim is data-verified.
2. Exemplary evidential honesty — the primary contrast is disclosed as sub-2-MC-SE and below the pre-registered bar (A1), every anomaly is addressed by ID, theory carries PROVEN/CONJECTURED labels, and the claim is explicitly transferred to the corroborating sweeps. This is exactly the discipline an Applications venue rewards.
3. Method–theory alignment: the deployed depth-floor (SD1–SD4) is *exactly* the semi-dense hypothesis the target functional-BvM theorem requires, and the RoBART non-transfer remark shows why the prior-side route is the right one.

### Top 3 Weaknesses (with suggested improvement)
1. No real-data DHS case study for an Applications-venue target — caps Impact and Completeness. (Not polish-fixable; Phase-1/2 scope. The planned Kenya/Nigeria deployment would lift Impact to 4.)
2. Headline within-family effect is statistically modest (primary Δcov80 below 2 MC SE) — caps Performance. (Not polish-fixable; needs the planned ≥60–100-rep production run.)
3. Deployed-regime theory is conjectural (Thm 2 CONJ; dimension wall s_a>p/2 fails at p≈49; binomial-logit transfer open) — caps Correctness/Novelty. (Not polish-fixable; live mathematics.)

### Fixable Issues (for paper-fixer)
- **None blocking.** All five round-0 fixable issues are resolved and the build is clean.
- Minor, optional clarity polish (would not change the grade): the combined-BvM "Status, route, and what is confirmed" paragraph (lines 511–526) and the honest-framing paragraph (lines 195–213) are long and qualification-dense; splitting one or two sentences would aid first-pass readability. This is cosmetic, not required.

### "So What?" Test
A week later I would remember two things: (1) the calibration-ratio framing that exposes every competitor as calibrated-at-one-scale, and (2) that a single √(log N) tree-prior depth-floor — borrowed from CART theory and provably the theorem's own hypothesis — recenters the one soft spot at zero compute cost. The modest primary effect would also stick, but as a candid limitation rather than an overclaim.

### Verdict
A well-executed, unusually honest methods paper with a clean cross-disciplinary contribution and verified, internally-consistent numbers; the round-1 polish closed every fixable round-0 defect. Against top-tier statistics-journal standards (JASA A&CS / AOAS) it is a **revise**: the score caps are structural (no real-data deployment, conjectural deployed-regime theory, sub-2-MC-SE primary effect, ~100x compute, DEFF non-robustness) and are Phase-1/2 properties, not polish defects. Weighted 41/50 (B), one point below the 42 target; the gap is not closable by further prose editing and would require the planned production run and/or real-data case study.
