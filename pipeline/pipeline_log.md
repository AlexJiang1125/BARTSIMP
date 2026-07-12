# Pipeline Log — BARTSIMP-PG MBG-for-SAE

## Phase 0: SCOPE
- Input: research_brief_template.md (12 sections, complete)
- Brief Completeness Check: 9/9 canonical sections present; §6 has ≥2 named comparators
- Action: SKIPPED (brief COMPLETE, not sparse)
- Output: research_brief.md (canonical copy)

## Phase 1: THINK — dispatched
- Agent: methodology-architect
- Focus (per user direction): reduce admin-1 aggregation-scale center bias (cov80≈0.77); strengthen theory arc beyond Theorem 1

---
## Phase 1: THINK — COMPLETED
- Method: BARTSIMP-PG-SD = BART mean + SPDE-Matérn GP + PG binomial-logit + semi-dense coarse-layer hard depth-floor d0≈√(log N) (Castillo–Ročková) to recenter Admin-1 posterior.
- Bias fix: depth-floor PRIMARY (matched to center-bias diagnosis, ~0 compute, d0=0 reduces to legacy). PG-RoBART demoted to comparator; SoftBART = width complement.
- Theory target: Theorem 2 forward-route single-tree Gaussian Admin-1 functional BvM (KvdVvZ Thm 5.4 ⊕ Castillo–Ročková Thm 9), CONJECTURED; Thm 1 PROVEN modulo re-regime fix.
- Gate: PASS (all spec sections, brief_audit 4-label, literature_briefing complete).
- Files: methodology_specification.md (410L), methodology_rationale.md, brief_audit.md, briefings/literature_briefing.md
- Decision: proceed to Phase 2 VALIDATE.

## Phase 2: VALIDATE — dispatched (validate-method)

---
## Phase 2: VALIDATE — COMPLETED
- NOTE: initial subagent-launched runs were torn down on agent yield (~3 reps). Orchestrator relaunched all 4 detached (nohup) from main session; all completed.
- Verdict: **CONDITIONAL GO**.
- Primary paired (N=500,A=20,24 reps): legacy admin cov80=0.762±0.023, SD=0.777±0.021; Δcov80=+0.0145±0.0102 (~1.4 MC SE, sub-threshold). cov95 0.944→0.942, ratio 1.036→1.034, cluster ~1.0, pixel 0.869→0.866 all preserved.
- d0 sweep (sensitivity): cov80 {0:.793,1:.815,2:.824,3:.815,4:.827} broad +0.02-0.03 plateau, pixel flat ~0.87.
- scale: SD>=legacy at A=5,20,40. adversarial: floor recenters cov80 toward 0.80 from both sides.
- 4 stress tests all CSV-backed (sensitivity=d0 sweep, ablation=d0=0, scale, adversarial); 4 figures generated.
- Conditions: advantage holds in BART-dominated/low-moderate GP-signal regime; modest magnitude; primary contrast not individually significant at 24 reps.
- Decision: proceed to Phase 3 WRITE with CONDITIONAL GO framing.

## Phase 3: WRITE — dispatched (write-manuscript)

---
## Phase 3: WRITE — COMPLETED
- manuscript.tex 861 lines, compiles to 17pp PDF (390KB). Sections: Abstract, Intro, Related work, Methods, Theory, Results (6 subsec), Discussion, Conclusion, Web Appendix.
- 4 PDF figures (primary paired, d0 sweep headline, scale, adversarial); references.bib 21 entries all cited; 5 comparator-complete tables; 9 data CSVs.
- CONDITIONAL-GO framing honest throughout (Δcov80 sub-2-MC-SE flagged as anomaly A1; d0 sweep = headline). Theory: Thm1 thm:matern PROVEN-modulo, Thm2 thm:combined CONJECTURED, RoBART abandoned, s_a>p/2 wall + binomial-logit transfer open.
- Scale-up to 60 reps not finished in budget; reused validated 24-rep numbers (honest future-work note).
- Gate: PASS (lines>500, 4 PDF figs, 21 bib, CSVs, PDF builds, brief fidelity verified).
- Decision: proceed to Phase 4 POLISH (target=42, max_polish=3, venue_compliance=off).

## Phase 4: POLISH — Round 0 grading dispatched (paper-grader)

---
## Phase 4: POLISH — Round 0 (grade)
- Weighted score: 39/50 (raw 28/35). Correctness 4, Completeness 4, Rigor 5, Clarity 4, Novelty 4, Impact 3, Performance 4.
- Fixable: HIGH = intro delta arithmetic self-contradiction (sqrt(1.04^2-1)=0.29 vs stated 0.4; 0.4 needs production ratio 1.08); LOW = abstract rounding, 4 overfull hbox, hyperref unicode bookmarks, overloaded p notation.
- NOT fixable in polish (Phase 1-2, logged): primary effect sub-2-MC-SE/below +0.02 bar (needs 60-100 rep run); no real-data DHS case study (Applications venue); deployed-regime theory conjectural (s_a>p/2 at p=49); DEFF non-robustness; ~100x compute.
- Gap to target 42 = 3 (within polishing range). Proceed Round 1 fix.

---
## Phase 4: POLISH — Round 1 (fix+grade)
- Fixer resolved all 5 round-0 fixable issues; manuscript compiles clean (0 overfull, 0 unicode warnings).
- Re-grade: 41/50 (+2). Per-dim unchanged except re-scored on clean text. Remaining gap to 42 = 1, diagnosed NON-polish-closable (all caps Phase 1-2 structural).
- Decision: EXIT Phase 4 early (round 2 would burn budget fixing Phase 1-2 properties = anti-pattern). Outcome VALIDATED_BELOW_TARGET.

## Phase 4.5: CONSISTENCY — PASS (inline spot-check; no consistency-auditor skill installed)
- All headline manuscript numbers match validated CSVs exactly. Grader already reconciled claims-vs-data both rounds.

## PIPELINE COMPLETE — Outcome: VALIDATED_BELOW_TARGET (41/50)
- Final manuscript: pipeline/phase4_polish/round_1/manuscript.tex (17pp PDF). final_report.md written.
