# Research Pipeline Final Report — BARTSIMP-PG-SD (MBG-for-SAE)

## Outcome: VALIDATED_BELOW_TARGET (41/50; target 42)
## Total Phases Executed: 4 (+ Phase 4.5 consistency spot-check)
## Rethinks Used: 0

The 1-point gap to target is **structural, not polish-closable**: every score cap is a
Phase-1/2 property (no real-data DHS case study for an Applications venue; conjectural
deployed-regime theory; modest, sub-2-MC-SE primary effect; ~100× compute; DEFF
non-robustness) — exactly the limitations the brief anticipated. Reaching 42 requires
the planned production-scale run (≥60–100 reps) and/or a real-data deployment, not
further prose editing.

---

## Phase 1: THINK
- **Methodology**: BARTSIMP-PG-SD — BART mean + SPDE–Matérn GP + Pólya–Gamma binomial-logit
  likelihood, extended with a **semi-dense coarse-layer hard depth-floor** d₀≈√(log N)
  (Castillo–Ročková) that recenters the Admin-1 posterior to fix aggregation-scale
  under-coverage without disturbing fine-scale calibration.
- **Comparators**: GLM, Vanilla-BART, GLMM-iCAR, GLMM-SPDE, disaggregation, RF-GLS,
  SPARforest, BARTSIMP-PG (legacy d₀=0); calibration-ratio "earned vs laundered" diagnostic.
- **Key design decisions**: depth-floor chosen as PRIMARY (matched to the center-bias
  diagnosis: cov95 near-nominal, cov80 lagging ⇒ bias not width; ~0 compute; reduces to
  legacy at d₀=0). PG-RoBART demoted to an empirical comparator (theory route fatally holed).
  SoftBART kept as a width-only complement.

## Phase 2: VALIDATE
- **Verdict**: CONDITIONAL GO.
- **Primary advantage**: Admin-1 cov80 0.762±0.023 (legacy) → 0.777±0.021 (SD);
  paired Δcov80 = +0.0145±0.0102 (~1.4 MC SE — directional, below the +0.02 bar).
- **Headline calibration evidence**: d₀ sweep — cov80 rises to a broad plateau
  {0:0.793, 1:0.815, 2:0.824, 3:0.815, 4:0.827}, pixel flat ~0.87, cov95 ~0.94 throughout.
- **Corroboration**: SCALE (SD≥legacy at A=5,20,40); ADVERSARIAL (floor recenters cov80
  toward 0.80 from both sides). Fine-scale calibration (cov95, ratio, cluster~1.0, pixel)
  preserved, not improved.
- **Key limitation**: primary contrast not individually significant at 24 reps; magnitude modest.
- **Stress tests**: sensitivity (d₀ sweep), ablation (d₀=0 arm), scale, adversarial — all CSV-backed.
- **Iterations**: depth-floor edit → smoke test → 4 paired runs → finalize. Rethinks: 0.
- **Operational note**: initial subagent-launched runs were torn down on agent yield
  (~3 reps); orchestrator relaunched all 4 detached (nohup) from the main session — all completed.

## Phase 3: WRITE
- **Manuscript**: `pipeline/phase4_polish/round_1/manuscript.tex` (final; 871 lines), source
  originated at `pipeline/phase3_write/manuscript.tex` (861 lines). Compiles clean to 17pp PDF.
- **Figures**: 4 PDF (primary paired coverage; d₀-sweep headline; scale; adversarial).
- **Tables**: 5, all carrying the full comparator roster + calibration-ratio column.
- **Evaluation types**: paired A/B, d₀ sensitivity sweep, scale, adversarial.
- **Theory**: Theorem 1 (Matérn-only Admin-1 BvM) PROVEN modulo a stated KvdVvZ Thm 5.4
  re-regime fix; Theorem 2 (single-tree Gaussian combined BvM) CONJECTURED, no known fatal
  holes; RoBART route abandoned and disclosed; s_a>p/2 dimension wall (p≈49) and binomial-logit
  transfer stated open.

## Phase 4: POLISH
- **Starting score**: 39/50 (round 0).
- **Final score**: 41/50 (round 1).
- **Rounds**: 1 (exited early; remaining gap diagnosed non-polish-closable).
- **Issues fixed**: δ-arithmetic intro contradiction (now regime-labeled: validation 1.034⇒δ≈0.26,
  production 1.08⇒δ≈0.4); abstract/body d₀-plateau rounding; 4 overfull \hbox → 0;
  hyperref Unicode bookmark warnings → 0; disambiguated overloaded `p` in Theorem 1.

## Score Breakdown (round 1, final)
| Dimension    | Score |
|--------------|-------|
| Correctness  | 4     |
| Completeness | 4     |
| Rigor        | 5     |
| Clarity      | 4     |
| Novelty      | 4     |
| Impact       | 3     |
| Performance  | 4     |
| **Raw Total**      | **28/35** |
| **Weighted Total** | **41/50** |

## Phase 4.5: Consistency
- Spot-check (inline; no dedicated consistency-auditor skill installed): all headline numbers
  in the final manuscript match the validated CSVs exactly (cov80 0.762/0.777, Δ 0.0145,
  cov95 0.944/0.942, d₀ plateau 0.815–0.827, ratios 1.034/1.08). Paper-grader had already
  reconciled claims vs data in both rounds (formula_code_audit 16/16 MATCH, anomaly_log A1–A6).
  Verdict: PASS.

## Not-Fixable-In-Polish (Phase 1-2 caps — logged, honest)
1. Primary paired effect sub-2-MC-SE and below the +0.02 bar → needs ≥60–100-rep production run.
2. No real-data DHS (Kenya wasting) case study despite the Applications-venue target.
3. Deployed-regime theory conjectural; s_a>p/2 cubic-rate condition fails at p≈49.
4. Admin-level design-effect / overdispersion non-robustness (BART absorbs the cluster RE).
5. ~2-orders-of-magnitude compute cost vs vanilla BART.

## Lessons Learned
- **What worked**: the d₀-sweep is a cleaner, more persuasive calibration figure than the
  single primary paired contrast — robust broad plateau, no fine-scale cost. The honest
  CONDITIONAL-GO framing was credited as Rigor (5/5), not penalized.
- **What didn't / risk**: subagent-launched background R jobs are killed when the agent yields —
  long runs MUST be launched detached (nohup/setsid) from a durable session. This cost one
  ~1-hour false start before relaunch.
- **Next run**: (a) production-scale paired run (≥60–100 reps) to push Δcov80 above 2 MC SE;
  (b) add the Kenya DHS real-data deployment to lift Impact (the single biggest score cap for
  an Applications venue); (c) progress Theorem 2 from CONJECTURED toward proven (or strengthen
  the Gaussian single-tree anchor) to lift Rigor/Novelty.

## All Files Produced
- `pipeline/phase1_think/`: methodology_specification.md, methodology_rationale.md,
  brief_audit.md, briefings/literature_briefing.md
- `pipeline/phase2_validate/`: validation_report.md, iteration_log.md, comparator_inventory.md,
  validated_code/{rbart depth-floor edit, run_paired_dgp5.R, run_d0_sweep.R,
  run_scale_adversarial.R, make_figures.R, ...}, validated_results/{*.csv, figures/*.png}
- `pipeline/phase3_write/`: manuscript.tex, references.bib (21), figures/*.pdf (4),
  tab_*.tex (5), data/*.csv (9), briefings/{formula_code_audit, anomaly_log,
  literature_briefing, modeling_briefing}.md, decision_log.md
- `pipeline/phase4_polish/`: round_0/paper_grade.md, round_1/{manuscript.tex (final),
  manuscript.pdf, paper_grade.md, tab_*.tex, references.bib, figures/}
- `pipeline/`: pipeline_log.md, pipeline_state.md, final_report.md
