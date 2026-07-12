# Anomaly Log (paper-modeler-owned, canonical)

Anomalies detected in the validated results. Types: 1 = result contradicts a
claim; 2 = unexpected magnitude/sign; 3 = data-quality/harness artifact;
4 = benign/expected. Every Type 1–3 anomaly MUST be addressed in Discussion by ID.

| ID | Type | What | Where | Why / mechanism | When it bites | Meaning for the paper |
|----|------|------|-------|-----------------|---------------|-----------------------|
| A1 | 2 | Primary paired Δcov80 = +0.0145 ± 0.0102 is below 2 MC SE at 24 reps (one-sided p≈0.08); below the spec's +0.02 bar | `dgp5_paired_primary_deltas.csv` | 24 reps insufficient to resolve a genuinely modest effect; the per-rep cov80 is a noisy 0/1-average over admin units | small replicate count | Honest CONDITIONAL-GO limitation; the d₀ sweep (clears +0.02) and the scale/adversarial consistency carry the claim, not the single primary contrast. Stated in Limitations. |
| A2 | 3 | Pixel cov80 ≈ 0.866–0.873 (both arms), below the spec's ≥0.90 conservative band | all dgp5 CSVs (`pixel_cov80`) | the reduced-scale validation mesh/N lowers the absolute pixel level; it is a property of the validation DGP mesh, NOT the depth-floor | only in the light validation harness; production `sim_mbg_surface` reports pixel over-coverage ~0.97 | The load-bearing fact is the knob does NOT change it (0.869→0.866, Δ=−0.003). "No degradation of pixel calibration" holds. Discuss as mesh-driven, not a regression. |
| A3 | 4 | Adversarial: at σ_gp=1.0,1.5 the depth-floor moves cov80 *down* (0.820→0.812, 0.821→0.804) | `dgp5_adversarial_aggregated.csv` | where the GP already supplies the coarse signal, legacy sits at/above 0.80; the floor recenters *toward* 0.80 from above rather than inflating | GP-dominated regime | Expected, desirable: recentering not monotone inflation. Confirms the mechanism (the floor is inert/gently-correcting where the bias is absent). |
| A4 | 4 | Cluster cov80 ≈ 0.9997–1.0 for both BARTSIMP arms (forced over-coverage) | `dgp5_paired_primary_aggregated.csv` | binomial discreteness floor at small n_i (Theorem 9.X / cluster-floor); a feature, not a defect | small n_i clusters | The fine scale must NOT be "fixed"; the floor leaves it intact. Motivates the "calibrated jointly, not uniformly nominal" framing. |
| A5 | 4 | A=40 advantage is small (+0.006) vs A=20 (+0.034) | `dgp5_scale_aggregated.csv` | at large A legacy already sits near 0.80; less bias to recenter | coarse aggregation | Advantage SHRINKS where unnecessary — consistent with a targeted recentering knob, not a universal inflator. |
| A6 | 4 | RF-GLS / SPARforest carried at low reps (roster-context only) | `comparator_inventory.md` | ~480 s/rep makes full paired-rep budget infeasible | always | Reported as context rows with reduced-rep caveat; does not affect the d₀ gate (the paired arm-to-beat is BARTSIMP-PG d₀=0, run full-rep). |

No Type 1 anomalies (no result contradicts a stated claim). A1–A2 are the
honesty-critical items the Discussion and Limitations must address by ID.
