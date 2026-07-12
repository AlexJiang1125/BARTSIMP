# BARTSIMP-PG sampler — developer / performance notes

Runtime engineering log for the pure-R Pólya-Gamma sampler
(`dev/rbart.R` + `dev/bartsimp_pg.R`). Records the time audit, the current
bottlenecks, what optimizations were tried (and how each was verified), and —
most importantly — **how to roll back** if a change misbehaves.

Keep this current whenever sampler *performance* behavior changes. It is the
companion to `notes/workflow.tex` (which tracks *modeling* changes); pure
speed work that does not alter results lives here.

---

## 1. File topology — read before editing or rolling back

The sampler exists as **two on-disk copies** that are meant to stay in sync.
Neither is git-tracked as code (this notes file *is* git-tracked).

| role | path |
|---|---|
| **cluster canonical** (shipped to HPC by `scripts/sync_to_cluster.sh`) | `…/Claude/Projects/BARTSIMP_glm/bartsimp_cluster/dev/{rbart.R, bartsimp_pg.R}` |
| **local twin** (sourced by local `model/04b_predict_aggregate.R`) | `…/GitHub/BARTSIMP/dev/{rbart.R, bartsimp_pg.R}` |

**Current sync state (2026-06-13):**
- `rbart.R` — **byte-identical** in both copies. ✓
- `bartsimp_pg.R` — **byte-identical** in both copies. ✓ The earlier drift
  (cluster copy carried a "Route 3 DEFF-tempering" `deff_rho` feature the local
  twin lacked) was reconciled this session by backing up both copies, copying
  cluster (canonical) → local, and confirming `diff -q`. The `chol_perm` flag
  (§4d) was then landed in the cluster copy and re-synced, so both now carry it.

**Rule when editing the sampler:** edit the **cluster canonical** copy, then
reconcile cluster → local (don't hand-edit the stale twin), then
`bash scripts/sync_to_cluster.sh` to push to the HPC. Confirm with
`diff -q <local> <cluster>` — expect **identical** for both files.

**Always use the CRAN framework R**, not the conda/miniforge R that is first
on `PATH` (the latter segfaults on INLA/Rcpp from an ABI mismatch):

```
/Library/Frameworks/R.framework/Resources/bin/Rscript <script>.R
```

---

## 2. Runtime audit

Config notation: `n_cl` = #clusters, `m_trees` = #BART trees, `n_mesh` = #SPDE
mesh nodes. Production fit (`model/03_fit_bartsimp_pg.R`): `n_iter=5000`,
`burn=2000`, `m_trees=100`, `bd_mode="conditional"`, `chol_perm=TRUE`
(fill-reducing Cholesky; default since 2026-06-13 — see §4d), no in-sampler
prediction.

### True per-iteration wall time (no Rprof)

| config | n_mesh | `chol_perm=F` s/iter | `chol_perm=T` s/iter | speedup |
|---|---|---|---|---|
| n_cl=600, m_trees=50 (audit baseline) | 672 | 0.077–0.085 | **0.057** | ~1.35× |
| n_cl=1380, m_trees=100 (production scale) | 730 | 0.146–0.152 | **0.114–0.117** | **~1.30×** |

(Production-scale numbers above are post-4e. Pre-4e they were `chol_perm=F`
≈ 0.156–0.159, `chol_perm=T` ≈ 0.123–0.126.)

Pre-optimization the baseline config was ≈ **0.27 s/iter**. The landed tree
work (4a+4b) gave ~3.2× at `chol_perm=FALSE`; `chol_perm=TRUE` (4d) gives a
further ~1.3× on top at production scale, and the `tree_leaf_assign`
allocation tuning (4e) another ~1.08× — **~4.5× cumulative**. These
per-component multipliers were measured on **short runs** (n_iter≈40–60,
immature BART trees) via `total/n_iter`; the end-to-end *steady-state*
speedup is substantially larger and scale-dependent — see the next subsection.

### Old-vs-new end-to-end benchmark (steady-state, two-point method)

Full original sampler (scalar `rbart` `.bak_20260613_150257` + pre-`chol_perm`
sampler `.bak_20260613_154712`) vs current code, three arms interleaved,
**two-point per-iter** `(t140 − t40)/100` so fixed mesh-build/setup cancels and
the slope reflects the mature-tree marginal cost (what a 5000-iter run actually
pays). `bd_mode="conditional"`, same seed/data across arms.

| scale | arm | steady-state s/iter | speedup vs OLD |
|---|---|---|---|
| n_cl=600, m_trees=100 | OLD (scalar) | 0.760 | — |
| | NEW `chol_perm=F` | 0.103 | **7.4×** |
| | NEW `chol_perm=T` | 0.079 | **9.6×** |
| n_cl=1380, m_trees=100 | OLD (scalar) | 1.660 | — |
| | NEW `chol_perm=F` | 0.166 | **10.0×** |
| | NEW `chol_perm=T` | 0.138 | **12.0×** |

- **Bit-identical default — confirmed at BOTH scales.** NEW `chol_perm=FALSE`
  reproduces OLD draws exactly: `max|new−old| = 0` across f/z/sigma_m/rho/u at
  the same seed. (NEW-T differs only by the seeded GP-draw realization — same
  distribution, different draw.)
- **Speedup grows with run length.** OLD's per-iter *rises* as BART trees
  mature (n_cl=1380: 1.04 s/iter n40-avg-incl-setup → 1.66 s/iter 140-iter
  marginal) because scalar descent cost scales with tree size; NEW's vectorized
  descent stays ~flat. Tree *shapes* are identical between arms (the zero-diff
  proves it) — purely a descent-speed gap. The short-run `total/n_iter` table
  understates the steady-state win; the 4.5× above is a **floor**.
- **Speedup grows with `n_cl` (counter to first intuition).** The SPDE mesh is
  nearly `n_cl`-invariant (`n_mesh` 672→730 from 600→1380 clusters), so the
  GP/Cholesky cost is roughly fixed while tree-descent cost scales with the
  number of rows. At production scale trees dominate more, so the
  tree-vectorization win matters *more*, not less.
- **Production projection (5000 iters, steady-state, n_cl=1380):** OLD ≈
  **138 min** (~2.3 h) → NEW `chol_perm=T` ≈ **11.5 min**. `chol_perm`'s
  marginal share shrinks at steady state (1.20× here vs the 1.30× short-run
  table) because mature trees crowd out the Cholesky.

Bench harness (throwaway, in `/tmp`): `bench_arm.R` (one arm per process to
avoid version name-collisions), interleaved bash driver, `bench_summarize.R`
(n_cl=600) / `bench_sum1380.R` (n_cl=1380). Method is fully described above if
the scripts are gone.

> ⚠️ **Rprof pitfall.** At a fine sampling interval (e.g. 0.005s) with the
> old scalar tree loops, Rprof inflated wall time ~40× (it thrashes on the
> many small allocations). **Rprof *proportions* are valid; absolute Rprof
> time is not.** Always cross-check absolute timing with a no-Rprof
> `Sys.time()` run. The audit below uses `interval=0.01`.

### Rprof allocation breakdown (n_cl=600, m_trees=50, post tree-optimization)

By total %:

| component | total % | notes |
|---|---|---|
| `margprob_woodbury` (Woodbury GP marginal, ×2/iter) | **43.9%** | now the #1 cost |
| &nbsp;&nbsp;└ `Matrix::Cholesky` | 31.1% | sparse factorization of `Qu = AᵀΩA + Q_psi` |
| &nbsp;&nbsp;└ `determinant` | 9.2% | read off the cached factor |
| `ensemble_sweep` (all BART tree machinery) | **32.5%** | was ~74% pre-opt |
| &nbsp;&nbsp;└ `tree_leaf_assign` | 24.7% (9.9% self) | was ~74% self pre-opt |
| &nbsp;&nbsp;└ `tree_draw_mu` | 11.6% | |
| &nbsp;&nbsp;└ `tree_bd` (grow + prune) | 10.9% | |
| `rpg` (Pólya-Gamma draws, BayesLogit) | **18.2%** | `.C` compiled |
| `make_Qpsi` / `inla.spde.precision` | 4.0% | |

By self %: `.Call` 41% (CHOLMOD/Matrix + INLA), `.C` 18% (`rpg`),
`tree_leaf_assign` 9.9% + its helpers `which`/`ifelse`/`cbind` (~14% combined).

### Rprof allocation at production scale (n_cl=1380, m_trees=100): F vs T

Total %; the key story is **Cholesky's share collapsing** once the
fill-reducing reorder is on (`/tmp/bartsimp_pg_audit3.R`):

| component | `chol_perm=F` total % | `chol_perm=T` total % |
|---|---|---|
| `ensemble_sweep` (all BART tree machinery) | 49.9% | **57.8%** |
| &nbsp;&nbsp;└ `tree_leaf_assign` | 39.2% (16.3% self) | 44.8% (21.1% self) |
| &nbsp;&nbsp;└ `tree_draw_mu` | 17.8% | 20.4% |
| &nbsp;&nbsp;└ `tree_bd` (grow + prune) | 18.1% | 21.7% |
| `rpg` (Pólya-Gamma, BayesLogit `.C`) | 22.0% | **27.4%** |
| `margprob_woodbury` (×2/iter) | 25.3% | 10.3% |
| &nbsp;&nbsp;└ `Matrix::Cholesky` | **17.4%** | **1.7%** |
| &nbsp;&nbsp;└ `determinant` | 6.3% | 7.1% |

(At `m_trees=100` the trees take a bigger share than the `m_trees=50` baseline
above, so Cholesky starts at 17.4% rather than 31.1%; the absolute per-call
28× win is the same.)

---

## 3. Current bottleneck picture

With `chol_perm=TRUE` (landed), the per-iteration cost is decisively
**BART trees > Pólya-Gamma draws ≫ linear algebra**:
- `ensemble_sweep` ~50–58%, dominated by `tree_leaf_assign` (~42% total after
  4e) — still the largest *pure-R* residual.
- `rpg` ~27% — compiled `.C`, effectively at its floor.
- `margprob_woodbury` down to ~10%, now mostly `determinant` (~7%), not the
  Cholesky (1.7%). **The Cholesky is no longer a bottleneck.**

With `chol_perm=FALSE` (default), the picture is the pre-4d one: Cholesky is
~17% (m_trees=100) to ~31% (m_trees=50) of the iter and is the #2 cost behind
trees.

`tree_leaf_assign` is now allocation-tuned (4e). The remaining levers on it are
**structural** — cut the ~3 calls/tree/sweep by reusing the leaf assignment
across `tree_bd`→`tree_draw_mu` on rejected moves (fiddly, needs 4b-style
bit-verification) — or **compiled** (Rcpp descent), which cuts against the
deliberately pure-R design. Both are bigger projects than 4e for a similar or
smaller marginal return; not currently pursued.

---

## 4. Optimizations tried & verified

### 4a. LANDED — `tree_leaf_assign` vectorization (`dev/rbart.R`)

Replaced the per-row scalar tree walk with a vectorized level-synchronous
descent (strict `<` at each split, requires `X` to be a matrix).
- **Effect:** ~12× faster in isolation; dominant contributor to the ~3.2×
  end-to-end speedup. `tree_leaf_assign` self-time 73.6% → 9.9%.
- **Verified bit-identical:** 400 BD-grown trees, 0 mismatches vs the original
  scalar walk (`/tmp/post_edit_verify.R`).

### 4b. LANDED — `tree_draw_mu` cleanup (`dev/rbart.R`)

(1) Vectorized the per-leaf μ Gibbs draw with `split()` over leaf assignment;
(2) returns the leaf assignment so `ensemble_sweep`/`ensemble_sweep_marginal`
build the new fit by `mu[leaf]` instead of re-walking the tree
(`tree_predict`) — removing one redundant tree walk per tree per sweep.
- **RNG invariant preserved:** `rnorm(L, mean_vec, sd_vec)` consumes the
  stream in the same order as `L` scalar `rnorm(1, …)` calls, so seeded output
  is unchanged.
- **Verified bit-identical:** vs the pre-edit backup at seed 123,
  `max|new − old| = 0` for `f_draws`, `z_draws`, `sigma_m_draws`,
  `rho_draws`, `u_draws` (`/tmp/verify_drawmu.R`).

### 4c. REJECTED — caching the symbolic Cholesky factor (`Matrix::update`)

Idea: `margprob_woodbury` factorizes `Qu = AᵀΩA + Q_psi` twice per iteration
(current + MH proposal), and `Qu`'s sparsity pattern is invariant across all
calls (A fixed; only `Q_psi`/`Ω` *values* change). So cache one symbolic
factor and `Matrix::update()` it numerically instead of a fresh
`Matrix::Cholesky()`.
- **Bit-identical:** yes — solves and log-determinant match to 0.000e+00.
- **Speed:** **1.01× — useless.** On the real production `Qu` (n_mesh=730),
  fresh `Cholesky` = 14.5 ms vs `update()` = 14.3 ms. With `perm=FALSE` the
  cost is the *numeric* factorization; the *symbolic* analysis that `update()`
  skips is ~1% of it. (Probe: `/tmp/chol_update_probe2.R`.)
- **Conclusion:** not worth implementing.

### 4d. LANDED (opt-in, default off) — `perm=TRUE` fill-reducing Cholesky

Investigating 4c revealed the real issue: `margprob_woodbury` uses
`Matrix::Cholesky(Qu, perm=FALSE, super=FALSE)`. The natural mesh ordering
gives catastrophic fill-in.

| variant | per-call time | nnz(L) |
|---|---|---|
| `perm=F super=F` (**current**) | 14.1–14.5 ms | 191,719 |
| `perm=T super=F` | 0.70 ms | 31,338 |
| `perm=T super=NA` (auto) | **0.45–0.51 ms** | 31,338 |

A fill-reducing permutation makes the factorization **~28–32× faster** (6× less
fill). Cholesky is ~70% of `margprob_woodbury`, so this cuts that function ~3×;
estimated **~1.2–1.4× end-to-end** (larger at small-n_cl configs where
Cholesky is a bigger share). Probes: `/tmp/chol_perm_probe.R`,
`/tmp/perm_compare.R`.

**Evaluation — fill-reducing is NOT "coarser", it is exact** (perm_compare.R,
real production Qu, n_mesh=730):

- Posterior mean `m_u`: `max|permF − permT| = 1.0e-15`; `log|Qu|`: `9.1e-13`.
  ⇒ the MH acceptance math is identical to machine precision (same mixing).
- GP draw covariance vs the exact `Qu⁻¹` (built the draw operator explicitly,
  no Monte-Carlo): perm=F `1.7e-15`, **perm=T + Pt-fix `7.8e-16`** — both exact.
  So perm=T loses **no** accuracy.
- The only real difference: under a fixed seed the two correct draws are
  **different realizations** (`max|u_permF − u_permT| ≈ 3.0`) of the **same**
  distribution.

**Correctness gotcha (must-fix).** The draw `u = m_u + L⁻ᵀη` via
`solve(Fchol, η, system="Lt")` is only valid when `perm=FALSE`. With `perm=TRUE`
the factor is of `P Qu Pᵀ`, so the draw must also apply `system="Pt"`. The naive
flip *without* the Pt step gives `max|Cov − Qu⁻¹| = 0.61` — a **silent
wrong-covariance bug**. The unified correct form
`m_u + solve(F, solve(F, η, "Lt"), "Pt")` is bit-identical to the current code
when `perm=FALSE` (Pt is the identity permutation → no-op), so it is safe as a
single code path.

**LANDED 2026-06-13 — `chol_perm` flag, default FALSE in the sampler (opt-in, no forced re-run).**

> **Production wiring (2026-06-13):** `model/03_fit_bartsimp_pg.R` exposes a
> `--chol_perm` optparse option **defaulting to `TRUE`**, so production fits take
> the fill-reducing path out of the box. Pass `--chol_perm FALSE` to bit-reproduce
> frozen `chol_perm=FALSE` results. The chosen value is echoed in the run's startup
> log line. The sampler default itself stays `FALSE` so direct `run_sampler_*` calls
> (e.g. `model/04b`, benchmarks) are unchanged.

Implementation (6 edits in `bartsimp_pg.R`, both copies byte-identical):
- `margprob_woodbury(yhat, Omega, A, Q_psi, perm = FALSE)` → `Matrix::Cholesky(…, perm = perm)`.
- `run_sampler_bartsimp_pg(…, chol_perm = FALSE)` threads `perm = chol_perm`
  into both `margprob_woodbury` call sites (current + MH proposal); only the
  `bd_mode="conditional"` path uses Woodbury, so only it is affected.
- `draw_z_woodbury` and `draw_z_woodbury_rsr` rewritten to the unified
  Pt-correct two-solve form (`w <- solve(F, η, "Lt"); m_u + solve(F, w, "Pt")`),
  with the `rnorm()` call left in place so the RNG stream is unchanged.
- `chol_perm` echoed into the output list for provenance.
- Marginal-path `Qu_chol` in `rbart.R` left at `perm=FALSE` (not the
  conditional hot path) — unchanged.

Verification (`/tmp/verify_chol_perm.R`, framework R, seed=123):
- **`chol_perm=FALSE` is byte-identical** to the pre-perm backup
  (`bartsimp_pg.R.bak_20260613_154712`): `max|curF−bak| = 0.000e+00` for
  `f`, `z`, `sigma_m`, `rho`, `u` draws. The Pt-no-op refactor changed nothing.
- **`chol_perm=TRUE` is the same distribution, different realization**: draws
  genuinely differ (so perm=T is truly engaged, not a silent fallback), but
  per-cluster posterior-mean prevalence correlates `0.994` with the FALSE run,
  and summaries sit within MC noise of a 60-draw chain. Distributional
  correctness is guaranteed by the exact-covariance proof above (`7.8e-16`).

Post-fix audit (`/tmp/bartsimp_pg_audit3.R`, production scale): **1.30×
end-to-end** at n_cl=1380/m_trees=100 (0.159→0.123 s/iter); Cholesky's Rprof
share **17.4%→1.7%** (see §2 table). Bottleneck now firmly on the BART trees.

Switching the *default* to TRUE later would require refreshing frozen seeded
references and re-running Sim 7 (conclusions unchanged within MC error).

### 4e. LANDED — `tree_leaf_assign` allocation tuning (`dev/rbart.R`)

Once 4d made Cholesky cheap, `tree_leaf_assign` was the largest single cost
(~45–50% total, the dominant *pure-R* residual). 4a had vectorized the descent;
4e tunes its inner loop, which ran ~3× per tree per sweep (`tree_predict` +
`tree_bd` + `tree_draw_mu`), i.e. ~300 calls/sweep at `m_trees=100`. Two
allocation-level changes, both **bit-identical**:
- `X[cbind(idx, sv)]` → column-major linear index `X[idx + (sv-1L)*n]`: drops a
  `|idx|×2` index-matrix allocation at every level of every descent.
- `ifelse(goleft, left[node], right[node])` → masked assignment
  (`nxt <- left[node]; nxt[!goleft] <- right[node][!goleft]`): avoids `ifelse`'s
  both-branch evaluation + NA bookkeeping (safe — split `X` is NA-free).

- **Verified bit-identical (3 ways, `/tmp/verify_leafassign_microopt.R`):**
  (i) kernel unit test over 300 random tree shapes — **0 mismatches** vs a
  brute-force per-row scalar walk *and* vs the old `cbind`/`ifelse` kernel;
  (ii) end-to-end sampler draws at seed 123 — `max|new−old| = 0.000e+00` for
  `f`/`z`/`sigma_m`/`rho`/`u` vs the pre-edit backup; (iii) same `0.000e+00`
  vs the **original scalar** rbart (`.bak_20260613_150257`), so the full
  4a+4b+4e chain is faithful to the pre-optimization implementation.
- **Speed (interleaved A/B, `/tmp/time_microopt_ab.R`, production scale):**
  **~1.07–1.09×** per-iter at `chol_perm=TRUE` (0.126→0.117 s/iter), ~1.04× at
  FALSE. Modest — 4a already captured the 12× win; this is the allocation
  residual. `tree_leaf_assign` total Rprof share ~50%→42%.
  > Note: its Rprof *self*% appears to *rise* (23%→38%) — an artifact of
  > attribution, not a regression. Work formerly in child frames `cbind`/
  > `ifelse`/`which` is now inline base-R `[`/`[<-`/Ops counted as self. The
  > honest metric is total% (fell) and wall-clock (fell).

---

## 5. Rollback procedures

### Revert the LANDED tree optimizations (4a + 4b)

A single pre-edit backup reverts **both** (4a and 4b were made after the same
backup). Backup = the original 27,698-byte `rbart.R`.

```sh
# in the GitHub repo:
cp dev/rbart.R.bak_20260613_150257 dev/rbart.R
# in the cluster project:
cp <…>/bartsimp_cluster/dev/rbart.R.bak_20260613_150257 \
   <…>/bartsimp_cluster/dev/rbart.R
# push to HPC:
bash <…>/bartsimp_cluster/scripts/sync_to_cluster.sh
# confirm both copies match again:
diff -q dev/rbart.R <…>/bartsimp_cluster/dev/rbart.R   # -> identical
```

Backups (identical content, 27,698 bytes) live at:
- `…/GitHub/BARTSIMP/dev/rbart.R.bak_20260613_150257`
- `…/bartsimp_cluster/dev/rbart.R.bak_20260613_150257`

Restoring this backup reverts **4a + 4b + 4e** at once (back to the scalar
walk), since they all live in `rbart.R`.

### Revert only the LANDED `tree_leaf_assign` micro-opt (4e)

To drop just 4e and keep 4a+4b, restore the pre-4e backup
`rbart.R.bak_microopt_20260613_160821` (28,384 bytes, identical in both trees)
instead, then re-sync the same way. Because 4e is bit-identical, this changes
no results — only the inner-loop allocation pattern.

### Revert the LANDED `chol_perm` flag (4d)

The `chol_perm` work was made after the pre-perm backup
`bartsimp_pg.R.bak_20260613_154712`. Note: **the default (`chol_perm=FALSE`)
is byte-identical to that backup**, so a rollback is only needed to remove the
flag entirely, not to restore default behavior.

> ⚠️ **Restore from the CLUSTER backup only, then re-sync.** The two
> `bak_20260613_154712` files are **not** the same: the cluster one (38,969 B,
> has `deff_rho`) is the canonical pre-perm state; the local one (36,356 B,
> *no* `deff_rho`) is the pre-sync **drifted** twin. Restoring each copy from
> its own `.bak` would reintroduce the drift. Always restore the cluster
> canonical, then copy cluster → local.

```sh
# restore the cluster canonical, then propagate to the local twin:
cp <…>/bartsimp_cluster/dev/bartsimp_pg.R.bak_20260613_154712 \
   <…>/bartsimp_cluster/dev/bartsimp_pg.R
cp <…>/bartsimp_cluster/dev/bartsimp_pg.R  dev/bartsimp_pg.R   # cluster -> local
bash <…>/bartsimp_cluster/scripts/sync_to_cluster.sh
diff -q dev/bartsimp_pg.R <…>/bartsimp_cluster/dev/bartsimp_pg.R   # -> identical
```

The change spans `margprob_woodbury`, `draw_z_woodbury` (+ `_rsr`), the
signature/output-list, and both call sites — all within `bartsimp_pg.R` (the
marginal-path `Qu_chol` in `rbart.R` was left untouched). Any results produced
under `chol_perm=TRUE` are a different (valid) chain and must be regenerated if
you need to reproduce them bit-for-bit.

---

## 6. How to re-verify (recipe)

Use the framework R (§1). Pattern used by all verification scripts in `/tmp`:

1. **Bit-identity vs backup:** source the edited file, run the sampler at a
   fixed seed → `fit_new`; source the `.bak_*` original over the same
   functions, run again → `fit_old`; assert
   `max|fit_new[[d]] − fit_old[[d]]| == 0` for
   `d ∈ {f,z,sigma_m,rho,u}_draws`. (Templates: `/tmp/verify_drawmu.R` for
   rbart swaps; `/tmp/verify_chol_perm.R` for the `chol_perm` default path
   — the latter also checks `chol_perm=TRUE` matches summaries + prevalence.)
2. **Component bit-identity** (for linear-algebra swaps like 4c/4d): compare
   `solve(system="A")`, `solve(system="Lt")`, and `determinant(…, sqrt=FALSE)`
   on the real `Qu`; for `perm=TRUE` also check the draw operator covariance
   `M Mᵀ` against `Qu⁻¹`. (Templates: `/tmp/chol_update_probe2.R`,
   `/tmp/perm_compare.R`.)
3. **Timing:** always a no-Rprof `Sys.time()` per-iter run; use Rprof only for
   *proportions* (`interval=0.01`). (Template: `/tmp/bartsimp_pg_audit3.R` —
   compares both `chol_perm` settings at production scale.)
