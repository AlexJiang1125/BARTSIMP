# Formula ↔ Code Audit (paper-modeler-owned, canonical)

Every load-bearing equation in `methodology_specification.md` checked against the
validated implementation (`dev/bartsimp_pg.R`, `dev/rbart.R`, and the Phase-2
`validated_code/`). Status MATCH means the LaTeX in the "Reconciled LaTeX" column
is the version the Writer must copy verbatim. **Zero unresolved MISMATCH rows.**

| # | Quantity | Spec | Code location | Status | Reconciled LaTeX |
|---|----------|------|---------------|--------|------------------|
| F1 | Outcome model | `k_i ~ Binomial(n_i,p_i)`, `logit p_i = η_i = α₀ + f(x_i) + z(s_i)[+ξ_i]` | bartsimp_pg.R:631 `eta = binary_offset + f + z + xi` | MATCH | `k_i \mid p_i \sim \mathrm{Binomial}(n_i,p_i),\quad \operatorname{logit} p_i = \eta_i = \alpha_0 + f(x_i) + z(s_i) + \xi_i` |
| F2 | Intercept/offset | `α₀ = binary_offset = logit(p̄)`, `p̄=Σk_i/Σn_i` | bartsimp_pg.R:459 `binary_offset <- qlogis(p_bar)` | MATCH | `\alpha_0 = \operatorname{logit}\bar p,\qquad \bar p = \sum_i k_i / \sum_i n_i` |
| F3 | BART mean | `f(x)=Σ_{t=1}^m g(x;T_t,M_t)` | rbart.R:457 `ensemble_predict = Reduce(+, tree_predict)` | MATCH | `f(x)=\sum_{t=1}^{m} g(x;T_t,M_t)` |
| F4 | Leaf prior | `β_{t,ℓ}\midT ~ N(0,τ²)`, `τ²=(3/(k_BART√m))²` | bartsimp_pg.R:448 `tau2 <- (3/(k_bart*sqrt(m_trees)))^2` | MATCH | `\beta_{t,\ell}\mid T_t \sim N(0,\tau^2),\quad \tau^2=\bigl(3/(k_{\mathrm{BART}}\sqrt m)\bigr)^2` |
| F5 | Chipman branch prior | `P(split at depth d)=α_tree(1+d)^{−β_tree}` | rbart.R:12 / log_prior_node | MATCH | `\mathbb P(\text{split at depth }d)=\alpha_{\mathrm{tree}}(1+d)^{-\beta_{\mathrm{tree}}}` |
| F6 | SPDE field | `z(s)=Σ_j A(s)_j u_j`, `u~N(0,Q_ψ^{-1})` | bartsimp_pg.R GP step (Woodbury, A,Qu) | MATCH | `z(s)=\sum_j A(s)_j u_j,\quad u\sim N(0,Q_\psi^{-1})` |
| F7 | PG augmentation | `ω_i~PG(n_i,η_i)`; `ȳ_i=(k_i−n_i/2)/ω_i`; `ȳ_i\midω_i,η_i~N(η_i,1/ω_i)` | bartsimp_pg.R:331-332,635 `rpg(h=w*n_i,...)`, `kappa=w*(k_i-n_i/2)` | MATCH | `\omega_i\sim\mathrm{PG}(n_i,\eta_i),\quad \bar y_i=(k_i-n_i/2)/\omega_i,\quad \bar y_i\mid\omega_i,\eta_i\sim N(\eta_i,1/\omega_i)` |
| F8 | DEFF tempering | `w_i=1/[1+(n_i−1)ρ_DEFF]`; replace `(n_i,k_i)→(w_i n_i, w_i(k_i−n_i/2))` | bartsimp_pg.R:635-636 `w_deff` lines (collapse at ρ=0) | MATCH | `w_i = 1/[1+(n_i-1)\rho_{\mathrm{DEFF}}]` |
| F9 | Depth-floor SD1 (complete init to d₀) | trees init complete to depth d₀ | methods_bartsimp.R d0 passed; rbart ensemble init | MATCH | `\gamma_{\ell k}=1\ \forall k,\ \ell\le d_0` (complete coarse layers) |
| F10 | Depth-floor SD2/SD3 (prune guard) | reject prune of NOG at depth ≤ d₀; depth>d₀ moves unchanged | rbart.R:331-332,408-414 `nogs[depth>=d0]`; reverse-NOG count floor-filtered (272,306) | MATCH | reject prune at `d_{\mathrm{NOG}} < d_0`; moves at depth `>d_0` unchanged; `d_0=0` recovers legacy |
| F11 | Depth-floor d₀ rule | `d₀=⌈√(log N)⌉` | run_paired_dgp5.R:53 `ceiling(sqrt(log(N)))` | MATCH | `d_0 = \lceil \sqrt{\log N}\,\rceil` |
| F12 | Area functional (estimand) | `χ^a=Σ_{i∈a} w_i p_i`, `w_i=n_i/Σ_{j∈a}n_j` | coverage.R `compute_admin_cov` (in-sample n_i weighting) | MATCH | `\chi^a=\sum_{i\in a} w_i\, p_i,\quad w_i = n_i/\textstyle\sum_{j\in a} n_j` |
| F13 | Calibration ratio | `ratio = RMSE/post_sd` | run_paired_dgp5.R:90 `admin_ratio = rmse/psd` | MATCH | `\text{ratio}=\mathrm{RMSE}/\overline{\mathrm{post\_sd}}` |
| F14 | Standardized bias from ratio | `δ=√(ratio²−1)` (ratio≈1.08⇒δ≈0.41) | diagnostic identity (report) | MATCH | `\delta = \sqrt{\mathrm{ratio}^2 - 1}` |
| F15 | In-sample efficient variance (theory) | `V^a_smp = σ²/π_a` | theory draft (Gaussian sibling) | MATCH (theory-only) | `V^a_{\mathrm{smp}} = \sigma^2/\pi_a` |
| F16 | MC standard error | `SE = sd/√B` | run_paired_dgp5.R:151 `mc_se <- sd/sqrt(n)` | MATCH | `\mathrm{SE} = \mathrm{sd}/\sqrt{B}` |

**Resolution note.** Two spec equations needed a regime/label fix, neither a
code bug: (i) Theorem 1 must cite KvdVvZ **Thm 5.4** (nominal) not 5.3
(conservative) — a citation/regime correction internal to the theory prose,
carried into the manuscript with the corrected `β⋆ > ν/2 + 1/4` form. (ii) The
theory variance is the **in-sample** `V^a_smp = σ²/π_a` (the posterior conditions
on the design), not the superpopulation `V_0^a`; the manuscript states the
in-sample target as primary and the superpopulation add-on as a remark. Both are
exposition fixes on the theory side, not implementation MISMATCHes.
