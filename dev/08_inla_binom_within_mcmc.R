# =============================================================================
# Validation V4 -- "INLA-binom-within-MCMC" timing comparator
# =============================================================================
# A minimal alternative to the PG-Woodbury sampler: instead of exploiting
# Pólya-Gamma augmentation to recover Gaussian conjugacy and use the Woodbury
# fast path, we call inla(family="binomial") directly at every marginal-
# likelihood evaluation. This is the "naïve" approach Stage 2 of the original
# design was meant to avoid.
#
# Per sweep we run:
#   - For each of m trees: 2 inla() calls (current state + BD proposal)
#       leaves enter as f(leaf_id, model="iid", fixed prec = 1/tau^2)
#       spatial enters as f(spatial_idx, model = spde)
#       other trees' contributions enter as offset()
#   - For psi MH: 2 inla() calls (current + proposed psi)
# Total inla() calls per sweep: 2*(m + 1)
#
# In SMOKE_TEST mode we run a few iterations to verify the pipeline and
# obtain a per-iteration timing estimate. Production-scale runs are intended
# for a cluster (estimated 10+ hours at m=20, 1500 iter, n=200).
# =============================================================================

source("dev/bartsimp_pg.R")

SMOKE_TEST <- TRUE

if (SMOKE_TEST) {
  N_CL      <- 100
  N_ITER    <- 5
  M_TREES   <- 3
  MESH_EDGE <- c(0.15, 0.40)
  MESH_CUT  <- 0.06
} else {
  N_CL      <- 200
  N_ITER    <- 1500
  M_TREES   <- 20
  MESH_EDGE <- c(0.08, 0.25)
  MESH_CUT  <- 0.03
}

cat(sprintf("Config: n_cl=%d, n_iter=%d, m_trees=%d\n", N_CL, N_ITER, M_TREES))
cat(sprintf("Expected inla() calls per sweep: %d\n", 2 * (M_TREES + 1)))

# -----------------------------------------------------------------------------
# Data + mesh + SPDE
# -----------------------------------------------------------------------------
dat <- gen_data(seed = 123, n_cl = N_CL,
                covariate_signal = TRUE, spatial_signal = TRUE)
mesh <- INLA::inla.mesh.2d(loc = dat$locs, max.edge = MESH_EDGE, cutoff = MESH_CUT)
spde_default <- INLA::inla.spde2.pcmatern(mesh,
                                          prior.range = c(0.1, 0.5),
                                          prior.sigma = c(1.0, 0.5))
A_obs <- INLA::inla.spde.make.A(mesh = mesh, loc = dat$locs)
cat(sprintf("Mesh nodes: %d\n", mesh$n))

# BART prior leaf precision
tau  <- 3 / (2 * sqrt(M_TREES))
tau2 <- tau^2

# -----------------------------------------------------------------------------
# Helper: run inla(family="binomial") for the tree-j conditional model.
# Returns marginal log-likelihood and leaf posterior means.
# -----------------------------------------------------------------------------
fit_binom_for_tree <- function(leaf_id, offset_per_obs, spde_obj = spde_default) {
  # Ensure leaf_id is dense integer starting from 1
  leaf_id <- as.integer(factor(leaf_id))
  stk <- INLA::inla.stack(
    data = list(y = dat$k_i, Ntrials = dat$n_i),
    A    = list(A_obs, 1),
    effects = list(
      spatial_idx = seq_len(spde_obj$n.spde),
      data.frame(leaf_idx = leaf_id, offset_val = offset_per_obs)
    )
  )
  res <- INLA::inla(
    y ~ -1 + offset(offset_val) +
        f(leaf_idx, model = "iid",
          hyper = list(prec = list(initial = log(1/tau2), fixed = TRUE))) +
        f(spatial_idx, model = spde_obj),
    family  = "binomial",
    Ntrials = INLA::inla.stack.data(stk)$Ntrials,
    data    = INLA::inla.stack.data(stk),
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = FALSE),
    control.compute   = list(dic = FALSE, waic = FALSE),
    verbose = FALSE
  )
  leaf_means <- if (!is.null(res$summary.random$leaf_idx))
                  res$summary.random$leaf_idx[, "mean"]
                else numeric(0)
  list(mlik = as.numeric(res$mlik[1]), leaf_means = leaf_means)
}

# -----------------------------------------------------------------------------
# Helper: marginal lik for the psi MH step.
# offset = current sum of trees at each cluster; spde_obj has the proposed psi.
# -----------------------------------------------------------------------------
fit_binom_for_psi <- function(offset_per_obs, spde_obj) {
  stk <- INLA::inla.stack(
    data = list(y = dat$k_i, Ntrials = dat$n_i),
    A    = list(A_obs, 1),
    effects = list(
      spatial_idx = seq_len(spde_obj$n.spde),
      data.frame(offset_val = offset_per_obs)
    )
  )
  res <- INLA::inla(
    y ~ -1 + offset(offset_val) + f(spatial_idx, model = spde_obj),
    family  = "binomial",
    Ntrials = INLA::inla.stack.data(stk)$Ntrials,
    data    = INLA::inla.stack.data(stk),
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = FALSE),
    control.compute   = list(dic = FALSE, waic = FALSE),
    verbose = FALSE
  )
  as.numeric(res$mlik[1])
}

# -----------------------------------------------------------------------------
# Build SPDE with given (sigma_m, rho) point-fixed prior  (for psi MH)
# -----------------------------------------------------------------------------
make_spde_at <- function(sigma_m, rho) {
  INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(rho,     NA),
    prior.sigma = c(sigma_m, NA)
  )
}

# -----------------------------------------------------------------------------
# Sampler
# -----------------------------------------------------------------------------
ensemble <- ensemble_init(M_TREES, init_val = 0)
sigma_m  <- 0.8
rho      <- 0.5
log_pc_prior <- make_log_pc_prior()

inla_call_count <- 0
sweep_times     <- numeric(N_ITER)

t_total_start <- Sys.time()
for (it in seq_len(N_ITER)) {
  t_sweep <- Sys.time()

  # -------- BART updates --------
  for (j in seq_len(M_TREES)) {
    others_fit <- if (M_TREES == 1) rep(0, N_CL) else
      as.numeric(Reduce(`+`, lapply(ensemble[-j],
                                    function(t) tree_predict(t, dat$X))))

    # current state mlik
    leaf_assign_curr <- tree_leaf_assign(ensemble[[j]], dat$X)
    curr <- fit_binom_for_tree(leaf_assign_curr, others_fit)
    inla_call_count <- inla_call_count + 1

    # propose BD: reuse tree_bd's structure, but ignore its acceptance
    # decision (we'll re-do MH using INLA mliks).
    bd_attempt <- tree_bd(ensemble[[j]], dat$X,
                          r = rep(0, N_CL), Omega = rep(1, N_CL),
                          tau2 = tau2)
    proposed_tree <- bd_attempt$tree
    if (identical(proposed_tree, ensemble[[j]])) {
      # No proposal was actually made (e.g. failed grow). Keep current leaf means.
      leaves <- which(ensemble[[j]]$is_leaf)
      if (length(curr$leaf_means) == length(leaves)) {
        ensemble[[j]]$mu[leaves] <- curr$leaf_means
      }
      next
    }

    # Proposed state mlik
    leaf_assign_prop <- tree_leaf_assign(proposed_tree, dat$X)
    prop <- fit_binom_for_tree(leaf_assign_prop, others_fit)
    inla_call_count <- inla_call_count + 1

    # MH ratio: log lik ratio from INLA's mliks, plus prior+transition ratio
    # already embedded in bd_attempt$log_alpha. The bd_attempt$log_alpha
    # contains log_lik_ratio (from PG-style leaf marginal) + log_prior_ratio +
    # log_trans_ratio. To get the INLA-based MH ratio we subtract the PG lik
    # and add the INLA lik:
    pg_lik_ratio  <- bd_attempt$log_alpha    # contains lik + prior + trans
    # We want: log_alpha_INLA = (prop_mlik - curr_mlik) + log_prior_trans
    # bd_attempt$log_alpha = pg_lik + log_prior_trans
    # So log_prior_trans = bd_attempt$log_alpha - pg_lik_ratio_only.
    # Since pg_lik_ratio_only was already added to bd_attempt$log_alpha, but
    # we don't store it separately, we just use the INLA mliks plus the prior
    # and transition from the bd machinery's *prior+transition* component.
    # For the smoke timing run we approximate by using INLA mliks alone:
    log_alpha_inla <- prop$mlik - curr$mlik

    if (is.finite(log_alpha_inla) && log(runif(1)) < log_alpha_inla) {
      ensemble[[j]] <- proposed_tree
      leaves <- which(ensemble[[j]]$is_leaf)
      if (length(prop$leaf_means) == length(leaves)) {
        ensemble[[j]]$mu[leaves] <- prop$leaf_means
      }
    } else {
      leaves <- which(ensemble[[j]]$is_leaf)
      if (length(curr$leaf_means) == length(leaves)) {
        ensemble[[j]]$mu[leaves] <- curr$leaf_means
      }
    }
  }

  # -------- psi MH --------
  full_offset <- as.numeric(Reduce(`+`, lapply(ensemble,
                                               function(t) tree_predict(t, dat$X))))
  curr_spde  <- make_spde_at(sigma_m, rho)
  mlik_curr  <- fit_binom_for_psi(full_offset, curr_spde)
  inla_call_count <- inla_call_count + 1

  sigma_m_p <- exp(log(sigma_m) + rnorm(1, 0, 0.15))
  rho_p     <- exp(log(rho)     + rnorm(1, 0, 0.15))
  prop_spde <- make_spde_at(sigma_m_p, rho_p)
  mlik_prop <- fit_binom_for_psi(full_offset, prop_spde)
  inla_call_count <- inla_call_count + 1

  log_alpha_psi <- (mlik_prop + log_pc_prior(rho_p, sigma_m_p)) -
                   (mlik_curr + log_pc_prior(rho,   sigma_m))   +
                   (log(sigma_m_p) + log(rho_p)) -
                   (log(sigma_m)   + log(rho))
  if (is.finite(log_alpha_psi) && log(runif(1)) < log_alpha_psi) {
    sigma_m <- sigma_m_p; rho <- rho_p
  }

  sweep_times[it] <- as.numeric(difftime(Sys.time(), t_sweep, units = "secs"))
  cat(sprintf("iter %3d  sweep=%.2fs  sigma_m=%.3f  rho=%.3f  cum_inla_calls=%d\n",
              it, sweep_times[it], sigma_m, rho, inla_call_count))
}
t_total <- as.numeric(difftime(Sys.time(), t_total_start, units = "secs"))

# -----------------------------------------------------------------------------
# Timing report
# -----------------------------------------------------------------------------
cat(sprintf("\n=== Timing summary ===\n"))
cat(sprintf("Total elapsed: %.1f s\n", t_total))
cat(sprintf("Mean per sweep: %.2f s (range %.2f-%.2f)\n",
            mean(sweep_times), min(sweep_times), max(sweep_times)))
cat(sprintf("Total inla() calls: %d\n", inla_call_count))
cat(sprintf("Mean per inla() call: %.3f s\n", t_total / inla_call_count))

# Project to production
prod_iter  <- 1500
prod_trees <- 20
prod_calls_per_sweep <- 2 * (prod_trees + 1)
mean_call <- t_total / inla_call_count
projected_total_sec <- prod_calls_per_sweep * mean_call * prod_iter
cat(sprintf("\nProjected production runtime (n_iter=%d, m_trees=%d):\n",
            prod_iter, prod_trees))
cat(sprintf("  calls per sweep: %d\n", prod_calls_per_sweep))
cat(sprintf("  projected wall time: %.1f h\n", projected_total_sec / 3600))

# Compare to PG-Woodbury Proto 2b's observed 2.5 min for n_iter=1500, m=20:
woodbury_total_sec_at_prod <- 2.5 * 60
cat(sprintf("\nProto 2b (PG-Woodbury) for same config: ~%.1f min\n",
            woodbury_total_sec_at_prod / 60))
cat(sprintf("Speedup factor (Woodbury vs INLA-binom): ~%.0fx\n",
            projected_total_sec / woodbury_total_sec_at_prod))

saveRDS(list(sweep_times = sweep_times, inla_call_count = inla_call_count,
             total_sec = t_total, mean_call_sec = mean_call,
             config = list(N_CL = N_CL, N_ITER = N_ITER, M_TREES = M_TREES,
                           SMOKE_TEST = SMOKE_TEST)),
        sprintf("dev/08_inla_binom_within_mcmc_%s.rds",
                if (SMOKE_TEST) "smoke" else "full"))
cat(sprintf("\nResults saved to dev/08_inla_binom_within_mcmc_%s.rds\n",
            if (SMOKE_TEST) "smoke" else "full"))
