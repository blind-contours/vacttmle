###############################################################################
# ITT/PP Estimator Wrappers
#
# Wrappers for all 7 methods:
# 1. VA-CT-TMLE (C-TMLE ON)   — proposed
# 2. VA-CT-TMLE (plug-in)     — ablation (no C-TMLE)
# 3. VA visit-node TMLE        — no within-interval targeting
# 4. VA-CT-GCOMP               — no targeting, bootstrap inference
# 5. IPCW-KM                   — mechanism-only baseline
# 6. RPSFTM/two-stage          — PP switching baseline
# 7. Monthly LTMLE             — artificial-grid baseline
###############################################################################

# --- 1. VA-CT-TMLE with C-TMLE ---
ittpp_va_ct_tmle_contrast <- function(obs_data, estimand = "PP", tau = obs_data$tau,
                                       g0 = 0.5, M = 50, B = 10, max_iter = 5,
                                       seed = NULL, true_rd = NA_real_,
                                       scenario = NA_character_,
                                       replicate_id = NA_integer_) {
  ittpp_estimate_contrast(
    obs_data, estimator = "ct_tmle", estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = max_iter,
    seed = seed, true_rd = true_rd, scenario = scenario,
    replicate_id = replicate_id, update_method = "adaptive",
    ctmle = TRUE
  )
}

# --- 2. VA-CT-TMLE plug-in (no C-TMLE) ---
ittpp_va_ct_tmle_plugin_contrast <- function(obs_data, estimand = "PP",
                                              tau = obs_data$tau,
                                              g0 = 0.5, M = 50, B = 10,
                                              max_iter = 5, seed = NULL,
                                              true_rd = NA_real_,
                                              scenario = NA_character_,
                                              replicate_id = NA_integer_) {
  res <- ittpp_estimate_contrast(
    obs_data, estimator = "ct_tmle", estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = max_iter,
    seed = seed, true_rd = true_rd, scenario = scenario,
    replicate_id = replicate_id, update_method = "adaptive",
    ctmle = FALSE
  )
  res[, estimator := "va_ct_tmle_plugin"]
  res
}

# --- 2b. VA-CT-TMLE with the full (L-dependent) coarsening model and 95th-
# percentile weight truncation. Under positivity stress this keeps the
# L-dependence (low bias) while capping the extreme near-violation tail
# (controlled variance, convergence) -- a more effective variance-control device
# than C-TMLE's collaborative model-selection, which over-simplifies. Truncation
# is negligible where weights are small, so it matches the others off-stress. ---
ittpp_va_ct_tmle_trunc_contrast <- function(obs_data, estimand = "PP",
                                            tau = obs_data$tau,
                                            g0 = 0.5, M = 50, B = 10,
                                            max_iter = 5, seed = NULL,
                                            true_rd = NA_real_,
                                            scenario = NA_character_,
                                            replicate_id = NA_integer_) {
  old <- getOption("ittpp.weight_trunc", NULL)
  options(ittpp.weight_trunc = "p95")
  on.exit(options(ittpp.weight_trunc = old), add = TRUE)
  res <- ittpp_estimate_contrast(
    obs_data, estimator = "ct_tmle", estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = max_iter,
    seed = seed, true_rd = true_rd, scenario = scenario,
    replicate_id = replicate_id, update_method = "adaptive",
    ctmle = FALSE
  )
  res[, estimator := "va_ct_tmle_trunc"]
  res
}

# --- 3. VA visit-node TMLE ---
ittpp_va_visit_tmle_contrast <- function(obs_data, estimand = "PP",
                                          tau = obs_data$tau,
                                          g0 = 0.5, M = 50, B = 10,
                                          max_iter = 5, seed = NULL,
                                          true_rd = NA_real_,
                                          scenario = NA_character_,
                                          replicate_id = NA_integer_) {
  res <- ittpp_estimate_contrast(
    obs_data, estimator = "visit_tmle", estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = max_iter,
    seed = seed, true_rd = true_rd, scenario = scenario,
    replicate_id = replicate_id, update_method = "adaptive"
  )
  res
}

# --- 4. VA-CT-GCOMP ---
ittpp_va_ct_gcomp_contrast <- function(obs_data, estimand = "PP",
                                        tau = obs_data$tau,
                                        g0 = 0.5, M = 50, B = 10,
                                        seed = NULL, true_rd = NA_real_,
                                        scenario = NA_character_,
                                        replicate_id = NA_integer_) {
  ittpp_estimate_contrast(
    obs_data, estimator = "ct_gcomp", estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = 0,
    seed = seed, true_rd = true_rd, scenario = scenario,
    replicate_id = replicate_id
  )
}

# --- 5. IPCW-KM ---
ittpp_ipcw_km_contrast <- function(obs_data, estimand = "PP",
                                    tau = obs_data$tau, g0 = 0.5,
                                    seed = NULL, true_rd = NA_real_,
                                    scenario = NA_character_,
                                    replicate_id = NA_integer_) {
  t0 <- proc.time()[["elapsed"]]
  nuisance <- fit_ittpp_nuisance(obs_data, tau = tau, g0 = g0, estimand = estimand)
  n <- ittpp_subject_n(obs_data)
  subj <- obs_data$subject

  results_arm <- list()
  for (arm in c(0L, 1L)) {
    # Subjects with A0 == arm, weighted by 1/g0
    long <- ittpp_build_cumulative_weights(obs_data, nuisance, arm, estimand)

    # Last row per subject
    last_rows <- long[order(id, k), .SD[.N], by = id]

    # Survival indicator: T > tau (failure-only)
    surv <- as.integer(subj$T_time > tau)

    w <- numeric(n)
    w[match(last_rows$id, subj$id)] <- last_rows$Cplus
    denom <- sum(w)
    psi <- if (denom > 0) sum(w * surv) / denom else NA_real_
    psi <- va_bound(psi)
    mean_w <- mean(w)
    eif <- if (is.finite(mean_w) && mean_w > 0) (w / mean_w) * (surv - psi) else rep(NA_real_, n)

    results_arm[[as.character(arm)]] <- list(
      psi_surv = psi, risk = 1 - psi, eif = eif,
      se_surv = stats::sd(eif, na.rm = TRUE) / sqrt(n),
      converged = is.finite(psi),
      weights = ittpp_weight_diagnostics(last_rows, arm)
    )
  }

  f0 <- results_arm[["0"]]
  f1 <- results_arm[["1"]]
  eif_delta <- f0$eif - f1$eif
  se <- stats::sd(eif_delta, na.rm = TRUE) / sqrt(n)
  rd <- f1$risk - f0$risk
  ci <- rd + c(-1, 1) * 1.96 * se
  elapsed <- proc.time()[["elapsed"]] - t0

  data.table(
    scenario = scenario, replicate_id = replicate_id, n = n, tau = tau,
    estimator = "ipcw_km", estimand = estimand, update_method = "none",
    psi_surv_d0 = f0$psi_surv, psi_surv_d1 = f1$psi_surv,
    risk_d0 = f0$risk, risk_d1 = f1$risk,
    risk_difference = rd, true_risk_difference = true_rd,
    bias = rd - true_rd, se = se,
    ci_lower = ci[1], ci_upper = ci[2],
    covered = if (is.finite(true_rd)) ci[1] <= true_rd && true_rd <= ci[2] else NA,
    runtime_seconds = elapsed, converged = f0$converged && f1$converged,
    n_outer_iterations = 0L,
    max_abs_eif_T = NA_real_, max_abs_eif_J = NA_real_,
    max_abs_eif_Q = NA_real_, max_abs_total_eif = abs(mean(eif_delta, na.rm = TRUE)),
    mean_weight_d0 = f0$weights$mean_weight, mean_weight_d1 = f1$weights$mean_weight,
    max_weight_d0 = f0$weights$max_weight, max_weight_d1 = f1$weights$max_weight,
    ess_d0 = f0$weights$ess, ess_d1 = f1$weights$ess
  )
}

# --- 6. RPSFTM / two-stage (PP only) ---
ittpp_rpsftm_contrast <- function(obs_data, estimand = "PP",
                                   tau = obs_data$tau, g0 = 0.5,
                                   seed = NULL, true_rd = NA_real_,
                                   scenario = NA_character_,
                                   replicate_id = NA_integer_) {
  # RPSFTM (Robins--Tsiatis): rank-preserving structural accelerated failure-time
  # model with g-estimation of the acceleration factor psi and recensoring.
  t0 <- proc.time()[["elapsed"]]
  n <- ittpp_subject_n(obs_data)
  subj <- obs_data$subject

  # Observed follow-up to failure or administrative/loss censoring (switch is NOT
  # an event; it changes the treatment received, not the at-risk status).
  cens_time <- pmin(subj$C_time, tau)                 # loss-to-follow-up or horizon
  fu        <- pmin(subj$T_time, cens_time)           # observed follow-up time
  fail      <- as.integer(subj$T_time <= cens_time)   # failure-only event indicator

  # Decompose observed follow-up into time on assigned regimen and time after switch
  Rt       <- subj$R_time
  time_on  <- pmin(Rt, fu)                            # before switch (or never switch)
  time_off <- pmax(0, fu - Rt)                        # after switch

  use_surv <- requireNamespace("survival", quietly = TRUE)

  # g-estimation: psi balances the counterfactual *untreated-scale* time across arms.
  # U(psi) = time_off + time_on * exp(psi); recensor at C*(psi) = C * min(1, exp(psi))
  # to avoid psi-induced informative censoring (Robins--Tsiatis recensoring).
  psi_grid <- seq(-1.5, 1.5, by = 0.02)
  z_of <- function(psi) {
    if (!use_surv) return(Inf)
    U     <- time_off + time_on * exp(psi)
    Cstar <- cens_time * min(1, exp(psi))
    delta <- as.integer(fail == 1L & U <= Cstar)
    Ucens <- pmin(U, Cstar)
    if (sum(delta) < 2L || length(unique(subj$A0)) < 2L) return(Inf)
    chi <- tryCatch(survival::survdiff(survival::Surv(Ucens, delta) ~ subj$A0)$chisq,
                    error = function(e) Inf)
    chi
  }
  stats <- vapply(psi_grid, z_of, numeric(1))
  best_psi <- if (all(!is.finite(stats))) 0 else psi_grid[which.min(stats)]

  # PP "assigned arm a and never switch": map every subject's observed time onto the
  # always-on-treatment scale by accelerating the off-treatment portion by exp(-psi).
  U_on  <- time_on + time_off * exp(-best_psi)
  Con   <- cens_time * max(1, exp(-best_psi))         # recensor on the on-treatment scale
  delta_on <- as.integer(fail == 1L & U_on <= Con)
  U_on_c   <- pmin(U_on, Con)

  surv_arm <- function(a) {
    sel <- subj$A0 == a
    if (use_surv && sum(delta_on[sel]) > 0) {
      kf <- tryCatch(
        survival::survfit(survival::Surv(U_on_c[sel], delta_on[sel]) ~ 1),
        error = function(e) NULL
      )
      if (!is.null(kf)) {
        s <- tryCatch(summary(kf, times = tau, extend = TRUE)$surv, error = function(e) NA_real_)
        if (length(s) && is.finite(s)) return(va_bound(s))
      }
    }
    va_bound(mean(U_on_c[sel] > tau, na.rm = TRUE))
  }
  surv_d0 <- surv_arm(0L)
  surv_d1 <- surv_arm(1L)

  rd <- (1 - surv_d1) - (1 - surv_d0)
  se <- NA_real_  # Proper RPSFTM inference requires bootstrap; reported as NA
  elapsed <- proc.time()[["elapsed"]] - t0

  data.table(
    scenario = scenario, replicate_id = replicate_id, n = n, tau = tau,
    estimator = "rpsftm", estimand = estimand, update_method = "none",
    psi_surv_d0 = surv_d0, psi_surv_d1 = surv_d1,
    risk_d0 = 1 - surv_d0, risk_d1 = 1 - surv_d1,
    risk_difference = rd, true_risk_difference = true_rd,
    bias = rd - true_rd, se = se,
    ci_lower = NA_real_, ci_upper = NA_real_,
    covered = NA, runtime_seconds = elapsed, converged = TRUE,
    n_outer_iterations = 0L,
    max_abs_eif_T = NA_real_, max_abs_eif_J = NA_real_,
    max_abs_eif_Q = NA_real_, max_abs_total_eif = NA_real_,
    mean_weight_d0 = NA_real_, mean_weight_d1 = NA_real_,
    max_weight_d0 = NA_real_, max_weight_d1 = NA_real_,
    ess_d0 = NA_real_, ess_d1 = NA_real_
  )
}

# --- 7. Monthly LTMLE ---
ittpp_monthly_ltmle_contrast <- function(obs_data, estimand = "PP",
                                          tau = obs_data$tau, g0 = 0.5,
                                          M = 50, seed = NULL,
                                          true_rd = NA_real_,
                                          scenario = NA_character_,
                                          replicate_id = NA_integer_) {
  # Monthly grid LTMLE: expand to person-month, LOCF covariates
  t0 <- proc.time()[["elapsed"]]
  n <- ittpp_subject_n(obs_data)
  subj <- obs_data$subject
  month_grid <- 0:floor(tau)

  # Expand to person-month
  pm_list <- list()
  pos <- 1L
  for (i in seq_len(n)) {
    first_event <- subj$first_event_time[i]
    # Get visit covariates for LOCF
    pi_rows <- obs_data$person_interval[id == i]
    L_locf <- if (nrow(pi_rows)) pi_rows$L_k[1] else 0

    for (m in seq_along(month_grid)[-1]) {
      t_start <- month_grid[m - 1]
      t_end <- month_grid[m]
      if (first_event <= t_start + 1e-12) break

      # LOCF: update L if there's a visit covariate for this month
      matching_row <- pi_rows[t_start >= pi_rows$t_start - 1e-8 &
                                t_start < pi_rows$t_end + 1e-8]
      if (nrow(matching_row)) L_locf <- matching_row$L_k[1]

      # Carry the subject's (T, R, C) times; ittpp_analysis_long recomputes the
      # estimand-specific event indicators / risk set on this monthly grid.
      pm_list[[pos]] <- data.table(
        id = i, k = m - 2L, t_start = t_start, t_end = t_end, ell = 1,
        W1 = subj$W1[i], W2 = subj$W2[i], A0 = subj$A0[i],
        L_k = L_locf, L_next = NA_real_,
        T_time = subj$T_time[i], R_time = subj$R_time[i], C_time = subj$C_time[i]
      )
      pos <- pos + 1L
    }
  }

  if (pos <= 1L) {
    elapsed <- proc.time()[["elapsed"]] - t0
    return(data.table(
      scenario = scenario, replicate_id = replicate_id, n = n, tau = tau,
      estimator = "monthly_ltmle", estimand = estimand, update_method = "none",
      psi_surv_d0 = NA_real_, psi_surv_d1 = NA_real_,
      risk_d0 = NA_real_, risk_d1 = NA_real_,
      risk_difference = NA_real_, true_risk_difference = true_rd,
      bias = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      covered = NA, runtime_seconds = elapsed, converged = FALSE,
      n_outer_iterations = 0L,
      max_abs_eif_T = NA_real_, max_abs_eif_J = NA_real_,
      max_abs_eif_Q = NA_real_, max_abs_total_eif = NA_real_,
      mean_weight_d0 = NA_real_, mean_weight_d1 = NA_real_,
      max_weight_d0 = NA_real_, max_weight_d1 = NA_real_,
      ess_d0 = NA_real_, ess_d1 = NA_real_
    ))
  }

  pm_data <- rbindlist(pm_list[seq_len(pos - 1L)])

  # Build obs_data-like structure for the monthly grid
  monthly_obs <- list(
    subject = subj,
    person_interval = pm_data,
    visit_L = obs_data$visit_L,
    visit_times = month_grid,
    tau = tau,
    params = obs_data$params
  )

  # Run the ITT/PP estimator on the monthly grid
  res <- tryCatch(
    ittpp_estimate_contrast(
      monthly_obs, estimator = "ct_tmle", estimand = estimand,
      tau = tau, g0 = g0, M = M, B = 5, max_iter = 3,
      seed = seed, true_rd = true_rd, scenario = scenario,
      replicate_id = replicate_id, update_method = "adaptive"
    ),
    error = function(e) NULL
  )

  elapsed <- proc.time()[["elapsed"]] - t0

  if (is.null(res)) {
    return(data.table(
      scenario = scenario, replicate_id = replicate_id, n = n, tau = tau,
      estimator = "monthly_ltmle", estimand = estimand, update_method = "none",
      psi_surv_d0 = NA_real_, psi_surv_d1 = NA_real_,
      risk_d0 = NA_real_, risk_d1 = NA_real_,
      risk_difference = NA_real_, true_risk_difference = true_rd,
      bias = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      covered = NA, runtime_seconds = elapsed, converged = FALSE,
      n_outer_iterations = 0L,
      max_abs_eif_T = NA_real_, max_abs_eif_J = NA_real_,
      max_abs_eif_Q = NA_real_, max_abs_total_eif = NA_real_,
      mean_weight_d0 = NA_real_, mean_weight_d1 = NA_real_,
      max_weight_d0 = NA_real_, max_weight_d1 = NA_real_,
      ess_d0 = NA_real_, ess_d1 = NA_real_
    ))
  }

  res[, estimator := "monthly_ltmle"]
  res[, runtime_seconds := elapsed]
  res
}
