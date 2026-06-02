###############################################################################
# G-computation and benchmark comparators.
###############################################################################

va_ct_gcomp_contrast <- function(obs_data, ...) {
  va_estimate_contrast(obs_data, estimator = "ct_gcomp", ...)
}

va_direct_estimators <- function() {
  c("va_discrete_gcomp", "iptw_km", "naive_cox")
}

va_fit_discrete_event_model <- function(obs_data, tau = obs_data$tau) {
  d <- va_analysis_long(obs_data, tau)
  d[, k_f := factor(k)]
  fallback <- va_bound(mean(d$event_comp, na.rm = TRUE))
  if (!nrow(d) || uniqueN(d$event_comp) < 2L) {
    return(list(type = "constant", value = fallback, k_levels = levels(d$k_f)))
  }
  f <- va_build_formula("event_comp", d, c("A_k", "L_k", "W1", "W2", "k_f"))
  fit <- tryCatch(
    suppressWarnings(glm(f, family = binomial(), data = d)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    list(type = "constant", value = fallback, k_levels = levels(d$k_f))
  } else {
    list(type = "glm", fit = fit, value = fallback, k_levels = levels(d$k_f))
  }
}

va_predict_discrete_event <- function(event_model, rows, d) {
  if (!nrow(rows)) return(numeric())
  if (identical(event_model$type, "constant")) {
    return(rep(event_model$value, nrow(rows)))
  }
  nd <- va_set_intervention(rows, d)
  nd[, k_f := factor(k, levels = event_model$k_levels)]
  pred <- tryCatch(
    as.numeric(predict(event_model$fit, newdata = nd, type = "response")),
    error = function(e) rep(event_model$value, nrow(nd))
  )
  va_bound(pred)
}

va_discrete_gcomp_survival_d <- function(obs_data, d, event_model = NULL,
                                         nuisance = NULL, tau = obs_data$tau,
                                         g0 = 0.5, M = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(event_model)) event_model <- va_fit_discrete_event_model(obs_data, tau)
  if (is.null(nuisance)) nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  long <- va_analysis_long(obs_data, tau)
  K <- max(long$k) + 1L
  q_models <- vector("list", K)
  values <- vector("list", K)
  names(values) <- as.character(0:(K - 1L))

  for (kk in rev(0:(K - 1L))) {
    rows <- copy(long[k == kk])
    if (!nrow(rows)) next
    M0 <- if (kk == K - 1L) {
      rep(1, nrow(rows))
    } else {
      va_mc_next_q(rows, d, nuisance, q_models[[kk + 2L]], M = M)
    }
    p_event <- va_predict_discrete_event(event_model, rows, d)
    S <- va_bound(1 - p_event)
    Q <- va_bound(S * M0)
    q_models[[kk + 1L]] <- va_fit_q_model(rows, Q)
    values[[as.character(kk)]] <- data.table(
      id = rows$id, k = kk, S = S, M = M0, Q = Q
    )
  }

  ids <- obs_data$subject$id
  q0 <- values[["0"]][match(ids, id), Q]
  q0[is.na(q0)] <- mean(values[["0"]]$Q)
  psi <- va_bound(mean(q0))
  rows <- va_analysis_long(obs_data, tau)
  list(
    estimator = "va_discrete_gcomp",
    d = d,
    psi_surv = psi,
    risk = 1 - psi,
    eif = q0 - psi,
    se_surv = stats::sd(q0 - psi) / sqrt(va_subject_n(obs_data)),
    converged = TRUE,
    n_outer_iterations = 0L,
    diagnostics = list(
      max_abs_eif_T = NA_real_,
      max_abs_eif_D = NA_real_,
      max_abs_eif_Q = NA_real_,
      max_abs_total_eif = abs(mean(q0 - psi))
    ),
    weights = va_weight_diagnostics(rows, d, nuisance),
    state = list(psi = psi, q_models = q_models, values = values),
    nuisance = nuisance
  )
}

va_subject_survival_data <- function(obs_data, tau = obs_data$tau) {
  subj <- copy(obs_data$subject)
  l0 <- obs_data$visit_L[, 1]
  if (is.null(l0)) {
    first_rows <- va_analysis_long(obs_data, tau)[k == 0L, .(id, L0 = L_k)]
    subj <- merge(subj, first_rows, by = "id", all.x = TRUE, sort = FALSE)
  } else {
    subj[, L0 := l0]
  }
  subj[, time := pmin(composite_time, tau)]
  subj[, event := as.integer(composite_event == 1L & composite_time <= tau)]
  subj
}

va_iptw_km_survival_d <- function(obs_data, d, tau = obs_data$tau, g0 = 0.5,
                                  nuisance = NULL, estimator_label = "iptw_km") {
  subj <- va_subject_survival_data(obs_data, tau)
  surv <- as.integer(subj$T_time > tau & subj$D_time > tau)
  if (is.null(nuisance)) nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  rows <- va_analysis_long(obs_data, tau)
  last_rows <- rows[order(id, k), .SD[.N], by = id]
  w <- numeric(nrow(subj))
  w[match(last_rows$id, subj$id)] <- va_cplus_observed(last_rows, d, nuisance)
  denom <- sum(w)
  psi <- if (denom > 0) sum(w * surv) / denom else NA_real_
  psi <- va_bound(psi)
  mean_w <- mean(w)
  eif <- if (is.finite(mean_w) && mean_w > 0) (w / mean_w) * (surv - psi) else rep(NA_real_, nrow(subj))
  list(
    estimator = estimator_label,
    d = d,
    psi_surv = psi,
    risk = 1 - psi,
    eif = eif,
    se_surv = stats::sd(eif, na.rm = TRUE) / sqrt(va_subject_n(obs_data)),
    converged = is.finite(psi),
    n_outer_iterations = 0L,
    diagnostics = list(
      max_abs_eif_T = NA_real_,
      max_abs_eif_D = NA_real_,
      max_abs_eif_Q = NA_real_,
      max_abs_total_eif = abs(mean(eif, na.rm = TRUE))
    ),
    weights = va_weight_diagnostics(rows, d, nuisance),
    state = list(psi = psi),
    nuisance = nuisance
  )
}

va_fit_naive_cox <- function(obs_data, tau = obs_data$tau) {
  subj <- va_subject_survival_data(obs_data, tau)
  if (!requireNamespace("survival", quietly = TRUE) || sum(subj$event) < 2L) {
    return(list(type = "fallback", subject = subj))
  }
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, event) ~ A0 + L0 + W1 + W2,
      data = subj, ties = "efron", model = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) list(type = "fallback", subject = subj) else list(type = "cox", fit = fit, subject = subj)
}

va_predict_naive_cox_survival <- function(cox_model, d, tau) {
  subj <- copy(cox_model$subject)
  if (!identical(cox_model$type, "cox")) {
    surv <- as.integer(subj$T_time > tau & subj$D_time > tau)
    return(rep(mean(surv[subj$A0 == as.integer(d)], na.rm = TRUE), nrow(subj)))
  }
  nd <- copy(subj)
  nd[, A0 := as.integer(d)]
  sf <- tryCatch(survival::survfit(cox_model$fit, newdata = nd), error = function(e) NULL)
  if (!is.null(sf)) {
    ss <- tryCatch(summary(sf, times = tau, extend = TRUE), error = function(e) NULL)
    if (!is.null(ss) && length(ss$surv) == nrow(nd)) return(va_bound(as.numeric(ss$surv)))
  }
  bh <- tryCatch(survival::basehaz(cox_model$fit, centered = FALSE), error = function(e) NULL)
  lp <- tryCatch(as.numeric(predict(cox_model$fit, newdata = nd, type = "lp", reference = "zero")),
                 error = function(e) NULL)
  if (!is.null(bh) && !is.null(lp) && nrow(bh) && length(lp) == nrow(nd)) {
    idx <- max(which(bh$time <= tau), 0L)
    H0 <- if (idx > 0L) bh$hazard[idx] else 0
    return(va_bound(exp(-H0 * exp(lp))))
  }
  surv <- as.integer(subj$T_time > tau & subj$D_time > tau)
  rep(mean(surv[subj$A0 == as.integer(d)], na.rm = TRUE), nrow(subj))
}

va_naive_cox_survival_d <- function(obs_data, d, cox_model = NULL,
                                    tau = obs_data$tau, g0 = 0.5,
                                    nuisance = NULL) {
  if (is.null(cox_model)) cox_model <- va_fit_naive_cox(obs_data, tau)
  pred <- va_bound(va_predict_naive_cox_survival(cox_model, d, tau))
  psi <- va_bound(mean(pred))
  if (is.null(nuisance)) nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  rows <- va_analysis_long(obs_data, tau)
  list(
    estimator = "naive_cox",
    d = d,
    psi_surv = psi,
    risk = 1 - psi,
    eif = pred - psi,
    se_surv = stats::sd(pred - psi) / sqrt(va_subject_n(obs_data)),
    converged = is.finite(psi),
    n_outer_iterations = 0L,
    diagnostics = list(
      max_abs_eif_T = NA_real_,
      max_abs_eif_D = NA_real_,
      max_abs_eif_Q = NA_real_,
      max_abs_total_eif = abs(mean(pred - psi))
    ),
    weights = va_weight_diagnostics(rows, d, nuisance),
    state = list(psi = psi, cox_type = cox_model$type),
    nuisance = nuisance
  )
}

va_make_direct_contrast_row <- function(obs_data, estimator, fit0, fit1,
                                        tau, true_rd = NA_real_,
                                        scenario = NA_character_,
                                        replicate_id = NA_integer_,
                                        elapsed = NA_real_) {
  eif_delta <- fit0$eif - fit1$eif
  se <- stats::sd(eif_delta, na.rm = TRUE) / sqrt(va_subject_n(obs_data))
  rd <- fit1$risk - fit0$risk
  ci <- rd + c(-1, 1) * 1.96 * se
  data.table(
    scenario = scenario,
    replicate_id = replicate_id,
    n = va_subject_n(obs_data),
    tau = tau,
    estimator = estimator,
    update_method = "none",
    psi_surv_d0 = fit0$psi_surv,
    psi_surv_d1 = fit1$psi_surv,
    risk_d0 = fit0$risk,
    risk_d1 = fit1$risk,
    risk_difference = rd,
    true_risk_difference = true_rd,
    bias = rd - true_rd,
    se = se,
    ci_lower = ci[1],
    ci_upper = ci[2],
    covered = if (is.finite(true_rd)) ci[1] <= true_rd && true_rd <= ci[2] else NA,
    runtime_seconds = elapsed,
    converged = fit0$converged && fit1$converged && is.finite(rd) && is.finite(se),
    n_outer_iterations = 0L,
    max_abs_eif_T = NA_real_,
    max_abs_eif_D = NA_real_,
    max_abs_eif_Q = NA_real_,
    max_abs_total_eif = abs(mean(eif_delta, na.rm = TRUE)),
    mean_weight_d0 = fit0$weights$mean_weight,
    mean_weight_d1 = fit1$weights$mean_weight,
    max_weight_d0 = fit0$weights$max_weight,
    max_weight_d1 = fit1$weights$max_weight,
    ess_d0 = fit0$weights$ess,
    ess_d1 = fit1$weights$ess,
    truncation_rate_d0 = fit0$weights$truncation_rate,
    truncation_rate_d1 = fit1$weights$truncation_rate
  )
}

va_estimate_direct_comparator_contrast <- function(obs_data, estimator, tau = obs_data$tau,
                                                   g0 = 0.5, M = 100, seed = NULL,
                                                   true_rd = NA_real_,
                                                   scenario = NA_character_,
                                                   replicate_id = NA_integer_) {
  estimator <- match.arg(estimator, va_direct_estimators())
  t0 <- proc.time()[["elapsed"]]
  nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  if (identical(estimator, "va_discrete_gcomp")) {
    seed0 <- seed %||% sample.int(.Machine$integer.max, 1L)
    if (!is.null(seed)) set.seed(seed)
    event_model <- va_fit_discrete_event_model(obs_data, tau)
    fit0 <- va_discrete_gcomp_survival_d(obs_data, 0L, event_model, nuisance, tau, g0, M, seed0)
    fit1 <- va_discrete_gcomp_survival_d(obs_data, 1L, event_model, nuisance, tau, g0, M, seed0 + 1L)
  } else if (identical(estimator, "iptw_km")) {
    fit0 <- va_iptw_km_survival_d(obs_data, 0L, tau, g0, nuisance)
    fit1 <- va_iptw_km_survival_d(obs_data, 1L, tau, g0, nuisance)
  } else {
    cox_model <- va_fit_naive_cox(obs_data, tau)
    fit0 <- va_naive_cox_survival_d(obs_data, 0L, cox_model, tau, g0, nuisance)
    fit1 <- va_naive_cox_survival_d(obs_data, 1L, cox_model, tau, g0, nuisance)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  va_make_direct_contrast_row(
    obs_data, estimator, fit0, fit1, tau, true_rd, scenario, replicate_id, elapsed
  )
}
