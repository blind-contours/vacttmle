###############################################################################
# Initial nuisance fits for the target-aligned estimators.
###############################################################################

fit_va_nuisance <- function(obs_data, tau = obs_data$tau, g0 = 0.5) {
  d <- va_analysis_long(obs_data, tau)
  candidates <- c("A_k", "L_k", "W1", "W2", "k_f")
  hz_T <- va_fit_poisson_rate(d, "event_T", candidates)
  hz_D <- va_fit_poisson_rate(d, "event_D", candidates)

  trans_data <- d[k < max(k) & Y_end == 1L & !is.na(L_next)]
  trans_data[, k_f := factor(k)]
  trans_fit <- NULL
  trans_sigma <- 0
  if (nrow(trans_data) > 5 && uniqueN(trans_data$L_next) > 1) {
    f <- va_build_formula("L_next", trans_data, c("L_k", "A_k", "W1", "W2", "k_f"))
    trans_fit <- tryCatch(lm(f, data = trans_data), error = function(e) NULL)
    if (!is.null(trans_fit)) trans_sigma <- max(sigma(trans_fit), 1e-6)
  }

  intervals <- d[, .(
    t_start = min(t_start), t_end = max(t_end), ell = max(ell),
    n_rows = .N, n_T = sum(event_T), n_D = sum(event_D)
  ), by = k][order(k)]

  list(
    hazard_T = hz_T,
    hazard_D = hz_D,
    transition_L = trans_fit,
    transition_sigma = trans_sigma,
    intervals = intervals,
    g0 = g0,
    tau = tau
  )
}

va_set_intervention <- function(rows, d) {
  nd <- copy(rows)
  nd[, A_k := as.integer(d)]
  nd
}

va_predict_hazard_pair <- function(nuisance, rows, d = NULL) {
  nd <- copy(rows)
  if (!is.null(d)) nd <- va_set_intervention(nd, d)
  data.table(
    lambda_T = va_predict_rate(nuisance$hazard_T, nd),
    lambda_D = va_predict_rate(nuisance$hazard_D, nd)
  )
}

va_predict_transition_mean <- function(nuisance, rows, d) {
  nd <- va_set_intervention(rows, d)
  nd[, k_f := factor(k)]
  if (is.null(nuisance$transition_L)) return(nd$L_k)
  pred <- tryCatch(
    as.numeric(predict(nuisance$transition_L, newdata = nd)),
    error = function(e) nd$L_k
  )
  pred
}

va_cplus_observed <- function(rows, d, nuisance, truncate = FALSE) {
  follow_col <- paste0("follow_d", as.integer(d))
  if (follow_col %in% names(rows) && "g_cum_observed" %in% names(rows)) {
    denom <- pmax(rows$g_cum_observed, 1e-8)
    follows <- !is.na(rows[[follow_col]]) & rows[[follow_col]]
    w <- ifelse(follows, 1 / denom, 0)
  } else {
    w <- ifelse(rows$A0 == as.integer(d), 1 / nuisance$g0, 0)
  }
  if (truncate && any(w > 0)) {
    cap <- as.numeric(stats::quantile(w[w > 0], 0.99, names = FALSE))
    w <- pmin(w, cap)
  }
  w
}

va_cplus_intervention <- function(rows, d, nuisance) {
  g_col <- paste0("g_cum_d", as.integer(d))
  if (g_col %in% names(rows)) return(1 / pmax(rows[[g_col]], 1e-8))
  rep(1 / nuisance$g0, nrow(rows))
}

va_weight_diagnostics <- function(rows, d, nuisance) {
  w <- va_cplus_observed(rows, d, nuisance)
  data.table(
    mean_weight = mean(w),
    max_weight = max(w),
    ess = va_effective_sample_size(w),
    truncation_rate = 0
  )
}
