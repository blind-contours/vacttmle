###############################################################################
# ITT/PP Nuisance Estimation
#
# Fits 3 Poisson GLMs: hazard_T, hazard_R (PP only), hazard_C.
# Adapts from R/03_fit_nuisance.R for the three-cause ITT/PP setting.
###############################################################################

ittpp_analysis_long <- function(obs_data, tau = obs_data$tau,
                                 estimand = c("PP", "ITT")) {
  estimand <- match.arg(estimand)
  d <- copy(obs_data$person_interval)
  d <- d[t_start < tau]
  d[, t_end := pmin(t_end, tau)]
  d[, ell := t_end - t_start]

  # Estimand-specific risk set end:
  #   ITT: switch is NOT terminal; follow-up ends at min(T, C, tau).
  #   PP : censor at the switch; follow-up ends at min(T, R, C, tau).
  if (estimand == "ITT") {
    d[, Xend := pmin(T_time, C_time, tau)]
  } else {
    d[, Xend := pmin(T_time, R_time, C_time, tau)]
  }
  d <- d[t_start < Xend - 1e-12]
  d[, time_at_risk := pmin(Xend, t_end) - t_start]
  d <- d[ell > 0 & time_at_risk > 1e-12]

  # Within-interval cause-specific event indicators
  if (estimand == "ITT") {
    # Switch is absorbed into the natural failure process; not a coarsening event.
    d[, event_T := as.integer(T_time > t_start & T_time <= t_end &
                                T_time <= C_time & T_time <= tau)]
    d[, event_R := 0L]
    d[, event_C := as.integer(C_time > t_start & C_time <= t_end &
                                C_time < T_time & C_time <= tau)]
    d[, event_time := fifelse(event_T == 1L, T_time - t_start,
                       fifelse(event_C == 1L, C_time - t_start, NA_real_))]
  } else {
    d[, event_T := as.integer(T_time > t_start & T_time <= t_end &
                                T_time <= R_time & T_time <= C_time & T_time <= tau)]
    d[, event_R := as.integer(R_time > t_start & R_time <= t_end &
                                R_time < T_time & R_time <= C_time & R_time <= tau)]
    d[, event_C := as.integer(C_time > t_start & C_time <= t_end &
                                C_time < T_time & C_time < R_time & C_time <= tau)]
    d[, event_time := fifelse(event_T == 1L, T_time - t_start,
                       fifelse(event_R == 1L, R_time - t_start,
                        fifelse(event_C == 1L, C_time - t_start, NA_real_)))]
  }

  d[, Y_end := as.integer(Xend > t_end + 1e-12)]
  d[Y_end == 0L, L_next := NA_real_]
  setorder(d, id, k)
  d
}

ittpp_subject_n <- function(obs_data) nrow(obs_data$subject)

fit_ittpp_nuisance <- function(obs_data, tau = obs_data$tau, g0 = 0.5,
                                estimand = c("PP", "ITT"),
                                misspec = NULL) {
  estimand <- match.arg(estimand)
  d <- ittpp_analysis_long(obs_data, tau, estimand)
  candidates <- c("A0", "L_k", "W1", "W2", "k_f")

  # Optionally misspecify for DR matrix (A3/A4)
  # misspec is a list with elements: hazard_T, hazard_R, hazard_C, transition
  #   each TRUE/FALSE indicating whether to misspecify
  hz_T_candidates <- candidates
  hz_R_candidates <- candidates
  hz_C_candidates <- candidates
  trans_candidates <- c("L_k", "A0", "W1", "W2", "k_f")

  if (!is.null(misspec)) {
    if (isTRUE(misspec$hazard_T)) hz_T_candidates <- c("k_f")
    if (isTRUE(misspec$hazard_R)) hz_R_candidates <- c("k_f")
    if (isTRUE(misspec$hazard_C)) hz_C_candidates <- c("k_f")
    if (isTRUE(misspec$transition)) trans_candidates <- c("k_f")
  }

  # Failure hazard
  hz_T <- va_fit_poisson_rate(d, "event_T", hz_T_candidates)

  # Switch hazard (fit for PP; also fit for ITT to allow augmentation check)
  hz_R <- if (any(d$event_R == 1L)) {
    va_fit_poisson_rate(d, "event_R", hz_R_candidates)
  } else {
    list(model = NULL, fallback = 1e-10, formula = NULL)
  }

  # Censoring hazard
  hz_C <- if (any(d$event_C == 1L)) {
    va_fit_poisson_rate(d, "event_C", hz_C_candidates)
  } else {
    list(model = NULL, fallback = 1e-10, formula = NULL)
  }

  # Covariate transition
  trans_data <- d[k < max(k) & Y_end == 1L & !is.na(L_next)]
  trans_data[, k_f := factor(k)]
  trans_fit <- NULL
  trans_sigma <- 0
  if (nrow(trans_data) > 5 && uniqueN(trans_data$L_next) > 1) {
    f <- va_build_formula("L_next", trans_data, trans_candidates)
    trans_fit <- tryCatch(lm(f, data = trans_data), error = function(e) NULL)
    if (!is.null(trans_fit)) trans_sigma <- max(sigma(trans_fit), 1e-6)
  }

  intervals <- d[, .(
    t_start = min(t_start), t_end = max(t_end), ell = max(ell),
    n_rows = .N, n_T = sum(event_T), n_R = sum(event_R), n_C = sum(event_C)
  ), by = k][order(k)]

  list(
    hazard_T      = hz_T,
    hazard_R      = hz_R,
    hazard_C      = hz_C,
    transition_L  = trans_fit,
    transition_sigma = trans_sigma,
    intervals     = intervals,
    g0            = g0,
    tau           = tau,
    estimand      = estimand
  )
}

ittpp_set_intervention <- function(rows, d) {
  nd <- copy(rows)
  nd[, A0 := as.integer(d)]
  nd
}

ittpp_predict_hazard <- function(nuisance, rows, cause = c("T", "R", "C"), d = NULL) {
  cause <- match.arg(cause)
  nd <- copy(rows)
  if (!is.null(d)) nd <- ittpp_set_intervention(nd, d)
  fit <- switch(cause,
    T = nuisance$hazard_T,
    R = nuisance$hazard_R,
    C = nuisance$hazard_C
  )
  va_predict_rate(fit, nd)
}

ittpp_predict_transition_mean <- function(nuisance, rows, d) {
  nd <- ittpp_set_intervention(rows, d)
  nd[, k_f := factor(k)]
  if (is.null(nuisance$transition_L)) return(nd$L_k)
  pred <- tryCatch(
    as.numeric(predict(nuisance$transition_L, newdata = nd)),
    error = function(e) nd$L_k
  )
  pred
}
