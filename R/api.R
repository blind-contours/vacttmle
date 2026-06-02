###############################################################################
# Trialist-facing package API.
###############################################################################

required_subject_cols <- function() {
  c("id", "W1", "W2", "A0", "T_time", "D_time", "composite_time",
    "composite_event")
}

required_interval_cols <- function() {
  c("id", "k", "t_start", "t_end", "ell", "W1", "W2", "A0", "A_k",
    "L_k", "time_at_risk", "event_T", "event_D", "event_comp", "Y_end")
}

stop_missing_cols <- function(x, cols, name) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(name, " is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

#' Build a visit-aligned analysis object
#'
#' `as_va_data()` converts subject-level and person-interval tables into the
#' analysis object used by the composite-endpoint VA-CT-TMLE estimators.
#'
#' @param subject_data Subject-level data. Required columns are `id`, `W1`, `W2`,
#'   `A0`, `T_time`, `D_time`, `composite_time`, and `composite_event`.
#' @param interval_data Person-interval data. Required columns are `id`, `k`,
#'   `t_start`, `t_end`, `ell`, `W1`, `W2`, `A0`, `A_k`, `L_k`,
#'   `time_at_risk`, `event_T`, `event_D`, `event_comp`, and `Y_end`.
#' @param visit_times Numeric vector of scheduled visit times.
#' @param tau Analysis horizon.
#' @param visit_L Optional subject-by-visit matrix of visit covariates.
#' @param visit_A Optional subject-by-interval matrix of visit treatment values.
#' @param visit_g Optional subject-by-interval matrix of known treatment
#'   probabilities for the observed treatment path.
#' @param params Optional list of DGM or metadata parameters.
#'
#' @return A list with class `va_data`.
#' @export
as_va_data <- function(subject_data, interval_data, visit_times, tau,
                       visit_L = NULL, visit_A = NULL, visit_g = NULL,
                       params = list()) {
  subject <- as.data.table(subject_data)
  interval <- as.data.table(interval_data)
  stop_missing_cols(subject, required_subject_cols(), "subject_data")
  stop_missing_cols(interval, required_interval_cols(), "interval_data")

  subject <- copy(subject)
  interval <- copy(interval)
  subject[, id := as.integer(id)]
  interval[, id := as.integer(id)]
  interval[, k := as.integer(k)]
  setorder(subject, id)
  setorder(interval, id, k)

  K <- length(visit_times) - 1L
  n <- nrow(subject)
  if (is.null(visit_L)) {
    visit_L <- matrix(NA_real_, nrow = n, ncol = length(visit_times))
    first_rows <- interval[, .SD[1L], by = id]
    visit_L[match(first_rows$id, subject$id), 1L] <- first_rows$L_k
    for (kk in seq_len(K)) {
      next_col <- paste0("L_", kk)
      if (next_col %in% names(subject)) visit_L[, kk + 1L] <- subject[[next_col]]
    }
  }
  if (is.null(visit_A)) {
    visit_A <- matrix(NA_integer_, nrow = n, ncol = K)
    for (kk in seq_len(K)) {
      vals <- interval[k == kk - 1L, .SD[1L], by = id][, .(id, A_k)]
      visit_A[match(vals$id, subject$id), kk] <- vals$A_k
    }
  }
  if (is.null(visit_g)) {
    visit_g <- matrix(1, nrow = n, ncol = K)
    visit_g[, 1L] <- 0.5
  }

  out <- list(
    subject = subject,
    person_interval = interval,
    visit_L = visit_L,
    visit_A = visit_A,
    visit_g = visit_g,
    visit_times = visit_times,
    tau = tau,
    params = params
  )
  class(out) <- c("va_data", class(out))
  validate_va_data(out)
  out
}

#' Validate a visit-aligned analysis object
#'
#' @param obs_data Object produced by `as_va_data()` or `simulate_va_trial()`.
#' @param endpoint Endpoint type. Only `"composite"` is currently supported.
#'
#' @return The validated object, invisibly.
#' @export
validate_va_data <- function(obs_data, endpoint = "composite") {
  endpoint <- match.arg(endpoint, "composite")
  if (!is.list(obs_data)) stop("obs_data must be a list-like va_data object.", call. = FALSE)
  needed <- c("subject", "person_interval", "visit_times", "tau")
  missing <- setdiff(needed, names(obs_data))
  if (length(missing)) {
    stop("obs_data is missing required components: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  subject <- as.data.table(obs_data$subject)
  interval <- as.data.table(obs_data$person_interval)
  stop_missing_cols(subject, required_subject_cols(), "obs_data$subject")
  stop_missing_cols(interval, required_interval_cols(), "obs_data$person_interval")
  if (anyDuplicated(subject$id)) stop("subject ids must be unique.", call. = FALSE)
  if (any(interval$time_at_risk < -1e-12, na.rm = TRUE)) {
    stop("time_at_risk must be non-negative.", call. = FALSE)
  }
  if (any(interval$event_T == 1L & interval$event_D == 1L, na.rm = TRUE)) {
    stop("event_T and event_D cannot both equal 1 in the same interval.", call. = FALSE)
  }
  if (!all(interval$event_comp == as.integer(interval$event_T == 1L | interval$event_D == 1L), na.rm = TRUE)) {
    stop("event_comp must equal event_T OR event_D for the composite endpoint.", call. = FALSE)
  }
  if (is.null(obs_data$tau) || length(obs_data$tau) != 1L || !is.finite(obs_data$tau)) {
    stop("obs_data$tau must be a single finite value.", call. = FALSE)
  }
  invisible(obs_data)
}

estimate_composite <- function(obs_data, estimator, d0, d1, tau, g0, M, B,
                               max_iter, seed, update_method, step_size,
                               min_step_size, max_line_iter) {
  validate_va_data(obs_data, endpoint = "composite")
  d0 <- as.integer(d0)
  d1 <- as.integer(d1)
  if (!identical(d0, 0L) || !identical(d1, 1L)) {
    stop("This first package version reports risk under d1 = 1 minus risk under d0 = 0.",
         call. = FALSE)
  }
  res <- va_estimate_contrast(
    obs_data = obs_data,
    estimator = estimator,
    tau = tau,
    g0 = g0,
    M = M,
    B = B,
    max_iter = max_iter,
    seed = seed %||% 1L,
    update_method = update_method,
    step_size = step_size,
    min_step_size = min_step_size,
    max_line_iter = max_line_iter
  )
  out <- list(
    call = match.call(),
    endpoint = "composite",
    estimator = estimator,
    result = res,
    risk_difference = res$risk_difference,
    se = res$se,
    ci = c(lower = res$ci_lower, upper = res$ci_upper),
    diagnostics = res[, .(
      converged, n_outer_iterations, max_abs_eif_T, max_abs_eif_D,
      max_abs_eif_Q, max_abs_total_eif, max_weight_d0, max_weight_d1,
      ess_d0, ess_d1
    )]
  )
  class(out) <- "vacttmle"
  out
}

#' Estimate a composite risk difference with VA-CT-TMLE
#'
#' @param obs_data Visit-aligned analysis object.
#' @param d0,d1 Static binary regimes to compare. The reported risk difference is
#'   risk under `d1` minus risk under `d0`.
#' @param tau Analysis horizon.
#' @param g0 Baseline randomization probability for treatment 1.
#' @param M Monte Carlo draws for visit-level integration.
#' @param B Number of sub-bins for hazard targeting.
#' @param max_iter Maximum targeting iterations.
#' @param seed Random seed.
#' @param endpoint Endpoint type. Only `"composite"` is supported.
#' @param update_method Targeting update method.
#' @param step_size,min_step_size,max_line_iter Adaptive update controls.
#'
#' @return A `vacttmle` object.
#' @export
va_ct_tmle <- function(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
                       g0 = 0.5, M = 100, B = 10, max_iter = 10,
                       seed = 1L, endpoint = "composite",
                       update_method = c("adaptive", "standard"),
                       step_size = 1, min_step_size = 1e-4,
                       max_line_iter = 10) {
  endpoint <- match.arg(endpoint, "composite")
  update_method <- match.arg(update_method)
  estimate_composite(obs_data, "ct_tmle", d0, d1, tau, g0, M, B, max_iter,
                     seed, update_method, step_size, min_step_size,
                     max_line_iter)
}

#' Estimate a composite risk difference with visit-node TMLE
#'
#' @inheritParams va_ct_tmle
#' @return A `vacttmle` object.
#' @export
va_visit_tmle <- function(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
                          g0 = 0.5, M = 100, B = 10, max_iter = 10,
                          seed = 1L, endpoint = "composite",
                          update_method = c("adaptive", "standard"),
                          step_size = 1, min_step_size = 1e-4,
                          max_line_iter = 10) {
  endpoint <- match.arg(endpoint, "composite")
  update_method <- match.arg(update_method)
  estimate_composite(obs_data, "visit_tmle", d0, d1, tau, g0, M, B, max_iter,
                     seed, update_method, step_size, min_step_size,
                     max_line_iter)
}

#' Estimate a composite risk difference with VA-CT g-computation
#'
#' @inheritParams va_ct_tmle
#' @return A `vacttmle` object.
#' @export
va_ct_gcomp <- function(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
                        g0 = 0.5, M = 100, B = 10, seed = 1L,
                        endpoint = "composite") {
  endpoint <- match.arg(endpoint, "composite")
  estimate_composite(obs_data, "ct_gcomp", d0, d1, tau, g0, M, B, 0L,
                     seed, "adaptive", 1, 1e-4, 10)
}

#' Simulate a visit-aligned trial
#'
#' This helper is intended for examples, tests, and teaching. It returns the same
#' `va_data` object accepted by the estimators.
#'
#' @param n Number of subjects.
#' @param scenario Scenario name. Currently `"S1"`, `"S3"`, `"S4c"`, and `"S6"`
#'   are available.
#' @param intervention Optional static intervention, `0` or `1`. Leave `NULL`
#'   to sample treatment from the observed treatment mechanism.
#' @param seed Optional random seed.
#'
#' @return A `va_data` object.
#' @export
simulate_va_trial <- function(n = 1000, scenario = "S1", intervention = NULL,
                              seed = NULL) {
  p <- va_get_scenario(scenario)
  out <- simulate_va_dgm(n, p, intervention = intervention, seed = seed)
  class(out) <- c("va_data", class(out))
  validate_va_data(out)
}

#' Run lightweight package toy checks
#'
#' @param n Analytic toy-simulation sample size.
#' @param seed Random seed.
#'
#' @return A data.table with toy-check estimates and pass/fail flags.
#' @export
run_vacttmle_toy_checks <- function(n = 5000, seed = 7001) {
  dat1 <- simulate_toy_constant(n, c(0, 12), lambda_T = 0.04, lambda_D = 0,
                                treatment = FALSE, seed = seed)
  fit1 <- va_estimate_survival_d(dat1, d = 0L, estimator = "ct_tmle",
                                 tau = 12, g0 = 1, M = 10, B = 10,
                                 max_iter = 1, seed = seed)
  truth1 <- exp(-0.04 * 12)

  dat2 <- simulate_toy_constant(n, c(0, 6, 12), lambda_T = c(0.04, 0.02),
                                lambda_D = 0, treatment = FALSE,
                                seed = seed + 1L)
  fit2 <- va_estimate_survival_d(dat2, d = 0L, estimator = "ct_tmle",
                                 tau = 12, g0 = 1, M = 10, B = 10,
                                 max_iter = 1, seed = seed + 1L)
  truth2 <- exp(-0.04 * 6) * exp(-0.02 * 6)

  dat3 <- simulate_toy_constant(n, c(0, 6, 12), lambda_T = 0.04,
                                lambda_D = 0.02, treatment = FALSE,
                                seed = seed + 2L)
  fit3 <- va_estimate_survival_d(dat3, d = 0L, estimator = "ct_tmle",
                                 tau = 12, g0 = 1, M = 10, B = 10,
                                 max_iter = 1, seed = seed + 2L)
  truth3 <- exp(-(0.04 + 0.02) * 12)

  estimate <- c(fit1$psi_surv, fit2$psi_surv, fit3$psi_surv)
  truth <- c(truth1, truth2, truth3)
  abs_error <- abs(estimate - truth)
  data.table(
    check = c("one_interval", "two_interval", "composite_hazards"),
    estimate = estimate,
    truth = truth,
    abs_error = abs_error,
    passed = abs_error < 0.02
  )
}

#' @export
print.vacttmle <- function(x, ...) {
  res <- x$result
  cat("Visit-aligned continuous-time estimator\n")
  cat("Endpoint: ", x$endpoint, "\n", sep = "")
  cat("Estimator: ", x$estimator, "\n", sep = "")
  cat(sprintf("Risk(d0): %.4f\n", res$risk_d0))
  cat(sprintf("Risk(d1): %.4f\n", res$risk_d1))
  cat(sprintf("Risk difference: %.4f\n", res$risk_difference))
  cat(sprintf("SE: %.4f\n", res$se))
  cat(sprintf("95%% CI: [%.4f, %.4f]\n", res$ci_lower, res$ci_upper))
  cat(sprintf("Converged: %s\n", ifelse(isTRUE(res$converged), "yes", "no")))
  invisible(x)
}
