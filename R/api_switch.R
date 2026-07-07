###############################################################################
# Trialist-facing API for the treatment-switching estimands (ITT / PP).
#
# Estimand map (ICH E9(R1)):
#   estimand = "ITT" -> treatment-policy: switching left natural, censoring
#                       removed. P(T^{a, no-C} > tau).
#   estimand = "PP"  -> hypothetical no-switching: switching AND censoring
#                       removed. P(T^{a, no-R, no-C} > tau).
# Both use the corrected efficient influence function (the intervened switch and
# censoring intensities enter only through the cumulative inverse weight; they
# carry no score terms of their own).
###############################################################################

required_switch_subject_cols <- function() {
  c("id", "A0", "W1", "W2", "T_time", "R_time", "C_time")
}

#' Build a visit-aligned analysis object for the switching estimands
#'
#' Converts a one-row-per-subject table of exact event times plus a
#' subject-by-visit matrix of the time-varying covariate into the object
#' consumed by [va_ct_switch()].
#'
#' The current version uses two baseline covariates (`W1`, `W2`) and a single
#' scheduled-visit covariate `L` (measured at each `visit_times`), matching the
#' estimator in McCoy (2026). Map your trial's covariates onto these columns.
#'
#' @param subject_data Data frame with one row per subject and columns `id`,
#'   `A0` (randomized arm, 0/1), `W1`, `W2` (baseline covariates), `T_time`
#'   (failure time), `R_time` (treatment-switch time; `Inf` if never),
#'   `C_time` (censoring/dropout time; `Inf` if never). Times are on the same
#'   scale as `visit_times` and `tau`; administrative censoring beyond `tau`
#'   should be coded as `C_time = Inf`.
#' @param visit_L Numeric matrix, one row per subject (row order matching
#'   `subject_data`) and one column per entry of `visit_times`, giving the
#'   visit covariate `L` measured at each scheduled visit. Values after a
#'   subject leaves observation may be `NA`.
#' @param visit_times Numeric vector of scheduled visit times, starting at 0.
#' @param tau Analysis horizon.
#'
#' @return A `va_switch_data` object.
#' @export
as_va_switch_data <- function(subject_data, visit_L, visit_times,
                              tau = max(visit_times)) {
  s <- as.data.frame(subject_data)
  stop_missing_cols(as.data.table(s), required_switch_subject_cols(),
                    "subject_data")
  visit_L <- as.matrix(visit_L)
  if (nrow(visit_L) != nrow(s)) {
    stop("visit_L must have one row per subject (", nrow(s), ").", call. = FALSE)
  }
  if (ncol(visit_L) != length(visit_times)) {
    stop("visit_L must have one column per visit time (", length(visit_times),
         ").", call. = FALSE)
  }
  out <- make_ittpp_long(
    W1 = as.numeric(s$W1), W2 = as.numeric(s$W2), A0 = as.integer(s$A0),
    L_mat = visit_L,
    T_time = as.numeric(s$T_time), R_time = as.numeric(s$R_time),
    C_time = as.numeric(s$C_time),
    visit_times = visit_times, tau = tau
  )
  attr(out, "original_ids") <- s$id
  class(out) <- c("va_switch_data", class(out))
  out
}

#' Estimate ITT or per-protocol failure-free survival with VA-CT-TMLE
#'
#' @param obs_data A `va_switch_data` object from [as_va_switch_data()] or
#'   [simulate_va_switch_trial()].
#' @param estimand `"ITT"` (treatment-policy: natural switching, censoring
#'   removed) or `"PP"` (hypothetical no-switching: switching and censoring
#'   removed).
#' @param tau Analysis horizon.
#' @param g0 Known randomization probability of arm 1.
#' @param weight_trunc Optional cumulative-weight control under positivity
#'   stress: a numeric absolute cap, or a string like `"p95"` for the 95th
#'   percentile of the positive weights (recommended when switching is heavy).
#'   `NULL` (default) applies no truncation.
#' @param M Monte Carlo draws for visit-level integration.
#' @param B Sub-bins for within-interval hazard targeting.
#' @param max_iter Maximum targeting iterations.
#' @param seed Random seed.
#'
#' @return A `vacttmle` object; the reported risk difference is
#'   risk under arm 1 minus risk under arm 0 (failure by `tau`).
#' @export
va_ct_switch <- function(obs_data, estimand = c("ITT", "PP"),
                         tau = obs_data$tau, g0 = 0.5, weight_trunc = NULL,
                         M = 100, B = 10, max_iter = 10, seed = 1L) {
  estimand <- match.arg(estimand)
  if (!inherits(obs_data, "va_switch_data")) {
    stop("obs_data must be a va_switch_data object (see as_va_switch_data()).",
         call. = FALSE)
  }
  if (!is.null(weight_trunc)) {
    old <- getOption("ittpp.weight_trunc"); on.exit(options(ittpp.weight_trunc = old))
    options(ittpp.weight_trunc = weight_trunc)
    res <- ittpp_va_ct_tmle_trunc_contrast(obs_data, estimand, tau, g0, M, B,
                                           max_iter, seed = seed)
  } else {
    res <- ittpp_va_ct_tmle_contrast(obs_data, estimand, tau, g0, M, B,
                                     max_iter, seed = seed)
  }
  out <- list(
    call = match.call(),
    endpoint = "switching",
    estimand = estimand,
    estimator = "va_ct_tmle",
    result = res,
    risk_difference = res$risk_difference,
    se = res$se,
    ci = c(lower = res$ci_lower, upper = res$ci_upper),
    diagnostics = res[, intersect(c("converged", "n_outer_iterations",
      "max_abs_total_eif", "max_weight_d0", "max_weight_d1", "ess_d0", "ess_d1"),
      names(res)), with = FALSE]
  )
  class(out) <- "vacttmle"
  out
}

#' Simulate a scheduled-visit trial with switching and censoring
#'
#' Convenience generator for examples, tests, and teaching; returns the same
#' `va_switch_data` object accepted by [va_ct_switch()].
#'
#' @param n Number of subjects.
#' @param scenario Scenario name: `"A0"`--`"A7"` from McCoy (2026). `"A2"`
#'   (crossover only) is a good default; `"A7"` is the oncology-crossover
#'   co-occurrence stress scenario.
#' @param seed Optional random seed.
#'
#' @return A `va_switch_data` object.
#' @export
simulate_va_switch_trial <- function(n = 1000, scenario = "A2", seed = NULL) {
  p <- ittpp_get_scenario(scenario)
  out <- simulate_ittpp_dgm(n, p, seed = seed)
  class(out) <- c("va_switch_data", class(out))
  out
}
