#' vacttmle: Visit-Aligned Continuous-Time TMLE for Trials with Treatment
#' Switching and Censoring
#'
#' Visit-aligned continuous-time targeted minimum loss estimation for
#' scheduled-visit randomized trials in which treatment switching and loss to
#' follow-up occur in continuous time while covariates are recorded only at
#' scheduled visits. The primary interface is [va_ct_switch()], which estimates
#' the ICH E9(R1) treatment-policy (ITT) and hypothetical no-switching
#' (per-protocol) failure-free survival contrasts from exact switch and dropout
#' dates. Composite-endpoint estimators ([va_ct_tmle()]) are also provided.
#'
#' @keywords internal
"_PACKAGE"
