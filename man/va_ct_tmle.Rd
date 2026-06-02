\name{va_ct_tmle}
\alias{va_ct_tmle}
\alias{va_visit_tmle}
\alias{va_ct_gcomp}
\title{Estimate a composite risk difference}
\usage{
va_ct_tmle(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
  g0 = 0.5, M = 100, B = 10, max_iter = 10, seed = 1L,
  endpoint = "composite", update_method = c("adaptive", "standard"),
  step_size = 1, min_step_size = 1e-4, max_line_iter = 10)

va_visit_tmle(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
  g0 = 0.5, M = 100, B = 10, max_iter = 10, seed = 1L,
  endpoint = "composite", update_method = c("adaptive", "standard"),
  step_size = 1, min_step_size = 1e-4, max_line_iter = 10)

va_ct_gcomp(obs_data, d0 = 0L, d1 = 1L, tau = obs_data$tau,
  g0 = 0.5, M = 100, B = 10, seed = 1L, endpoint = "composite")
}
\arguments{
\item{obs_data}{Visit-aligned analysis object.}
\item{d0, d1}{Static binary regimes to compare. The risk difference is risk
under \code{d1} minus risk under \code{d0}.}
\item{tau}{Analysis horizon.}
\item{g0}{Baseline randomization probability for treatment 1.}
\item{M}{Monte Carlo draws for visit-level integration.}
\item{B}{Number of sub-bins for hazard targeting.}
\item{max_iter}{Maximum targeting iterations.}
\item{seed}{Random seed.}
\item{endpoint}{Endpoint type. Only \code{"composite"} is supported.}
\item{update_method}{Targeting update method.}
\item{step_size, min_step_size, max_line_iter}{Adaptive update controls.}
}
\value{A \code{vacttmle} object containing estimates, standard error,
confidence interval, and diagnostics.}
\description{
Estimate a 12-month or user-specified horizon composite risk difference using
VA-CT-TMLE, visit-node TMLE, or visit-aligned continuous-time g-computation.
}
