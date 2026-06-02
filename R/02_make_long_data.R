###############################################################################
# Person-interval data helpers.
###############################################################################

va_analysis_long <- function(obs_data, tau = obs_data$tau) {
  d <- copy(obs_data$person_interval)
  d <- d[t_start < tau]
  d[, t_end := pmin(t_end, tau)]
  d[, ell := t_end - t_start]
  d <- d[ell > 0 & time_at_risk > 0]
  d[, time_at_risk := pmin(time_at_risk, ell)]
  setorder(d, id, k)
  d
}

va_interval_template <- function(obs_data, tau = obs_data$tau) {
  d <- va_analysis_long(obs_data, tau)
  d[, .(t_start = min(t_start), t_end = max(t_end), ell = max(ell)), by = k][order(k)]
}

va_subject_n <- function(obs_data) nrow(obs_data$subject)

