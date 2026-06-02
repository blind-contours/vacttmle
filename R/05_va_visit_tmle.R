###############################################################################
# Thin wrapper for the visit-node-only TMLE.
###############################################################################

va_visit_tmle_contrast <- function(obs_data, ...) {
  va_estimate_contrast(obs_data, estimator = "visit_tmle", ...)
}

