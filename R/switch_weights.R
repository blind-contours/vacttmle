###############################################################################
# ITT/PP Cumulative Inverse Treatment-and-Coarsening Weights
#
# Implements Eq. (weight) from the manuscript:
#   C_k^+(u;d) = I(A0=a)/g0 * prod_{J in Jset} 1/G_k^J(u;d)
#
# ITT: Jset = {C}; PP: Jset = {R, C}
###############################################################################

ittpp_coarsening_survival <- function(long, nuisance, cause = c("R", "C"), d) {
  # Compute G_k^J(ell_k) = exp(-integral lambda^J_k) for each person-interval
  cause <- match.arg(cause)
  nd <- ittpp_set_intervention(long, d)
  lambda_J <- ittpp_predict_hazard(nuisance, nd, cause = cause)
  # Within-interval survival: exp(-lambda^J * time_at_risk)
  # For the full interval boundary: exp(-lambda^J * ell)
  G_interval <- exp(-lambda_J * long$ell)
  G_interval
}

ittpp_cumulative_G <- function(long, nuisance, cause, d) {
  # Compute cumulative G through interval k for each subject
  # Returns a data.table with (id, k, G_cum)
  nd <- ittpp_set_intervention(long, d)
  lambda_J <- ittpp_predict_hazard(nuisance, nd, cause = cause)
  G_interval <- exp(-lambda_J * long$ell)

  # Cumulative product by subject
  dt <- data.table(id = long$id, k = long$k, G_int = G_interval)
  setorder(dt, id, k)
  dt[, G_cum := cumprod(G_int), by = id]
  dt
}

ittpp_cplus <- function(rows, d, nuisance, estimand = c("PP", "ITT"),
                         mode = c("projection", "exact"), B = 10) {
  estimand <- match.arg(estimand)
  mode <- match.arg(mode)

  n <- nrow(rows)
  g0 <- nuisance$g0

  # Treatment indicator: I(A0 = a) / g0
  treat_weight <- ifelse(rows$A0 == as.integer(d), 1 / g0, 0)

  if (all(treat_weight == 0)) return(rep(0, n))

  # Determine Jset based on estimand
  Jset <- if (estimand == "PP") c("R", "C") else "C"

  # Compute cumulative coarsening weight at interval boundary
  coarsening_weight <- rep(1, n)

  for (J in Jset) {
    # Check if this cause exists in data
    fit <- switch(J, R = nuisance$hazard_R, C = nuisance$hazard_C)
    if (is.null(fit$model) && fit$fallback <= 1e-9) next

    nd <- ittpp_set_intervention(rows, d)
    lambda_J <- ittpp_predict_hazard(nuisance, nd, cause = J)

    if (mode == "projection") {
      # Use interval boundary value G_k^J(ell_k)
      # Need cumulative G up to and including this interval
      G_this <- exp(-lambda_J * rows$ell)
      coarsening_weight <- coarsening_weight * G_this
    } else {
      # Exact: use G at the time_at_risk point
      G_this <- exp(-lambda_J * rows$time_at_risk)
      coarsening_weight <- coarsening_weight * G_this
    }
  }

  # Also need cumulative G from prior intervals
  # This is handled by pre-computing cumulative weights
  # For now, use interval-level weight (cumulative done in targeting)

  w <- treat_weight / pmax(coarsening_weight, 1e-8)
  w
}

ittpp_build_cumulative_weights <- function(obs_data, nuisance, d,
                                            estimand = c("PP", "ITT"),
                                            hazard_eps = NULL) {
  estimand <- match.arg(estimand)
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  Jset <- if (estimand == "PP") c("R", "C") else "C"
  g0 <- nuisance$g0

  # Treatment weight
  treat_weight <- ifelse(long$A0 == as.integer(d), 1 / g0, 0)

  # Compute interval-specific coarsening survivals
  G_products <- rep(1, nrow(long))
  for (J in Jset) {
    fit <- switch(J, R = nuisance$hazard_R, C = nuisance$hazard_C)
    if (is.null(fit$model) && fit$fallback <= 1e-9) next

    nd <- ittpp_set_intervention(long, d)
    lambda_J <- ittpp_predict_hazard(nuisance, nd, cause = J)

    # Apply hazard eps if targeting has updated this cause
    if (!is.null(hazard_eps) && !is.null(hazard_eps[[J]])) {
      eps_vec <- hazard_eps[[J]]
      eps_k <- eps_vec[as.character(long$k)]
      eps_k[is.na(eps_k)] <- 0
      # The eps correction is not applied to the coarsening survival directly
      # It's applied in the clever covariate computation
    }

    G_interval <- exp(-lambda_J * long$ell)
    G_products <- G_products * G_interval
  }

  # Cumulative across intervals per subject
  dt <- data.table(id = long$id, k = long$k,
                    treat_w = treat_weight, G_int = G_products)
  setorder(dt, id, k)
  dt[, G_cum := cumprod(G_int), by = id]
  dt[, Cplus := treat_w / pmax(G_cum, 1e-8)]

  # Optional weight truncation to control variance under positivity stress while
  # keeping a correctly specified (L-dependent) coarsening model. Set via
  # options(ittpp.weight_trunc = cap) for an absolute cap, or "p99" for the 99th
  # percentile of the positive weights.
  trunc <- getOption("ittpp.weight_trunc", NULL)
  if (!is.null(trunc)) {
    pos <- dt$Cplus[dt$Cplus > 0]
    cap <- if (is.character(trunc) && grepl("^p", trunc)) {
      q <- as.numeric(sub("^p", "", trunc)) / 100
      if (length(pos)) as.numeric(stats::quantile(pos, q, names = FALSE)) else Inf
    } else as.numeric(trunc)
    if (is.finite(cap)) dt[Cplus > cap, Cplus := cap]
  }

  long[, Cplus := dt$Cplus]
  long[, G_cum := dt$G_cum]
  long
}

ittpp_weight_diagnostics <- function(long_with_weights, d) {
  w <- long_with_weights$Cplus
  w_active <- w[w > 0]
  if (!length(w_active)) {
    return(data.table(
      mean_weight = 0, max_weight = 0, ess = 0, truncation_rate = 0
    ))
  }
  data.table(
    mean_weight = mean(w_active),
    max_weight  = max(w_active),
    ess         = va_effective_sample_size(w_active),
    truncation_rate = 0
  )
}
