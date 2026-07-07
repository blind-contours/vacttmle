###############################################################################
# ITT/PP EIF Components
#
# Implements the efficient influence function from Theorem 1 of the manuscript.
# Key differences from composite (R/04_va_ct_tmle.R va_eif_components):
#   - h^T_k(u) = -C_k^+(u) * Sbar_T(u) * M_k   (failure, NEGATIVE sign)
#   - Sbar_T uses failure-only forward survival (NOT all-cause)
#
# F1 CORRECTION (2026-07-02): the coarsening scores int h^J dM^J with
# h^J = +C^+ Sbar_T M are NOT part of the EIF: with the failure score in
# martingale form and the visit-level residuals conditionally centered, adding
# them breaks the gradient property (verified analytically and numerically in
# scripts/itt_pp/diag_gradient_check.R; see the math note, Finding F1). They
# are still COMPUTED here (means_J) because the coarsening fluctuation steps
# use them as targeting directions -- solving those extra mean-zero equations
# is a harmless stabilization device -- but they are EXCLUDED from the EIF
# total used for inference (SE) and from the remainder identity.
###############################################################################

ittpp_make_subbins <- function(rows, d, nuisance, state, cause = c("T", "R", "C"),
                                estimand = c("PP", "ITT"), B = 10,
                                use_updated = FALSE, hazard_eps = NULL) {
  cause <- match.arg(cause)
  estimand <- match.arg(estimand)
  if (!nrow(rows)) return(data.table())

  val <- state$values[[as.character(rows$k[1])]]
  rows <- merge(copy(rows), val[, .(id, M)], by = "id", all.x = TRUE, sort = FALSE)

  # Predict failure hazard under intervention d
  base_d <- ittpp_predict_hazard(nuisance, ittpp_set_intervention(rows, d), cause = "T")
  # Predict cause-specific hazard under observed data
  base_obs <- ittpp_predict_hazard(nuisance, rows, cause = cause)

  # Cumulative weight C_k^+ (using observed-data weight for EIF)
  Cobs <- rows$Cplus
  if (is.null(Cobs) || all(is.na(Cobs))) Cobs <- rep(0, nrow(rows))

  n <- nrow(rows)
  idx <- rep(seq_len(n), each = B)
  b_idx <- rep(seq_len(B), times = n)
  ell <- rows$ell[idx]
  lo <- (b_idx - 1) * ell / B
  hi <- b_idx * ell / B
  mid <- (lo + hi) / 2
  pt <- pmax(0, pmin(rows$time_at_risk[idx], hi) - lo)
  keep <- pt > 0
  if (!any(keep)) return(data.table())

  idx <- idx[keep]
  lo  <- lo[keep]
  hi  <- hi[keep]
  mid <- mid[keep]
  ell <- ell[keep]
  pt  <- pt[keep]

  # Sbar_T: failure-only forward survival from u to ell_k
  # Sbar_T(u) = exp(-(ell_k - u) * lambda^T_k)
  sbar_T <- exp(-(ell - mid) * base_d[idx])

  # Within-interval weight: projection (default) uses the interval-boundary
  # C^+_k(ell_k); exact (option a, Section 8.2) uses the elapsed-time C^+_k(u),
  # which equals the boundary weight times exp(-sum_J lambda^J_k (ell_k - u)).
  Cw <- Cobs[idx]
  if (identical(getOption("ittpp.weight_mode", "projection"), "exact")) {
    Jset_w <- if (estimand == "PP") c("R", "C") else "C"
    nd_w <- ittpp_set_intervention(rows, d)
    adj <- numeric(length(rows$id))
    for (Jw in Jset_w) {
      fitJ <- switch(Jw, R = nuisance$hazard_R, C = nuisance$hazard_C)
      if (is.null(fitJ$model) && fitJ$fallback <= 1e-9) next
      adj <- adj + ittpp_predict_hazard(nuisance, nd_w, cause = Jw)
    }
    Cw <- Cw * exp(-adj[idx] * (ell - mid))
  }

  # Clever covariate sign:
  #   Failure T: h = -Cplus * Sbar_T * M  (NEGATIVE)
  #   Coarsening J: h = +Cplus * Sbar_T * M  (POSITIVE)
  sign_factor <- if (cause == "T") -1 else +1
  h <- sign_factor * Cw * sbar_T * rows$M[idx]

  # Event indicator in this sub-bin
  event_here <- !is.na(rows$event_time[idx]) &
    rows$event_time[idx] > lo - 1e-12 & rows$event_time[idx] <= hi + 1e-12

  event_col <- paste0("event_", cause)
  count <- as.integer(rows[[event_col]][idx] == 1L & event_here)

  lambda0 <- base_obs[idx]

  # Apply hazard eps if using updated estimates
  eps_vec <- if (!is.null(hazard_eps)) hazard_eps[[cause]] else NULL
  eps <- (eps_vec %||% numeric())[as.character(rows$k[idx])]
  eps[is.na(eps)] <- 0

  mu0 <- lambda0 * pt
  mu <- mu0 * exp(if (use_updated) eps * h else 0)

  data.table(
    id = rows$id[idx], k = rows$k[idx],
    h = h, count = count, mu0 = mu0, mu = mu
  )
}

ittpp_eif_components <- function(obs_data, nuisance, d, state,
                                  estimand = c("PP", "ITT"), B = 10) {
  estimand <- match.arg(estimand)
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  n <- ittpp_subject_n(obs_data)
  ids <- obs_data$subject$id
  total <- numeric(n)
  # Per-subject contributions of each EIF component, for the variance decomposition
  comp_base <- numeric(n); comp_T <- numeric(n); comp_J <- numeric(n); comp_Q <- numeric(n)

  # Baseline term: Q_0 - psi
  q0 <- state$values[["0"]][match(ids, id), Q]
  q0[is.na(q0)] <- mean(state$values[["0"]]$Q)
  psi <- mean(q0)
  total <- total + q0 - psi
  comp_base <- comp_base + (q0 - psi)

  K <- max(long$k) + 1L
  Jset <- if (estimand == "PP") c("R", "C") else "C"

  means_T <- setNames(rep(0, K), as.character(0:(K - 1L)))
  means_J <- list()
  for (J in Jset) {
    means_J[[J]] <- setNames(rep(0, K), as.character(0:(K - 1L)))
  }
  means_Q <- setNames(rep(0, K), as.character(0:(K - 1L)))

  # Ensure cumulative weights are on the long data
  if (!"Cplus" %in% names(long)) {
    long <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                            hazard_eps = state$hazard_eps)
  }

  for (kk in 0:(K - 1L)) {
    rows <- long[k == kk]

    # Failure score
    sb_T <- ittpp_make_subbins(rows, d, nuisance, state, "T", estimand, B,
                                use_updated = TRUE, hazard_eps = state$hazard_eps)
    if (nrow(sb_T)) {
      contr <- sb_T[, .(v = sum(h * (count - mu))), by = id]
      total[match(contr$id, ids)] <- total[match(contr$id, ids)] + contr$v
      comp_T[match(contr$id, ids)] <- comp_T[match(contr$id, ids)] + contr$v
      means_T[as.character(kk)] <- sum(contr$v) / n
    }

    # Coarsening scores (for each J in Jset)
    for (J in Jset) {
      event_col <- paste0("event_", J)
      if (!event_col %in% names(rows) || !any(rows[[event_col]] == 1L, na.rm = TRUE)) {
        # No events of this cause — check if hazard is nonzero
        fit <- switch(J, R = nuisance$hazard_R, C = nuisance$hazard_C)
        if (is.null(fit$model) && fit$fallback <= 1e-9) next
      }
      sb_J <- ittpp_make_subbins(rows, d, nuisance, state, J, estimand, B,
                                  use_updated = TRUE, hazard_eps = state$hazard_eps)
      if (nrow(sb_J)) {
        contr <- sb_J[, .(v = sum(h * (count - mu))), by = id]
        # F1: coarsening scores are targeting directions only -- NOT added to
        # the EIF total (they are not part of the canonical gradient).
        comp_J[match(contr$id, ids)] <- comp_J[match(contr$id, ids)] + contr$v
        means_J[[J]][as.character(kk)] <- sum(contr$v) / n
      }
    }

    # Visit-level score: C_k^+(ell_k) * I(Y_k(ell_k)=1) * (Q_{k+1} - M_k)
    if (kk < K - 1L) {
      rows_q <- copy(long[k == kk & Y_end == 1L])
      if (nrow(rows_q)) {
        cur <- state$values[[as.character(kk)]][, list(id = id, M = M)]
        nxt <- state$values[[as.character(kk + 1L)]][, list(id = id, Q_next = Q)]
        dd <- merge(rows_q, cur, by = "id", all.x = TRUE)
        dd <- merge(dd, nxt, by = "id", all.x = TRUE)
        dd <- dd[is.finite(M) & is.finite(Q_next)]
        if (nrow(dd)) {
          v <- dd$Cplus * (dd$Q_next - dd$M)
          total[match(dd$id, ids)] <- total[match(dd$id, ids)] + v
          comp_Q[match(dd$id, ids)] <- comp_Q[match(dd$id, ids)] + v
          means_Q[as.character(kk)] <- sum(v) / n
        }
      }
    }
  }

  sigma <- stats::sd(total)
  # Stopping criterion for the EIF score. The natural scale for asymptotic
  # linearity is the estimate's own standard error sigma/sqrt(n): once every
  # component score is below it, the residual bias from incomplete targeting is
  # dominated by the sampling error. The earlier sigma/(sqrt(n) log n) rule is an
  # order of magnitude tighter and is rarely met by the visit-level ICE score in
  # finite samples even though the estimate is unchanged (see the insensitivity
  # check in scripts/itt_pp/diag_convergence.R / convergence_report.md). The
  # multiplier is overridable via options(ittpp.tol_factor=) for that check.
  tol_factor <- getOption("ittpp.tol_factor", 1.0)
  tol <- tol_factor * sigma / sqrt(n)

  # Aggregate coarsening means
  max_abs_eif_J <- 0
  for (J in Jset) {
    max_abs_eif_J <- max(max_abs_eif_J, max(abs(means_J[[J]])))
  }

  list(
    eif            = total,
    means_T        = means_T,
    means_J        = means_J,
    means_Q        = means_Q,
    max_abs_eif_T  = max(abs(means_T)),
    max_abs_eif_J  = max_abs_eif_J,
    max_abs_eif_Q  = max(abs(means_Q)),
    max_abs_total_eif = abs(mean(total)),
    sigma          = sigma,
    tolerance      = tol,
    # Per-subject EIF component contributions (for the variance decomposition)
    var_total      = stats::var(total),
    var_base       = stats::var(comp_base),
    var_T          = stats::var(comp_T),
    var_J          = stats::var(comp_J),
    var_Q          = stats::var(comp_Q)
  )
}

ittpp_eif_component_objective <- function(diag) {
  max(diag$max_abs_eif_T, diag$max_abs_eif_J, diag$max_abs_eif_Q)
}
