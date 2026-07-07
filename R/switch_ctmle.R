###############################################################################
# ITT/PP Collaborative TMLE (C-TMLE)
#
# Collaborative selection of coarsening mechanism complexity.
# Builds candidate fits of increasing complexity for the coarsening
# intensities and selects the one minimizing the targeted outcome residual.
###############################################################################

ittpp_ctmle_candidate_formulas <- function() {
  # Ordered by increasing complexity
  # Three-rung ladder of increasing complexity. The middle and full rungs span the
  # bias--variance tradeoff in the coarsening mechanism; the intercept rung is the
  # maximally-stabilized weight. Three rungs keep collaborative selection affordable
  # at simulation scale without losing the qualitative C-TMLE behavior.
  list(
    list(name = "intercept_only", candidates = c("k_f")),
    list(name = "treatment_L",   candidates = c("A0", "L_k", "k_f")),
    list(name = "full",          candidates = c("A0", "L_k", "W1", "W2", "k_f"))
  )
}

ittpp_ctmle_ladder <- function(obs_data, nuisance, d, state,
                                estimand = c("PP", "ITT"), B = 10) {
  estimand <- match.arg(estimand)
  Jset <- if (estimand == "PP") c("R", "C") else "C"
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  formulas <- ittpp_ctmle_candidate_formulas()

  candidates <- list()
  for (idx in seq_along(formulas)) {
    f <- formulas[[idx]]
    nuis_candidate <- nuisance

    # Re-fit coarsening hazards at this complexity level
    for (J in Jset) {
      event_col <- paste0("event_", J)
      if (any(long[[event_col]] == 1L, na.rm = TRUE)) {
        nuis_candidate[[paste0("hazard_", J)]] <-
          va_fit_poisson_rate(long, event_col, f$candidates)
      }
    }

    candidates[[idx]] <- list(
      name      = f$name,
      nuisance  = nuis_candidate,
      complexity = idx
    )
  }
  candidates
}

ittpp_ctmle_select <- function(obs_data, candidates, d, state,
                                estimand = c("PP", "ITT"), B = 10,
                                M = 100, seed = NULL) {
  estimand <- match.arg(estimand)
  K <- max(ittpp_analysis_long(obs_data, candidates[[1]]$nuisance$tau, estimand)$k) + 1L
  Jset <- if (estimand == "PP") c("R", "C") else "C"

  # Collaborative criterion: among coarsening-mechanism candidates, choose the one
  # that minimizes the *targeted outcome residual* (the empirical failure + ICE
  # score, max(|P_n D^T|, |P_n D^Q|)) -- NOT the coarsening score, which a more
  # complex coarsening model would trivially shrink. A parsimony margin breaks ties
  # toward the simpler (weight-stabilizing) candidate, which is the point of C-TMLE.
  margin <- 1.0  # accept a simpler candidate within this multiple of the best residual

  results <- vector("list", length(candidates))
  for (idx in seq_along(candidates)) {
    nuis <- candidates[[idx]]$nuisance
    hazard_eps <- list(T = setNames(rep(0, K), as.character(0:(K - 1L))))
    for (J in Jset) hazard_eps[[J]] <- setNames(rep(0, K), as.character(0:(K - 1L)))
    visit_eps <- setNames(rep(0, K), as.character(0:(K - 1L)))

    trial_state <- tryCatch(
      ittpp_backward_state(obs_data, nuis, d, hazard_eps, visit_eps,
                            estimand, M = M, B = B, seed = seed),
      error = function(e) NULL
    )
    if (is.null(trial_state)) next
    trial_diag <- tryCatch(
      ittpp_eif_components(obs_data, nuis, d, trial_state, estimand, B = B),
      error = function(e) NULL
    )
    if (is.null(trial_diag)) next

    # Outcome residual = failure + ICE only
    outcome_resid <- max(trial_diag$max_abs_eif_T, trial_diag$max_abs_eif_Q)
    if (!is.finite(outcome_resid)) next

    # Weight stress of this candidate (max cumulative weight under d)
    lw <- ittpp_build_cumulative_weights(obs_data, nuis, d, estimand)
    max_w <- suppressWarnings(max(lw$Cplus[lw$Cplus > 0], na.rm = TRUE))

    results[[idx]] <- list(idx = idx, resid = outcome_resid,
                           max_w = if (is.finite(max_w)) max_w else Inf,
                           complexity = idx)
  }

  results <- Filter(Negate(is.null), results)
  if (!length(results)) return(candidates[[length(candidates)]])

  resids <- vapply(results, function(r) r$resid, numeric(1))
  best_resid <- min(resids)
  # Candidates whose outcome residual is within the margin of the best; among those,
  # take the one with the smallest weight stress (simplest sufficient mechanism).
  eligible <- Filter(function(r) r$resid <= best_resid * (1 + margin) + 1e-12, results)
  best <- eligible[[which.min(vapply(eligible, function(r) r$max_w, numeric(1)))]]

  candidates[[best$idx]]
}

ittpp_estimate_survival_d_ctmle <- function(obs_data, d,
                                             estimand = c("PP", "ITT"),
                                             tau = obs_data$tau, g0 = 0.5,
                                             M = 100, B = 10, max_iter = 5,
                                             seed = NULL,
                                             update_method = "adaptive",
                                             step_size = 1,
                                             min_step_size = 1e-4,
                                             max_line_iter = 10) {
  estimand <- match.arg(estimand)

  # Initial nuisance with full complexity
  nuisance <- fit_ittpp_nuisance(obs_data, tau = tau, g0 = g0, estimand = estimand)

  # Build C-TMLE ladder
  K <- nrow(nuisance$intervals)
  hazard_eps <- list(T = setNames(rep(0, K), as.character(0:(K - 1L))))
  Jset <- if (estimand == "PP") c("R", "C") else "C"
  for (J in Jset) hazard_eps[[J]] <- setNames(rep(0, K), as.character(0:(K - 1L)))
  visit_eps <- setNames(rep(0, K), as.character(0:(K - 1L)))

  state <- ittpp_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps,
                                 estimand, M = M, B = B, seed = seed)

  # Build and select C-TMLE candidate
  candidates <- ittpp_ctmle_ladder(obs_data, nuisance, d, state, estimand, B)
  selected <- ittpp_ctmle_select(obs_data, candidates, d, state, estimand, B, M, seed)

  # Run the full estimator with selected nuisance
  ittpp_estimate_survival_d(
    obs_data, d, estimator = "ct_tmle",
    nuisance = selected$nuisance, estimand = estimand,
    tau = tau, g0 = g0, M = M, B = B, max_iter = max_iter,
    seed = seed, update_method = update_method,
    step_size = step_size, min_step_size = min_step_size,
    max_line_iter = max_line_iter, ctmle = FALSE
  )
}
