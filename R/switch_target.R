###############################################################################
# ITT/PP Triple Targeting + Backward State
#
# Adapts R/04_va_ct_tmle.R for the ITT/PP estimand.
# Key differences:
#   - Triple targeting: failure + coarsening(s) + ICE
#   - Failure-only forward survival S^T_k (not all-cause)
#   - Coarsening clever covariates have POSITIVE sign
#   - Convergence: max(|P_n D^T|, |P_n D^J|, |P_n D^Q|) <= sigma/(sqrt(n)*log(n))
###############################################################################

ittpp_fit_q_model <- function(rows, q) {
  dd <- copy(rows)
  dd[, Q := va_bound(q)]
  candidates <- c("L_k", "W1", "W2", "k_f")
  dd[, k_f := factor(k)]
  if (nrow(dd) < 3 || stats::sd(dd$Q) < 1e-10) {
    return(list(type = "constant", value = mean(dd$Q)))
  }
  f <- va_build_formula("Q", dd, candidates)
  fit <- tryCatch(lm(f, data = dd), error = function(e) NULL)
  if (is.null(fit)) {
    list(type = "constant", value = mean(dd$Q))
  } else {
    list(type = "lm", fit = fit, value = mean(dd$Q))
  }
}

ittpp_predict_q_model <- function(q_model, rows) {
  if (is.null(q_model)) return(rep(1, nrow(rows)))
  if (identical(q_model$type, "constant")) return(rep(q_model$value, nrow(rows)))
  nd <- copy(rows)
  nd[, k_f := factor(k)]
  pred <- tryCatch(
    as.numeric(predict(q_model$fit, newdata = nd)),
    error = function(e) rep(q_model$value, nrow(rows))
  )
  va_bound(pred)
}

ittpp_mc_next_q <- function(rows, d, nuisance, q_next_model, M = 100) {
  if (!nrow(rows)) return(numeric())
  mu <- ittpp_predict_transition_mean(nuisance, rows, d)
  sig <- nuisance$transition_sigma %||% 0
  n <- nrow(rows)
  acc <- numeric(n)
  for (m in seq_len(M)) {
    nd <- copy(rows)
    nd[, k := k + 1L]
    nd[, L_k := if (sig > 1e-10) mu + rnorm(n, 0, sig) else mu]
    nd[, A0 := as.integer(d)]
    acc <- acc + ittpp_predict_q_model(q_next_model, nd)
  }
  acc / M
}

ittpp_interval_survival_T <- function(rows, d, nuisance, M_k,
                                       hazard_eps = NULL, B = 10) {
  # Failure-only within-interval survival: exp(-integral lambda^T_k)
  # Uses ONLY failure hazard (not all-cause like composite)
  if (!nrow(rows)) return(numeric())
  nd <- ittpp_set_intervention(rows, d)
  lambda_T <- ittpp_predict_hazard(nuisance, nd, cause = "T")
  ell <- nd$ell

  eps_T <- (hazard_eps$T %||% numeric())[as.character(nd$k)]
  eps_T[is.na(eps_T)] <- 0

  # Sub-bin integration for the targeting step
  haz_int <- numeric(nrow(nd))
  Cobs <- if ("Cplus" %in% names(rows)) rows$Cplus else rep(1, nrow(rows))

  # Exact within-interval weight (option a): elapsed-time coarsening survival
  exact <- identical(getOption("ittpp.weight_mode", "projection"), "exact")
  lam_J_sum <- numeric(nrow(nd))
  if (exact) {
    Jset_w <- if (identical(nuisance$estimand, "ITT")) "C" else c("R", "C")
    for (Jw in Jset_w) {
      fitJ <- switch(Jw, R = nuisance$hazard_R, C = nuisance$hazard_C)
      if (is.null(fitJ$model) && fitJ$fallback <= 1e-9) next
      lam_J_sum <- lam_J_sum + ittpp_predict_hazard(nuisance, nd, cause = Jw)
    }
  }

  for (b in seq_len(B)) {
    u <- (b - 0.5) / B * ell
    dt_b <- ell / B
    # Failure-only forward survival
    sbar_T <- exp(-(ell - u) * lambda_T)
    Cu <- if (exact) Cobs * exp(-lam_J_sum * (ell - u)) else Cobs
    # Clever covariate for failure: NEGATIVE sign
    h <- -Cu * sbar_T * M_k
    lt <- lambda_T * exp(eps_T * h)
    haz_int <- haz_int + lt * dt_b
  }
  exp(-haz_int)
}

ittpp_backward_state <- function(obs_data, nuisance, d,
                                  hazard_eps = NULL, visit_eps = NULL,
                                  estimand = c("PP", "ITT"),
                                  M = 100, B = 10, seed = NULL) {
  estimand <- match.arg(estimand)
  if (!is.null(seed)) set.seed(seed)
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)

  # Build cumulative weights
  long <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                          hazard_eps = hazard_eps)

  K <- max(long$k) + 1L
  q_models <- vector("list", K)
  values <- vector("list", K)
  names(values) <- as.character(0:(K - 1L))

  for (kk in rev(0:(K - 1L))) {
    rows <- copy(long[k == kk])
    if (!nrow(rows)) next

    if (kk == K - 1L) {
      M0 <- rep(1, nrow(rows))
    } else {
      M0 <- ittpp_mc_next_q(rows, d, nuisance, q_models[[kk + 2L]], M = M)
    }

    Cpred <- rows$Cplus
    eps_Q <- (visit_eps %||% numeric())[as.character(kk)]
    eps_Q <- ifelse(is.na(eps_Q), 0, eps_Q)
    Mstar <- if (kk == K - 1L) {
      rep(1, nrow(rows))
    } else {
      va_bound(va_expit(va_logit(M0) + eps_Q * Cpred))
    }

    # Failure-only survival (not all-cause!)
    S <- ittpp_interval_survival_T(rows, d, nuisance, Mstar,
                                    hazard_eps = hazard_eps, B = B)
    Q <- va_bound(S * Mstar)
    q_models[[kk + 1L]] <- ittpp_fit_q_model(rows, Q)
    values[[as.character(kk)]] <- data.table(
      id = rows$id, k = kk, S = S, M0 = M0, M = Mstar, Q = Q
    )
  }

  v0 <- values[["0"]]
  list(
    psi        = mean(v0$Q),
    q_models   = q_models,
    values     = values,
    hazard_eps = hazard_eps,
    visit_eps  = visit_eps
  )
}

ittpp_target_failure <- function(obs_data, nuisance, d, state,
                                  estimand = c("PP", "ITT"), B = 10) {
  # Poisson fluctuation of lambda^T with h^T (NEGATIVE sign)
  estimand <- match.arg(estimand)
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  long <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                          hazard_eps = state$hazard_eps)
  K <- max(long$k) + 1L
  eps_T <- state$hazard_eps$T %||% setNames(rep(0, K), as.character(0:(K - 1L)))

  for (kk in 0:(K - 1L)) {
    rows <- long[k == kk]
    sb_T <- ittpp_make_subbins(rows, d, nuisance, state, "T", estimand, B = B)
    if (nrow(sb_T)) {
      eps_T[as.character(kk)] <- eps_T[as.character(kk)] +
        va_solve_poisson_eps(sb_T$h, sb_T$count, sb_T$mu0)
    }
  }
  eps_T
}

ittpp_target_coarsening <- function(obs_data, nuisance, d, state,
                                     estimand = c("PP", "ITT"), B = 10) {
  # Poisson fluctuation of lambda^J with h^J (POSITIVE sign) for each J in Jset
  estimand <- match.arg(estimand)
  Jset <- if (estimand == "PP") c("R", "C") else "C"
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  long <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                          hazard_eps = state$hazard_eps)
  K <- max(long$k) + 1L

  eps_out <- list()
  for (J in Jset) {
    eps_J <- state$hazard_eps[[J]] %||% setNames(rep(0, K), as.character(0:(K - 1L)))
    for (kk in 0:(K - 1L)) {
      rows <- long[k == kk]
      sb_J <- ittpp_make_subbins(rows, d, nuisance, state, J, estimand, B = B)
      if (nrow(sb_J)) {
        eps_J[as.character(kk)] <- eps_J[as.character(kk)] +
          va_solve_poisson_eps(sb_J$h, sb_J$count, sb_J$mu0)
      }
    }
    eps_out[[J]] <- eps_J
  }
  eps_out
}

ittpp_target_visit <- function(obs_data, nuisance, d, state,
                                estimand = c("PP", "ITT")) {
  estimand <- match.arg(estimand)
  long <- ittpp_analysis_long(obs_data, nuisance$tau, estimand)
  long <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                          hazard_eps = state$hazard_eps)
  K <- max(long$k) + 1L
  eps_Q <- state$visit_eps %||% setNames(rep(0, K), as.character(0:(K - 1L)))

  if (K < 2L) return(eps_Q)

  for (kk in rev(0:(K - 2L))) {
    rows <- copy(long[k == kk & Y_end == 1L])
    if (!nrow(rows)) next
    cur <- state$values[[as.character(kk)]][, list(id = id, M0 = M0)]
    nxt <- state$values[[as.character(kk + 1L)]][, list(id = id, Q_next = Q)]
    dd <- merge(rows, cur, by = "id", all.x = TRUE)
    dd <- merge(dd, nxt, by = "id", all.x = TRUE)
    dd <- dd[is.finite(M0) & is.finite(Q_next)]
    if (!nrow(dd)) next

    Cobs <- dd$Cplus
    eps_add <- va_solve_logistic_eps(va_logit(dd$M0), Cobs, dd$Q_next)
    eps_Q[as.character(kk)] <- eps_Q[as.character(kk)] + eps_add
  }
  eps_Q
}

ittpp_blend_named_vector <- function(old, new, alpha) {
  nm <- union(names(old), names(new))
  old2 <- old[nm]; new2 <- new[nm]
  old2[is.na(old2)] <- 0; new2[is.na(new2)] <- 0
  out <- old2 + alpha * (new2 - old2)
  names(out) <- nm
  out
}

ittpp_blend_hazard_eps <- function(old, new, alpha) {
  all_names <- union(names(old), names(new))
  out <- list()
  for (nm in all_names) {
    out[[nm]] <- ittpp_blend_named_vector(
      old[[nm]] %||% numeric(), new[[nm]] %||% numeric(), alpha
    )
  }
  out
}

ittpp_estimate_survival_d <- function(obs_data, d,
                                       estimator = c("ct_tmle", "visit_tmle", "ct_gcomp"),
                                       nuisance = NULL,
                                       estimand = c("PP", "ITT"),
                                       tau = obs_data$tau, g0 = 0.5,
                                       M = 100, B = 10, max_iter = 5,
                                       seed = NULL,
                                       update_method = c("adaptive", "standard"),
                                       step_size = 1, min_step_size = 1e-4,
                                       max_line_iter = 10,
                                       ctmle = TRUE) {
  estimator <- match.arg(estimator)
  estimand <- match.arg(estimand)
  update_method <- match.arg(update_method)

  if (is.null(nuisance)) {
    nuisance <- fit_ittpp_nuisance(obs_data, tau = tau, g0 = g0, estimand = estimand)
  }

  K <- nrow(nuisance$intervals)
  Jset <- if (estimand == "PP") c("R", "C") else "C"

  # Initialize hazard eps for all causes
  init_eps <- function() setNames(rep(0, K), as.character(0:(K - 1L)))
  hazard_eps <- list(T = init_eps())
  for (J in Jset) hazard_eps[[J]] <- init_eps()
  visit_eps <- init_eps()

  # Collaborative TMLE: select coarsening-mechanism complexity that minimizes the
  # targeted outcome residual, before the targeting loop. Only for ct_tmle; this
  # is what distinguishes the proposed estimator from the plug-in ablation.
  ctmle_selected <- NA_character_
  if (isTRUE(ctmle) && estimator == "ct_tmle" &&
      exists("ittpp_ctmle_ladder") && exists("ittpp_ctmle_select")) {
    # Selection uses lighter B/M for affordability; the final targeting loop below
    # runs at full B/M with the selected coarsening mechanism.
    B_sel <- max(3L, min(B, 5L))
    M_sel <- max(10L, min(M, 20L))
    init_state <- ittpp_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps,
                                       estimand, M = M_sel, B = B_sel, seed = seed)
    cand <- tryCatch(
      ittpp_ctmle_ladder(obs_data, nuisance, d, init_state, estimand, B_sel),
      error = function(e) NULL
    )
    if (!is.null(cand)) {
      sel <- tryCatch(
        ittpp_ctmle_select(obs_data, cand, d, init_state, estimand, B_sel, M_sel, seed),
        error = function(e) NULL
      )
      if (!is.null(sel)) {
        nuisance <- sel$nuisance
        ctmle_selected <- sel$name
      }
    }
  }

  state <- ittpp_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps,
                                 estimand, M = M, B = B, seed = seed)
  diag <- ittpp_eif_components(obs_data, nuisance, d, state, estimand, B = B)
  diag0 <- diag           # initial (pre-targeting) diagnostics, for augmentation check
  n_iter <- 0L

  if (estimator != "ct_gcomp") {
    for (iter in seq_len(max_iter)) {
      n_iter <- iter

      if (update_method == "standard") {
        # --- Failure targeting ---
        if (estimator == "ct_tmle") {
          new_eps_T <- ittpp_target_failure(obs_data, nuisance, d, state, estimand, B)
          hazard_eps$T <- new_eps_T

          # --- Coarsening targeting ---
          new_eps_J <- ittpp_target_coarsening(obs_data, nuisance, d, state, estimand, B)
          for (J in Jset) hazard_eps[[J]] <- new_eps_J[[J]]

          state <- ittpp_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps,
                                         estimand, M = M, B = B)
        }

        # --- Visit-level ICE targeting ---
        visit_eps <- ittpp_target_visit(obs_data, nuisance, d, state, estimand)
        state <- ittpp_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps,
                                       estimand, M = M, B = B)
        diag <- ittpp_eif_components(obs_data, nuisance, d, state, estimand, B = B)

      } else {
        # Adaptive step-halving line search
        current_obj <- ittpp_eif_component_objective(diag)
        old_hazard_eps <- hazard_eps
        old_visit_eps  <- visit_eps
        old_state      <- state
        old_diag       <- diag
        accepted       <- FALSE
        alpha          <- step_size

        # Compute full step
        full_hazard_eps <- hazard_eps
        if (estimator == "ct_tmle") {
          full_hazard_eps$T <- ittpp_target_failure(obs_data, nuisance, d, old_state, estimand, B)
          full_eps_J <- ittpp_target_coarsening(obs_data, nuisance, d, old_state, estimand, B)
          for (J in Jset) full_hazard_eps[[J]] <- full_eps_J[[J]]
        }

        for (line_iter in seq_len(max_line_iter)) {
          hazard_try <- if (estimator == "ct_tmle") {
            ittpp_blend_hazard_eps(old_hazard_eps, full_hazard_eps, alpha)
          } else {
            old_hazard_eps
          }
          state_h <- ittpp_backward_state(obs_data, nuisance, d, hazard_try, old_visit_eps,
                                           estimand, M = M, B = B)
          full_visit_eps <- ittpp_target_visit(obs_data, nuisance, d, state_h, estimand)
          visit_try <- ittpp_blend_named_vector(old_visit_eps, full_visit_eps, alpha)
          state_try <- ittpp_backward_state(obs_data, nuisance, d, hazard_try, visit_try,
                                             estimand, M = M, B = B)
          diag_try <- ittpp_eif_components(obs_data, nuisance, d, state_try, estimand, B = B)
          try_obj <- ittpp_eif_component_objective(diag_try)

          if (is.finite(try_obj) && (try_obj < current_obj || try_obj <= diag_try$tolerance)) {
            hazard_eps <- hazard_try
            visit_eps  <- visit_try
            state      <- state_try
            diag       <- diag_try
            accepted   <- TRUE
            break
          }
          alpha <- alpha / 2
          if (alpha < min_step_size) break
        }

        if (!accepted) {
          hazard_eps <- old_hazard_eps
          visit_eps  <- old_visit_eps
          state      <- old_state
          diag       <- old_diag
          break
        }
      }

      if (ittpp_eif_component_objective(diag) <= diag$tolerance) break
    }
  }

  # Weight diagnostics
  long_w <- ittpp_build_cumulative_weights(obs_data, nuisance, d, estimand,
                                            hazard_eps = hazard_eps)
  wdiag <- ittpp_weight_diagnostics(long_w, d)

  list(
    estimator = estimator,
    estimand  = estimand,
    d         = d,
    psi_surv  = va_bound(state$psi),
    risk      = 1 - va_bound(state$psi),
    eif       = diag$eif,
    se_surv   = diag$sigma / sqrt(ittpp_subject_n(obs_data)),
    ci_surv   = va_bound(state$psi) + c(-1, 1) * 1.96 * diag$sigma / sqrt(ittpp_subject_n(obs_data)),
    converged = ittpp_eif_component_objective(diag) <= diag$tolerance,
    n_outer_iterations = n_iter,
    diagnostics = diag,
    initial_diagnostics = diag0,
    weights   = wdiag,
    state     = state,
    nuisance  = nuisance,
    ctmle_selected = ctmle_selected
  )
}

ittpp_estimate_contrast <- function(obs_data,
                                     estimator = c("ct_tmle", "visit_tmle", "ct_gcomp"),
                                     estimand = c("PP", "ITT"),
                                     tau = obs_data$tau, g0 = 0.5,
                                     M = 100, B = 10, max_iter = 5,
                                     seed = NULL,
                                     true_rd = NA_real_,
                                     scenario = NA_character_,
                                     replicate_id = NA_integer_,
                                     update_method = c("adaptive", "standard"),
                                     step_size = 1, min_step_size = 1e-4,
                                     max_line_iter = 10,
                                     ctmle = TRUE,
                                     misspec = NULL) {
  estimator <- match.arg(estimator)
  estimand <- match.arg(estimand)
  update_method <- match.arg(update_method)
  t0 <- proc.time()[["elapsed"]]

  nuisance <- fit_ittpp_nuisance(obs_data, tau = tau, g0 = g0,
                                  estimand = estimand, misspec = misspec)

  fit0 <- ittpp_estimate_survival_d(
    obs_data, 0L, estimator, nuisance, estimand, tau, g0, M, B, max_iter, seed,
    update_method = update_method, step_size = step_size,
    min_step_size = min_step_size, max_line_iter = max_line_iter,
    ctmle = ctmle
  )
  fit1 <- ittpp_estimate_survival_d(
    obs_data, 1L, estimator, nuisance, estimand, tau, g0, M, B, max_iter,
    seed = if (!is.null(seed)) seed + 1L else NULL,
    update_method = update_method, step_size = step_size,
    min_step_size = min_step_size, max_line_iter = max_line_iter,
    ctmle = ctmle
  )

  eif_delta <- fit0$eif - fit1$eif
  se <- stats::sd(eif_delta) / sqrt(ittpp_subject_n(obs_data))
  rd <- fit1$risk - fit0$risk
  ci <- rd + c(-1, 1) * 1.96 * se
  elapsed <- proc.time()[["elapsed"]] - t0

  data.table(
    scenario      = scenario,
    replicate_id  = replicate_id,
    n             = ittpp_subject_n(obs_data),
    tau           = tau,
    estimator     = estimator,
    estimand      = estimand,
    update_method = update_method,
    psi_surv_d0   = fit0$psi_surv,
    psi_surv_d1   = fit1$psi_surv,
    risk_d0       = fit0$risk,
    risk_d1       = fit1$risk,
    risk_difference    = rd,
    true_risk_difference = true_rd,
    bias          = rd - true_rd,
    se            = se,
    ci_lower      = ci[1],
    ci_upper      = ci[2],
    covered       = if (is.finite(true_rd)) ci[1] <= true_rd && true_rd <= ci[2] else NA,
    runtime_seconds = elapsed,
    converged     = fit0$converged && fit1$converged,
    n_outer_iterations = max(fit0$n_outer_iterations, fit1$n_outer_iterations),
    max_abs_eif_T = max(fit0$diagnostics$max_abs_eif_T, fit1$diagnostics$max_abs_eif_T),
    max_abs_eif_J = max(fit0$diagnostics$max_abs_eif_J, fit1$diagnostics$max_abs_eif_J),
    max_abs_eif_Q = max(fit0$diagnostics$max_abs_eif_Q, fit1$diagnostics$max_abs_eif_Q),
    max_abs_total_eif = abs(mean(eif_delta)),
    mean_weight_d0 = fit0$weights$mean_weight,
    mean_weight_d1 = fit1$weights$mean_weight,
    max_weight_d0  = fit0$weights$max_weight,
    max_weight_d1  = fit1$weights$max_weight,
    ess_d0         = fit0$weights$ess,
    ess_d1         = fit1$weights$ess
  )
}
