###############################################################################
# VA-CT-TMLE core implementation for the composite-survival target.
###############################################################################

va_fit_q_model <- function(rows, q) {
  dd <- copy(rows)
  dd[, Q := va_bound(q)]
  candidates <- c("L_k", "W1", "W2", "k_f")
  dd[, k_f := factor(k)]
  if (nrow(dd) < 3 || stats::sd(dd$Q) < 1e-10) {
    return(list(type = "constant", value = mean(dd$Q)))
  }
  f <- va_build_formula("Q", dd, candidates)
  fit <- tryCatch(lm(f, data = dd), error = function(e) NULL)
  if (is.null(fit)) list(type = "constant", value = mean(dd$Q)) else list(type = "lm", fit = fit, value = mean(dd$Q))
}

va_predict_q_model <- function(q_model, rows) {
  if (is.null(q_model)) return(rep(1, nrow(rows)))
  if (identical(q_model$type, "constant")) return(rep(q_model$value, nrow(rows)))
  nd <- copy(rows)
  nd[, k_f := factor(k)]
  pred <- tryCatch(as.numeric(predict(q_model$fit, newdata = nd)),
                   error = function(e) rep(q_model$value, nrow(rows)))
  va_bound(pred)
}

va_mc_next_q <- function(rows, d, nuisance, q_next_model, M = 100) {
  if (!nrow(rows)) return(numeric())
  mu <- va_predict_transition_mean(nuisance, rows, d)
  sig <- nuisance$transition_sigma %||% 0
  n <- nrow(rows)
  acc <- numeric(n)
  for (m in seq_len(M)) {
    nd <- copy(rows)
    nd[, k := k + 1L]
    nd[, L_k := if (sig > 1e-10) mu + rnorm(n, 0, sig) else mu]
    nd[, A_k := as.integer(d)]
    acc <- acc + va_predict_q_model(q_next_model, nd)
  }
  acc / M
}

va_interval_survival <- function(rows, d, nuisance, M_k, hazard_eps = NULL,
                                 B = 10, observed_A = FALSE) {
  if (!nrow(rows)) return(numeric())
  nd <- if (observed_A) copy(rows) else va_set_intervention(rows, d)
  base <- va_predict_hazard_pair(nuisance, nd, d = NULL)
  ell <- nd$ell
  eps_T <- (hazard_eps$T %||% numeric())[as.character(nd$k)]
  eps_D <- (hazard_eps$D %||% numeric())[as.character(nd$k)]
  eps_T[is.na(eps_T)] <- 0
  eps_D[is.na(eps_D)] <- 0
  Cpred <- if (observed_A) va_cplus_observed(rows, d, nuisance) else va_cplus_intervention(rows, d, nuisance)

  haz_int <- numeric(nrow(nd))
  for (b in seq_len(B)) {
    u <- (b - 0.5) / B * ell
    dt <- ell / B
    sbar0 <- exp(-(ell - u) * (base$lambda_T + base$lambda_D))
    h <- -Cpred * sbar0 * M_k
    lt <- base$lambda_T * exp(eps_T * h)
    ld <- base$lambda_D * exp(eps_D * h)
    haz_int <- haz_int + (lt + ld) * dt
  }
  exp(-haz_int)
}

va_backward_state <- function(obs_data, nuisance, d, hazard_eps = NULL,
                              visit_eps = NULL, M = 100, B = 10,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  long <- va_analysis_long(obs_data, nuisance$tau)
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
      M0 <- va_mc_next_q(rows, d, nuisance, q_models[[kk + 2L]], M = M)
    }
    Cpred <- va_cplus_intervention(rows, d, nuisance)
    eps_Q <- (visit_eps %||% numeric())[as.character(kk)]
    eps_Q <- ifelse(is.na(eps_Q), 0, eps_Q)
    Mstar <- if (kk == K - 1L) rep(1, nrow(rows)) else va_bound(va_expit(va_logit(M0) + eps_Q * Cpred))
    S <- va_interval_survival(rows, d, nuisance, Mstar, hazard_eps = hazard_eps, B = B)
    Q <- va_bound(S * Mstar)
    q_models[[kk + 1L]] <- va_fit_q_model(rows, Q)
    values[[as.character(kk)]] <- data.table(
      id = rows$id, k = kk, S = S, M0 = M0, M = Mstar, Q = Q
    )
  }
  v0 <- values[["0"]]
  list(
    psi = mean(v0$Q),
    q_models = q_models,
    values = values,
    hazard_eps = hazard_eps,
    visit_eps = visit_eps
  )
}

va_make_subbins <- function(rows, d, nuisance, state, cause = c("T", "D"),
                            B = 10, use_updated = FALSE) {
  cause <- match.arg(cause)
  if (!nrow(rows)) return(data.table())
  val <- state$values[[as.character(rows$k[1])]]
  rows <- merge(copy(rows), val[, .(id, M)], by = "id", all.x = TRUE, sort = FALSE)
  base_d <- va_predict_hazard_pair(nuisance, va_set_intervention(rows, d))
  base_obs <- va_predict_hazard_pair(nuisance, rows)
  Cobs <- va_cplus_observed(rows, d, nuisance)
  n <- nrow(rows)
  idx <- rep(seq_len(n), each = B)
  b <- rep(seq_len(B), times = n)
  ell <- rows$ell[idx]
  lo <- (b - 1) * ell / B
  hi <- b * ell / B
  mid <- (lo + hi) / 2
  pt <- pmax(0, pmin(rows$time_at_risk[idx], hi) - lo)
  keep <- pt > 0
  if (!any(keep)) return(data.table())
  idx <- idx[keep]
  lo <- lo[keep]
  hi <- hi[keep]
  mid <- mid[keep]
  ell <- ell[keep]
  pt <- pt[keep]
  sbar <- exp(-(ell - mid) * (base_d$lambda_T[idx] + base_d$lambda_D[idx]))
  h <- -Cobs[idx] * sbar * rows$M[idx]
  event_here <- !is.na(rows$event_time[idx]) &
    rows$event_time[idx] > lo - 1e-12 & rows$event_time[idx] <= hi + 1e-12
  count <- if (cause == "T") {
    as.integer(rows$event_T[idx] == 1L & event_here)
  } else {
    as.integer(rows$event_D[idx] == 1L & event_here)
  }
  lambda0 <- if (cause == "T") base_obs$lambda_T[idx] else base_obs$lambda_D[idx]
  eps_vec <- if (cause == "T") state$hazard_eps$T else state$hazard_eps$D
  eps <- (eps_vec %||% numeric())[as.character(rows$k[idx])]
  eps[is.na(eps)] <- 0
  mu0 <- lambda0 * pt
  mu <- mu0 * exp(if (use_updated) eps * h else 0)
  data.table(id = rows$id[idx], k = rows$k[idx], h = h, count = count, mu0 = mu0, mu = mu)
}

va_target_hazards <- function(obs_data, nuisance, d, state, B = 10) {
  long <- va_analysis_long(obs_data, nuisance$tau)
  K <- max(long$k) + 1L
  eps_T <- state$hazard_eps$T %||% setNames(rep(0, K), as.character(0:(K - 1L)))
  eps_D <- state$hazard_eps$D %||% setNames(rep(0, K), as.character(0:(K - 1L)))
  for (kk in 0:(K - 1L)) {
    rows <- long[k == kk]
    sb_T <- va_make_subbins(rows, d, nuisance, state, "T", B = B)
    sb_D <- va_make_subbins(rows, d, nuisance, state, "D", B = B)
    eps_T[as.character(kk)] <- eps_T[as.character(kk)] + va_solve_poisson_eps(sb_T$h, sb_T$count, sb_T$mu0)
    eps_D[as.character(kk)] <- eps_D[as.character(kk)] + va_solve_poisson_eps(sb_D$h, sb_D$count, sb_D$mu0)
  }
  list(T = eps_T, D = eps_D)
}

va_target_visit <- function(obs_data, nuisance, d, state) {
  long <- va_analysis_long(obs_data, nuisance$tau)
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
    Cobs <- va_cplus_observed(dd, d, nuisance)
    eps_add <- va_solve_logistic_eps(va_logit(dd$M0), Cobs, dd$Q_next)
    eps_Q[as.character(kk)] <- eps_Q[as.character(kk)] + eps_add
  }
  eps_Q
}

va_eif_components <- function(obs_data, nuisance, d, state, B = 10) {
  long <- va_analysis_long(obs_data, nuisance$tau)
  n <- va_subject_n(obs_data)
  ids <- obs_data$subject$id
  total <- numeric(n)
  q0 <- state$values[["0"]][match(ids, id), Q]
  q0[is.na(q0)] <- mean(state$values[["0"]]$Q)
  psi <- mean(q0)
  total <- total + q0 - psi

  K <- max(long$k) + 1L
  means_T <- means_D <- means_Q <- setNames(rep(0, K), as.character(0:(K - 1L)))
  for (kk in 0:(K - 1L)) {
    rows <- long[k == kk]
    sb_T <- va_make_subbins(rows, d, nuisance, state, "T", B = B, use_updated = TRUE)
    sb_D <- va_make_subbins(rows, d, nuisance, state, "D", B = B, use_updated = TRUE)
    if (nrow(sb_T)) {
      contr <- sb_T[, .(v = sum(h * (count - mu))), by = id]
      total[match(contr$id, ids)] <- total[match(contr$id, ids)] + contr$v
      means_T[as.character(kk)] <- sum(contr$v) / n
    }
    if (nrow(sb_D)) {
      contr <- sb_D[, .(v = sum(h * (count - mu))), by = id]
      total[match(contr$id, ids)] <- total[match(contr$id, ids)] + contr$v
      means_D[as.character(kk)] <- sum(contr$v) / n
    }
    if (kk < K - 1L) {
      rows_q <- copy(long[k == kk & Y_end == 1L])
      if (nrow(rows_q)) {
        cur <- state$values[[as.character(kk)]][, list(id = id, M = M)]
        nxt <- state$values[[as.character(kk + 1L)]][, list(id = id, Q_next = Q)]
        dd <- merge(rows_q, cur, by = "id", all.x = TRUE)
        dd <- merge(dd, nxt, by = "id", all.x = TRUE)
        dd <- dd[is.finite(M) & is.finite(Q_next)]
        if (nrow(dd)) {
          v <- va_cplus_observed(dd, d, nuisance) * (dd$Q_next - dd$M)
          total[match(dd$id, ids)] <- total[match(dd$id, ids)] + v
          means_Q[as.character(kk)] <- sum(v) / n
        }
      }
    }
  }
  sigma <- stats::sd(total)
  tol <- sigma / (sqrt(n) * log(n))
  list(
    eif = total,
    means_T = means_T,
    means_D = means_D,
    means_Q = means_Q,
    max_abs_eif_T = max(abs(means_T)),
    max_abs_eif_D = max(abs(means_D)),
    max_abs_eif_Q = max(abs(means_Q)),
    max_abs_total_eif = abs(mean(total)),
    sigma = sigma,
    tolerance = tol
  )
}

va_eif_component_objective <- function(diag) {
  max(diag$max_abs_eif_T, diag$max_abs_eif_D, diag$max_abs_eif_Q)
}

va_trace_row <- function(iteration, status, diag, alpha = NA_real_, line_iter = NA_integer_) {
  objective <- va_eif_component_objective(diag)
  data.table(
    iteration = as.integer(iteration),
    line_iter = as.integer(line_iter),
    status = status,
    alpha = alpha,
    max_abs_eif_T = diag$max_abs_eif_T,
    max_abs_eif_D = diag$max_abs_eif_D,
    max_abs_eif_Q = diag$max_abs_eif_Q,
    max_abs_total_eif = diag$max_abs_total_eif,
    sigma = diag$sigma,
    tolerance = diag$tolerance,
    objective = objective,
    objective_ratio = if (diag$tolerance > 0) objective / diag$tolerance else Inf,
    converged = objective <= diag$tolerance
  )
}

va_blend_named_vector <- function(old, new, alpha) {
  nm <- union(names(old), names(new))
  old2 <- old[nm]
  new2 <- new[nm]
  old2[is.na(old2)] <- 0
  new2[is.na(new2)] <- 0
  out <- old2 + alpha * (new2 - old2)
  names(out) <- nm
  out
}

va_blend_hazard_eps <- function(old, new, alpha) {
  list(
    T = va_blend_named_vector(old$T %||% numeric(), new$T %||% numeric(), alpha),
    D = va_blend_named_vector(old$D %||% numeric(), new$D %||% numeric(), alpha)
  )
}

va_estimate_survival_d <- function(obs_data, d, estimator = c("ct_tmle", "visit_tmle", "ct_gcomp"),
                                   nuisance = NULL, tau = obs_data$tau, g0 = 0.5,
                                   M = 100, B = 10, max_iter = 5,
                                   seed = NULL, update_method = c("standard", "adaptive"),
                                   step_size = 1, min_step_size = 1e-4,
                                   max_line_iter = 10) {
  estimator <- match.arg(estimator)
  update_method <- match.arg(update_method)
  if (is.null(nuisance)) nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  K <- nrow(nuisance$intervals)
  hazard_eps <- list(
    T = setNames(rep(0, K), as.character(0:(K - 1L))),
    D = setNames(rep(0, K), as.character(0:(K - 1L)))
  )
  visit_eps <- setNames(rep(0, K), as.character(0:(K - 1L)))
  state <- va_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps, M = M, B = B, seed = seed)
  diag <- va_eif_components(obs_data, nuisance, d, state, B = B)
  n_iter <- 0L
  trace <- va_trace_row(0L, "initial", diag)

  if (estimator != "ct_gcomp") {
    for (iter in seq_len(max_iter)) {
      n_iter <- iter
      if (update_method == "standard") {
        if (estimator == "ct_tmle") {
          hazard_eps <- va_target_hazards(obs_data, nuisance, d, state, B = B)
          state <- va_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps, M = M, B = B)
        }
        visit_eps <- va_target_visit(obs_data, nuisance, d, state)
        state <- va_backward_state(obs_data, nuisance, d, hazard_eps, visit_eps, M = M, B = B)
        diag <- va_eif_components(obs_data, nuisance, d, state, B = B)
        trace <- rbindlist(list(trace, va_trace_row(iter, "updated", diag, alpha = 1)), fill = TRUE)
      } else {
        current_ratio <- va_trace_row(iter - 1L, "current", diag)$objective_ratio
        old_hazard_eps <- hazard_eps
        old_visit_eps <- visit_eps
        old_state <- state
        old_diag <- diag
        accepted <- FALSE
        alpha <- step_size

        full_hazard_eps <- if (estimator == "ct_tmle") {
          va_target_hazards(obs_data, nuisance, d, old_state, B = B)
        } else {
          old_hazard_eps
        }

        for (line_iter in seq_len(max_line_iter)) {
          hazard_try <- if (estimator == "ct_tmle") {
            va_blend_hazard_eps(old_hazard_eps, full_hazard_eps, alpha)
          } else {
            old_hazard_eps
          }
          state_h <- va_backward_state(obs_data, nuisance, d, hazard_try, old_visit_eps, M = M, B = B)
          full_visit_eps <- va_target_visit(obs_data, nuisance, d, state_h)
          visit_try <- va_blend_named_vector(old_visit_eps, full_visit_eps, alpha)
          state_try <- va_backward_state(obs_data, nuisance, d, hazard_try, visit_try, M = M, B = B)
          diag_try <- va_eif_components(obs_data, nuisance, d, state_try, B = B)
          try_ratio <- va_trace_row(iter, "candidate", diag_try)$objective_ratio

          if (is.finite(try_ratio) && (try_ratio < current_ratio || va_eif_component_objective(diag_try) <= diag_try$tolerance)) {
            hazard_eps <- hazard_try
            visit_eps <- visit_try
            state <- state_try
            diag <- diag_try
            trace <- rbindlist(list(trace, va_trace_row(iter, "accepted", diag, alpha = alpha, line_iter = line_iter)), fill = TRUE)
            accepted <- TRUE
            break
          }

          trace <- rbindlist(list(
            trace,
            va_trace_row(iter, "rejected_increased_objective", diag_try, alpha = alpha, line_iter = line_iter)
          ), fill = TRUE)
          alpha <- alpha / 2
          if (alpha < min_step_size) break
        }

        if (!accepted) {
          hazard_eps <- old_hazard_eps
          visit_eps <- old_visit_eps
          state <- old_state
          diag <- old_diag
          trace <- rbindlist(list(trace, va_trace_row(iter, "no_accepted_step", diag, alpha = alpha)), fill = TRUE)
          break
        }
      }
      if (va_eif_component_objective(diag) <= diag$tolerance) break
    }
  }

  rows <- va_analysis_long(obs_data, tau)
  wdiag <- va_weight_diagnostics(rows, d, nuisance)
  list(
    estimator = estimator,
    d = d,
    psi_surv = va_bound(state$psi),
    risk = 1 - va_bound(state$psi),
    eif = diag$eif,
    se_surv = diag$sigma / sqrt(va_subject_n(obs_data)),
    ci_surv = va_bound(state$psi) + c(-1, 1) * 1.96 * diag$sigma / sqrt(va_subject_n(obs_data)),
    converged = max(diag$max_abs_eif_T, diag$max_abs_eif_D, diag$max_abs_eif_Q) <= diag$tolerance,
    n_outer_iterations = n_iter,
    diagnostics = diag,
    trace = trace,
    weights = wdiag,
    state = state,
    nuisance = nuisance
  )
}

va_estimate_contrast <- function(obs_data, estimator = c(
                                   "ct_tmle", "visit_tmle", "ct_gcomp",
                                   "va_discrete_gcomp", "iptw_km", "naive_cox"
                                 ),
                                 tau = obs_data$tau, g0 = 0.5, M = 100, B = 10,
                                 max_iter = 5, seed = NULL, true_rd = NA_real_,
                                 scenario = NA_character_, replicate_id = NA_integer_,
                                 update_method = c("standard", "adaptive"),
                                 step_size = 1, min_step_size = 1e-4,
                                 max_line_iter = 10) {
  estimator <- match.arg(estimator)
  update_method <- match.arg(update_method)
  if (estimator %in% va_direct_estimators()) {
    return(va_estimate_direct_comparator_contrast(
      obs_data = obs_data, estimator = estimator, tau = tau, g0 = g0, M = M,
      seed = seed, true_rd = true_rd, scenario = scenario,
      replicate_id = replicate_id
    ))
  }
  t0 <- proc.time()[["elapsed"]]
  nuisance <- fit_va_nuisance(obs_data, tau = tau, g0 = g0)
  fit0 <- va_estimate_survival_d(
    obs_data, 0L, estimator, nuisance, tau, g0, M, B, max_iter, seed,
    update_method = update_method, step_size = step_size,
    min_step_size = min_step_size, max_line_iter = max_line_iter
  )
  fit1 <- va_estimate_survival_d(
    obs_data, 1L, estimator, nuisance, tau, g0, M, B, max_iter, seed + 1L,
    update_method = update_method, step_size = step_size,
    min_step_size = min_step_size, max_line_iter = max_line_iter
  )
  eif_delta <- fit0$eif - fit1$eif
  se <- stats::sd(eif_delta) / sqrt(va_subject_n(obs_data))
  rd <- fit1$risk - fit0$risk
  ci <- rd + c(-1, 1) * 1.96 * se
  elapsed <- proc.time()[["elapsed"]] - t0
  data.table(
    scenario = scenario,
    replicate_id = replicate_id,
    n = va_subject_n(obs_data),
    tau = tau,
    estimator = estimator,
    update_method = update_method,
    psi_surv_d0 = fit0$psi_surv,
    psi_surv_d1 = fit1$psi_surv,
    risk_d0 = fit0$risk,
    risk_d1 = fit1$risk,
    risk_difference = rd,
    true_risk_difference = true_rd,
    bias = rd - true_rd,
    se = se,
    ci_lower = ci[1],
    ci_upper = ci[2],
    covered = if (is.finite(true_rd)) ci[1] <= true_rd && true_rd <= ci[2] else NA,
    runtime_seconds = elapsed,
    converged = fit0$converged && fit1$converged,
    n_outer_iterations = max(fit0$n_outer_iterations, fit1$n_outer_iterations),
    max_abs_eif_T = max(fit0$diagnostics$max_abs_eif_T, fit1$diagnostics$max_abs_eif_T),
    max_abs_eif_D = max(fit0$diagnostics$max_abs_eif_D, fit1$diagnostics$max_abs_eif_D),
    max_abs_eif_Q = max(fit0$diagnostics$max_abs_eif_Q, fit1$diagnostics$max_abs_eif_Q),
    max_abs_total_eif = abs(mean(eif_delta)),
    mean_weight_d0 = fit0$weights$mean_weight,
    mean_weight_d1 = fit1$weights$mean_weight,
    max_weight_d0 = fit0$weights$max_weight,
    max_weight_d1 = fit1$weights$max_weight,
    ess_d0 = fit0$weights$ess,
    ess_d1 = fit1$weights$ess,
    truncation_rate_d0 = fit0$weights$truncation_rate,
    truncation_rate_d1 = fit1$weights$truncation_rate
  )
}
