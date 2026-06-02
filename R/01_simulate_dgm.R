###############################################################################
# Composite-endpoint DGM for the revised VA-CT-TMLE manuscript.
###############################################################################

va_default_params <- function() {
  list(
    visit_times = c(0, 1, 3, 6, 12),
    tau = 12,
    dt = 0.01,
    alpha_0 = 0.5,
    alpha_A = -0.3,
    alpha_1 = 0.2,
    alpha_2 = 0.1,
    rho = 0.1,
    sigma = 0.3,
    sigma_L = 0.15,
    lambda0_T = 0.005,
    beta_A = -0.5,
    beta_L = 0.4,
    beta_W1 = 0.3,
    beta_W2 = 0.2,
    beta_Z = 0,
    lambda0_D = 0.003,
    gamma_A = 0,
    gamma_L = 0.5,
    gamma_W1 = 0.2,
    gamma_W2 = 0,
    gamma_Z = 0,
    scenario = "S1",
    treatment_model = "baseline_only",
    s6_defined = FALSE,
    s6_eta_0 = -0.2,
    s6_eta_prev = 1.0,
    s6_eta_L = -0.5,
    s6_eta_W1 = 0.15,
    s6_eta_W2 = 0.2,
    s6_g_bound = 0.15
  )
}

va_get_scenario <- function(name) {
  p <- va_default_params()
  p$scenario <- name
  if (name == "S1") {
    p$beta_Z <- 0
    p$gamma_Z <- 0
    p$lambda0_D <- 0.003
  } else if (name == "S3") {
    p$beta_Z <- 0
    p$gamma_Z <- 0
    p$lambda0_D <- 0.004
  } else if (name == "S4c") {
    p$beta_Z <- 0.3
    p$gamma_Z <- 0.3
  } else if (name == "S6") {
    p$beta_Z <- 0
    p$gamma_Z <- 0
    p$lambda0_D <- 0.003
    p$gamma_L <- 0.5
    p$treatment_model <- "visit_stochastic"
    p$s6_defined <- TRUE
  } else if (grepl("^GZ_", name)) {
    p$gamma_Z <- as.numeric(sub("^GZ_", "", name))
  } else if (grepl("^JZ_", name)) {
    z <- as.numeric(sub("^JZ_", "", name))
    p$beta_Z <- z
    p$gamma_Z <- z
  } else {
    stop("Unknown scenario: ", name)
  }
  p
}

va_visit_treatment_prob <- function(L_k, A_prev, W1, W2, params) {
  if (!identical(params$treatment_model %||% "baseline_only", "visit_stochastic")) {
    return(as.numeric(A_prev))
  }
  eta <- params$s6_eta_0 +
    params$s6_eta_prev * A_prev +
    params$s6_eta_L * L_k +
    params$s6_eta_W1 * W1 +
    params$s6_eta_W2 * W2
  p <- va_expit(eta)
  b <- params$s6_g_bound %||% 0
  pmin(pmax(p, b), 1 - b)
}

va_static_regime_g_cum <- function(d, W1, W2, L_mat, params) {
  n <- length(W1)
  K <- ncol(L_mat) - 1L
  out <- matrix(NA_real_, n, K)
  if (K < 1L) return(out)
  out[, 1L] <- 0.5
  if (K > 1L) {
    for (kk in 2:K) {
      Lk <- L_mat[, kk]
      p1 <- va_visit_treatment_prob(Lk, rep(as.integer(d), n), W1, W2, params)
      gd <- if (as.integer(d) == 1L) p1 else 1 - p1
      out[, kk] <- out[, kk - 1L] * gd
    }
  }
  out
}

simulate_va_dgm <- function(n, params = va_default_params(), intervention = NULL,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  vt <- p$visit_times
  tau <- p$tau
  dt <- p$dt
  K <- length(vt) - 1L
  t_grid <- seq(0, tau, by = dt)
  visit_idx <- match(round(vt / dt) * dt, round(t_grid / dt) * dt)

  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  A0 <- if (is.null(intervention)) rbinom(n, 1, 0.5) else rep(as.integer(intervention), n)
  A_current <- A0
  A_mat <- matrix(NA_integer_, n, K)
  g_mat <- matrix(NA_real_, n, K)
  A_mat[, 1L] <- A0
  g_mat[, 1L] <- if (is.null(intervention)) 0.5 else 1

  Z <- W1
  L_current <- Z + rnorm(n, 0, p$sigma_L)
  L_mat <- matrix(NA_real_, n, length(vt))
  L_mat[, 1] <- L_current

  U_T <- -log(runif(n))
  U_D <- -log(runif(n))
  H_T <- rep(0, n)
  H_D <- rep(0, n)
  T_time <- rep(Inf, n)
  D_time <- rep(Inf, n)
  alive <- rep(TRUE, n)

  next_visit <- 2L
  for (j in 2:length(t_grid)) {
    t_prev <- t_grid[j - 1L]
    t_now <- t_grid[j]
    idx <- which(alive)
    if (!length(idx)) break

    drift <- (p$alpha_0 + p$alpha_A * A_current[idx] +
      p$alpha_1 * W1[idx] + p$alpha_2 * W2[idx] - p$rho * Z[idx]) * dt
    Z[idx] <- Z[idx] + drift + p$sigma * sqrt(dt) * rnorm(length(idx))

    lambda_T <- p$lambda0_T * exp(
      p$beta_A * A_current[idx] + p$beta_L * L_current[idx] +
        p$beta_Z * Z[idx] + p$beta_W1 * W1[idx] + p$beta_W2 * W2[idx]
    )
    lambda_D <- p$lambda0_D * exp(
      p$gamma_A * A_current[idx] + p$gamma_L * L_current[idx] +
        p$gamma_Z * Z[idx] + p$gamma_W1 * W1[idx] + p$gamma_W2 * W2[idx]
    )
    H_T[idx] <- H_T[idx] + lambda_T * dt
    H_D[idx] <- H_D[idx] + lambda_D * dt

    new_T <- idx[H_T[idx] >= U_T[idx]]
    new_D <- idx[H_D[idx] >= U_D[idx]]
    if (length(new_T) || length(new_D)) {
      event_ids <- union(new_T, new_D)
      for (ii in event_ids) {
        if (!alive[ii]) next
        t_cross <- ii %in% new_T
        d_cross <- ii %in% new_D
        if (t_cross && d_cross) {
          pr_t <- lambda_T[match(ii, idx)] / (lambda_T[match(ii, idx)] + lambda_D[match(ii, idx)])
          if (runif(1) < pr_t) T_time[ii] <- t_now else D_time[ii] <- t_now
        } else if (t_cross) {
          T_time[ii] <- t_now
        } else {
          D_time[ii] <- t_now
        }
        alive[ii] <- FALSE
      }
    }

    if (next_visit <= length(vt) && j == visit_idx[next_visit]) {
      idx2 <- which(alive)
      if (length(idx2)) {
        L_current[idx2] <- Z[idx2] + rnorm(length(idx2), 0, p$sigma_L)
        L_mat[idx2, next_visit] <- L_current[idx2]
        if (next_visit <= K) {
          if (is.null(intervention)) {
            if (identical(p$treatment_model %||% "baseline_only", "visit_stochastic")) {
              p1 <- va_visit_treatment_prob(
                L_current[idx2], A_current[idx2], W1[idx2], W2[idx2], p
              )
              A_next <- rbinom(length(idx2), 1, p1)
              g_next <- ifelse(A_next == 1L, p1, 1 - p1)
            } else {
              A_next <- A_current[idx2]
              g_next <- rep(1, length(idx2))
            }
          } else {
            A_next <- rep(as.integer(intervention), length(idx2))
            g_next <- rep(1, length(idx2))
          }
          A_current[idx2] <- A_next
          A_mat[idx2, next_visit] <- A_next
          g_mat[idx2, next_visit] <- g_next
        }
      }
      next_visit <- next_visit + 1L
    }
  }

  make_long_from_times(
    W1 = W1, W2 = W2, A0 = A0, L_mat = L_mat, T_time = T_time,
    D_time = D_time, visit_times = vt, tau = tau, params = p,
    A_mat = A_mat, g_mat = g_mat
  )
}

simulate_va_truth_subject <- function(n, params = va_default_params(), intervention,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  vt <- p$visit_times
  tau <- p$tau
  dt <- p$dt
  t_grid <- seq(0, tau, by = dt)
  visit_idx <- match(round(vt / dt) * dt, round(t_grid / dt) * dt)

  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  A_current <- rep(as.integer(intervention), n)
  Z <- W1
  L_current <- Z + rnorm(n, 0, p$sigma_L)
  U_T <- -log(runif(n))
  U_D <- -log(runif(n))
  H_T <- rep(0, n)
  H_D <- rep(0, n)
  T_time <- rep(Inf, n)
  D_time <- rep(Inf, n)
  alive <- rep(TRUE, n)

  next_visit <- 2L
  for (j in 2:length(t_grid)) {
    t_now <- t_grid[j]
    idx <- which(alive)
    if (!length(idx)) break
    drift <- (p$alpha_0 + p$alpha_A * A_current[idx] +
      p$alpha_1 * W1[idx] + p$alpha_2 * W2[idx] - p$rho * Z[idx]) * dt
    Z[idx] <- Z[idx] + drift + p$sigma * sqrt(dt) * rnorm(length(idx))
    lambda_T <- p$lambda0_T * exp(
      p$beta_A * A_current[idx] + p$beta_L * L_current[idx] +
        p$beta_Z * Z[idx] + p$beta_W1 * W1[idx] + p$beta_W2 * W2[idx]
    )
    lambda_D <- p$lambda0_D * exp(
      p$gamma_A * A_current[idx] + p$gamma_L * L_current[idx] +
        p$gamma_Z * Z[idx] + p$gamma_W1 * W1[idx] + p$gamma_W2 * W2[idx]
    )
    H_T[idx] <- H_T[idx] + lambda_T * dt
    H_D[idx] <- H_D[idx] + lambda_D * dt
    new_T <- idx[H_T[idx] >= U_T[idx]]
    new_D <- idx[H_D[idx] >= U_D[idx]]
    if (length(new_T) || length(new_D)) {
      event_ids <- union(new_T, new_D)
      for (ii in event_ids) {
        if (!alive[ii]) next
        t_cross <- ii %in% new_T
        d_cross <- ii %in% new_D
        if (t_cross && d_cross) {
          pos <- match(ii, idx)
          pr_t <- lambda_T[pos] / (lambda_T[pos] + lambda_D[pos])
          if (runif(1) < pr_t) T_time[ii] <- t_now else D_time[ii] <- t_now
        } else if (t_cross) {
          T_time[ii] <- t_now
        } else {
          D_time[ii] <- t_now
        }
        alive[ii] <- FALSE
      }
    }
    if (next_visit <= length(vt) && j == visit_idx[next_visit]) {
      idx2 <- which(alive)
      if (length(idx2)) L_current[idx2] <- Z[idx2] + rnorm(length(idx2), 0, p$sigma_L)
      next_visit <- next_visit + 1L
    }
  }
  data.table(T_time = T_time, D_time = D_time)
}

make_long_from_times <- function(W1, W2, A0, L_mat, T_time, D_time,
                                 visit_times, tau, params = NULL,
                                 A_mat = NULL, g_mat = NULL) {
  n <- length(W1)
  K <- length(visit_times) - 1L
  if (is.null(params)) params <- va_default_params()
  if (is.null(A_mat)) {
    A_mat <- matrix(rep(A0, K), nrow = n, ncol = K)
  }
  if (is.null(g_mat)) {
    g_mat <- matrix(1, nrow = n, ncol = K)
    g_mat[, 1L] <- 0.5
  }
  g_cum_obs <- matrix(NA_real_, n, K)
  follow_d0 <- matrix(FALSE, n, K)
  follow_d1 <- matrix(FALSE, n, K)
  g_cum_obs[, 1L] <- ifelse(is.finite(g_mat[, 1L]), g_mat[, 1L], 1)
  follow_d0[, 1L] <- A_mat[, 1L] == 0L
  follow_d1[, 1L] <- A_mat[, 1L] == 1L
  if (K > 1L) {
    for (kk in 2:K) {
      g_cum_obs[, kk] <- g_cum_obs[, kk - 1L] * ifelse(is.finite(g_mat[, kk]), g_mat[, kk], 1)
      follow_d0[, kk] <- follow_d0[, kk - 1L] & A_mat[, kk] == 0L
      follow_d1[, kk] <- follow_d1[, kk - 1L] & A_mat[, kk] == 1L
    }
  }
  g_cum_d0 <- va_static_regime_g_cum(0L, W1, W2, L_mat, params)
  g_cum_d1 <- va_static_regime_g_cum(1L, W1, W2, L_mat, params)
  comp_time <- pmin(T_time, D_time)
  comp_event <- as.integer(comp_time <= tau)
  cause <- ifelse(T_time <= D_time & T_time <= tau, "T",
    ifelse(D_time < T_time & D_time <= tau, "D", "none")
  )
  subj <- data.table(
    id = seq_len(n), W1 = W1, W2 = W2, A0 = A0,
    T_time = T_time, D_time = D_time,
    composite_time = pmin(comp_time, tau),
    composite_event = comp_event,
    cause = cause
  )
  long <- vector("list", n * K)
  pos <- 1L
  for (i in seq_len(n)) {
    for (k0 in 0:(K - 1L)) {
      t0 <- visit_times[k0 + 1L]
      t1 <- visit_times[k0 + 2L]
      if (comp_time[i] <= t0 + 1e-12) next
      Lk <- L_mat[i, k0 + 1L]
      if (!is.finite(Lk)) next
      t_obs <- min(comp_time[i], t1, tau)
      y_end <- as.integer(comp_time[i] > t1)
      ev_T <- as.integer(T_time[i] > t0 && T_time[i] <= t1 && T_time[i] <= D_time[i] && T_time[i] <= tau)
      ev_D <- as.integer(D_time[i] > t0 && D_time[i] <= t1 && D_time[i] < T_time[i] && D_time[i] <= tau)
      long[[pos]] <- data.table(
        id = i, k = k0, t_start = t0, t_end = t1, ell = t1 - t0,
        W1 = W1[i], W2 = W2[i], A0 = A0[i],
        A_prev = if (k0 == 0L) NA_integer_ else A_mat[i, k0],
        A_k = A_mat[i, k0 + 1L],
        g_k = g_mat[i, k0 + 1L],
        g_cum_observed = g_cum_obs[i, k0 + 1L],
        follow_d0 = follow_d0[i, k0 + 1L],
        follow_d1 = follow_d1[i, k0 + 1L],
        g_cum_d0 = g_cum_d0[i, k0 + 1L],
        g_cum_d1 = g_cum_d1[i, k0 + 1L],
        L_k = Lk,
        L_next = if (y_end && k0 < K - 1L) L_mat[i, k0 + 2L] else NA_real_,
        time_at_risk = pmax(t_obs - t0, 0),
        event_T = ev_T, event_D = ev_D,
        event_comp = as.integer(ev_T || ev_D),
        event_time = if (ev_T) T_time[i] - t0 else if (ev_D) D_time[i] - t0 else NA_real_,
        Y_end = y_end
      )
      pos <- pos + 1L
    }
  }
  pi <- rbindlist(long[seq_len(pos - 1L)])
  list(subject = subj, person_interval = pi, visit_L = L_mat,
       visit_A = A_mat, visit_g = g_mat,
       visit_times = visit_times, tau = tau, params = params)
}

simulate_toy_constant <- function(n, visit_times, lambda_T, lambda_D = 0,
                                  treatment = FALSE, lambda_T_by_A = NULL,
                                  lambda_D_by_A = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(visit_times) - 1L
  tau <- max(visit_times)
  lt_vec <- rep(lambda_T, length.out = K)
  ld_vec <- rep(lambda_D, length.out = K)
  W1 <- rep(0, n)
  W2 <- rep(0, n)
  A0 <- if (treatment) rbinom(n, 1, 0.5) else rep(0L, n)
  L_mat <- matrix(0, n, length(visit_times))
  T_time <- rep(Inf, n)
  D_time <- rep(Inf, n)
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      if (is.finite(T_time[i]) || is.finite(D_time[i])) break
      a <- A0[i]
      lt <- if (!is.null(lambda_T_by_A)) lambda_T_by_A[as.character(a)] else lt_vec[k]
      ld <- if (!is.null(lambda_D_by_A)) lambda_D_by_A[as.character(a)] else ld_vec[k]
      rate <- lt + ld
      if (rate <= 0) next
      wait <- rexp(1, rate = rate)
      ell <- visit_times[k + 1L] - visit_times[k]
      if (wait <= ell) {
        if (runif(1) < lt / rate) T_time[i] <- visit_times[k] + wait else D_time[i] <- visit_times[k] + wait
      }
    }
  }
  make_long_from_times(W1, W2, A0, L_mat, T_time, D_time, visit_times, tau)
}

compute_truth_mc <- function(scenarios = c("S1", "S3", "S4c"),
                             n_truth = 100000, seed = 88001) {
  out <- list()
  row <- 1L
  for (sc in scenarios) {
    p <- va_get_scenario(sc)
    if (identical(sc, "S6") && !isTRUE(p$s6_defined)) {
      out[[row]] <- data.table(
        scenario = sc, n_truth = 0L, psi_surv_d0 = NA_real_,
        psi_surv_d1 = NA_real_, risk_d0 = NA_real_, risk_d1 = NA_real_,
        true_risk_difference = NA_real_,
        note = "Skipped: no stochastic post-baseline treatment/deviation mechanism is defined."
      )
      row <- row + 1L
      next
    }
    dat0 <- simulate_va_truth_subject(n_truth, p, intervention = 0L, seed = seed + row)
    dat1 <- simulate_va_truth_subject(n_truth, p, intervention = 1L, seed = seed + row + 1000L)
    s0 <- mean(dat0$T_time > p$tau & dat0$D_time > p$tau)
    s1 <- mean(dat1$T_time > p$tau & dat1$D_time > p$tau)
    out[[row]] <- data.table(
      scenario = sc, n_truth = n_truth,
      psi_surv_d0 = s0, psi_surv_d1 = s1,
      risk_d0 = 1 - s0, risk_d1 = 1 - s1,
      true_risk_difference = (1 - s1) - (1 - s0),
      note = if (n_truth >= 1000000) "final truth" else "pilot truth"
    )
    row <- row + 1L
  }
  rbindlist(out, fill = TRUE)
}
