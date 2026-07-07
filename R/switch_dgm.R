###############################################################################
# ITT/PP Data-Generating Mechanism
#
# Three cause-specific hazards: failure T, switch R, censoring C.
# Adapts the composite DGM (R/01_simulate_dgm.R) for the ITT/PP estimand.
###############################################################################

ittpp_default_params <- function() {
  list(
    visit_times = c(0, 1, 3, 6, 12),
    tau         = 12,
    dt          = 0.01,
    # OU latent process
    alpha_0  = 0.5,
    alpha_A  = -0.3,
    alpha_1  = 0.2,
    alpha_2  = 0.1,
    rho      = 0.1,
    sigma    = 0.3,
    sigma_L  = 0.15,
    # Failure hazard: lambda^T_0 * exp(beta . X)
    lambda0_T = 0.012,
    beta_A    = -0.5,
    beta_L    = 0.4,
    beta_W1   = 0.3,
    beta_W2   = 0.2,
    beta_Z    = 0,
    # Switch hazard: lambda^R_0 * exp(eta . X)
    lambda0_R = 0.003,
    eta_A     = 0.4,
    eta_L     = 0.3,
    eta_W1    = 0.1,
    eta_W2    = 0.1,
    eta_Z     = 0,
    # Censoring hazard: lambda^C_0 * exp(kappa . X)
    lambda0_C = 0.002,
    kappa_A   = 0,
    kappa_L   = 0.3,
    kappa_W1  = 0.1,
    kappa_W2  = 0.1,
    kappa_Z   = 0,
    # Scenario label
    scenario  = "A0"
  )
}

ittpp_get_scenario <- function(name) {
  p <- ittpp_default_params()
  p$scenario <- name

  if (name == "A0") {
    # Reduction check: no switch, no censoring
    p$lambda0_R <- 0
    p$lambda0_C <- 0
  } else if (name == "A1") {
    # Informative censoring only (ITT)
    p$lambda0_R <- 0
    p$lambda0_C <- 0.006
    p$kappa_L   <- 0.5
  } else if (name == "A2") {
    # Crossover only (PP)
    p$lambda0_C <- 0
    p$lambda0_R <- 0.008
    p$eta_L     <- 0.5
  } else if (name == "A3") {
    # Both present, DR matrix — mechanism side tests
    p$lambda0_R <- 0.006
    p$eta_L     <- 0.4
    p$lambda0_C <- 0.005
    p$kappa_L   <- 0.4
  } else if (name == "A4") {
    # Both present, DR matrix — Q-side tests
    p$lambda0_R <- 0.006
    p$eta_L     <- 0.4
    p$lambda0_C <- 0.005
    p$kappa_L   <- 0.4
  } else if (name == "A5") {
    # Positivity stress + C-TMLE
    p$lambda0_R <- 0.008
    p$eta_L     <- 0.8
    p$eta_W1    <- 0.3
    p$lambda0_C <- 0.006
    p$kappa_L   <- 0.7
    p$kappa_W1  <- 0.3
  } else if (name == "A6") {
    # Visit-sufficiency violation: latent Z drives switching/censoring
    p$lambda0_R <- 0.006
    p$eta_L     <- 0.3
    p$eta_Z     <- 0.5
    p$lambda0_C <- 0.005
    p$kappa_L   <- 0.3
    p$kappa_Z   <- 0.4
    # sigma_L is swept externally in {0.15, 0.5, 1.0}
  } else if (name == "A7") {
    # Co-occurrence stress (oncology crossover): high, L-driven crossover, a
    # strong treatment effect on failure (so the arm that "owns" a same-bin
    # death matters), and rapid post-crossover failure -- switches cluster near
    # failures, so a coarse discrete grid routinely bins switch + death together.
    # Motivated by PROFILE-1014/1007 (70-89% crossover) + short post-progression
    # survival (~4.7 mo). Used by the ordering-sensitivity experiment.
    p$lambda0_R <- 0.030      # high crossover intensity
    p$eta_L     <- 0.6        # L-driven switching (time-varying confounding)
    p$beta_A    <- -1.0       # strong protective on-treatment effect
    p$beta_L    <- 0.5
    p$lambda0_T <- 0.020      # higher event rate -> more co-occurrence
    p$lambda0_C <- 0.004
    p$kappa_L   <- 0.3
    # visit_times unchanged (0,1,3,6,12); the discrete comparator coarsens.
  } else {
    stop("Unknown ITT/PP scenario: ", name)
  }
  p
}

simulate_ittpp_dgm <- function(n, params = ittpp_default_params(),
                                intervention = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p  <- params
  vt <- p$visit_times
  tau <- p$tau
  dt <- p$dt
  K  <- length(vt) - 1L
  t_grid <- seq(0, tau, by = dt)
  visit_idx <- match(round(vt / dt) * dt, round(t_grid / dt) * dt)

  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  A0 <- if (is.null(intervention)) rbinom(n, 1, 0.5) else rep(as.integer(intervention), n)

  Z <- W1
  L_current <- Z + rnorm(n, 0, p$sigma_L)
  L_mat <- matrix(NA_real_, n, length(vt))
  L_mat[, 1] <- L_current

  U_T <- -log(runif(n))
  U_R <- -log(runif(n))
  U_C <- -log(runif(n))
  H_T <- rep(0, n)
  H_R <- rep(0, n)
  H_C <- rep(0, n)
  T_time <- rep(Inf, n)
  R_time <- rep(Inf, n)
  C_time <- rep(Inf, n)
  # Current treatment, flips once at the switch time. Switching is a treatment
  # CHANGE (crossover to the other arm), not a terminal event: the subject is
  # still followed for failure afterward, but with the OTHER arm's hazard. This
  # is what makes ITT (follow through switch) differ from PP (no switch).
  A_curr  <- A0
  on_prot <- rep(TRUE, n)                  # not yet switched
  alive   <- rep(TRUE, n)                  # still under observation: T,C not yet fired

  next_visit <- 2L
  for (j in 2:length(t_grid)) {
    t_now <- t_grid[j]
    idx <- which(alive)
    if (!length(idx)) break

    # OU latent process (drift uses current treatment)
    drift <- (p$alpha_0 + p$alpha_A * A_curr[idx] +
                p$alpha_1 * W1[idx] + p$alpha_2 * W2[idx] -
                p$rho * Z[idx]) * dt
    Z[idx] <- Z[idx] + drift + p$sigma * sqrt(dt) * rnorm(length(idx))

    # Failure hazard uses CURRENT treatment (changes at switch)
    lambda_T <- p$lambda0_T * exp(
      p$beta_A * A_curr[idx] + p$beta_L * L_current[idx] +
        p$beta_Z * Z[idx] + p$beta_W1 * W1[idx] + p$beta_W2 * W2[idx]
    )
    H_T[idx] <- H_T[idx] + lambda_T * dt

    # Switch hazard fires only while on-protocol; conditions on assigned arm A0
    op <- idx[on_prot[idx]]
    if (length(op) && p$lambda0_R > 0) {
      lambda_R <- p$lambda0_R * exp(
        p$eta_A * A0[op] + p$eta_L * L_current[op] +
          p$eta_Z * Z[op] + p$eta_W1 * W1[op] + p$eta_W2 * W2[op]
      )
      H_R[op] <- H_R[op] + lambda_R * dt
    }

    # Censoring hazard (terminal for observation)
    if (p$lambda0_C > 0) {
      lambda_C <- p$lambda0_C * exp(
        p$kappa_A * A0[idx] + p$kappa_L * L_current[idx] +
          p$kappa_Z * Z[idx] + p$kappa_W1 * W1[idx] + p$kappa_W2 * W2[idx]
      )
      H_C[idx] <- H_C[idx] + lambda_C * dt
    }

    # --- Switch event (non-terminal): flip current treatment ---
    sw <- op[H_R[op] >= U_R[op]]
    if (length(sw)) {
      R_time[sw] <- t_now
      A_curr[sw] <- 1L - A0[sw]
      on_prot[sw] <- FALSE
    }

    # --- Terminal events: failure T and censoring C (vectorized) ---
    cross_T <- H_T[idx] >= U_T[idx]
    cross_C <- (p$lambda0_C > 0) & (H_C[idx] >= U_C[idx])
    if (any(cross_T) || any(cross_C)) {
      both <- cross_T & cross_C
      chose_T <- logical(length(idx))
      if (any(both)) {
        # resolve simultaneous failure/censoring by current cause-specific rates
        lt <- p$lambda0_T * exp(
          p$beta_A * A_curr[idx][both] + p$beta_L * L_current[idx][both] +
            p$beta_Z * Z[idx][both] + p$beta_W1 * W1[idx][both] + p$beta_W2 * W2[idx][both])
        lc <- p$lambda0_C * exp(
          p$kappa_A * A0[idx][both] + p$kappa_L * L_current[idx][both] +
            p$kappa_Z * Z[idx][both] + p$kappa_W1 * W1[idx][both] + p$kappa_W2 * W2[idx][both])
        chose_T[both] <- runif(sum(both)) < lt / (lt + lc)
      }
      assign_T <- (cross_T & !cross_C) | (both & chose_T)
      assign_C <- (cross_C & !cross_T) | (both & !chose_T)
      if (any(assign_T)) T_time[idx[assign_T]] <- t_now
      if (any(assign_C)) C_time[idx[assign_C]] <- t_now
      alive[idx[cross_T | cross_C]] <- FALSE
    }

    # Visit boundary: update covariates for those still under observation
    if (next_visit <= length(vt) && j == visit_idx[next_visit]) {
      idx2 <- which(alive)
      if (length(idx2)) {
        L_current[idx2] <- Z[idx2] + rnorm(length(idx2), 0, p$sigma_L)
        L_mat[idx2, next_visit] <- L_current[idx2]
      }
      next_visit <- next_visit + 1L
    }
  }

  make_ittpp_long(
    W1 = W1, W2 = W2, A0 = A0, L_mat = L_mat,
    T_time = T_time, R_time = R_time, C_time = C_time,
    visit_times = vt, tau = tau, params = p
  )
}

make_ittpp_long <- function(W1, W2, A0, L_mat, T_time, R_time, C_time,
                             visit_times, tau, params = NULL) {
  n <- length(W1)
  K <- length(visit_times) - 1L
  if (is.null(params)) params <- ittpp_default_params()

  # Switching is NOT terminal: observation ends at min(T, C, tau). Post-switch
  # failures are observed (and used by the ITT analysis). The PP analysis later
  # censors at the switch; that truncation happens in ittpp_analysis_long().
  obs_end <- pmin(T_time, C_time, tau)

  subj <- data.table(
    id = seq_len(n), W1 = W1, W2 = W2, A0 = A0,
    T_time = T_time, R_time = R_time, C_time = C_time,
    first_event_time = obs_end,
    failure_observed = as.integer(T_time <= pmin(C_time, tau)),
    switch_observed  = as.integer(R_time <= pmin(T_time, C_time, tau)),
    censor_observed  = as.integer(C_time < T_time & C_time <= tau)
  )

  # Build person-interval rows over the ITT (full follow-up) risk set. The
  # estimand-specific event indicators, time-at-risk, and risk set are computed
  # in ittpp_analysis_long() from the carried (T, R, C) times.
  long <- vector("list", n * K)
  pos <- 1L
  for (i in seq_len(n)) {
    for (k0 in 0:(K - 1L)) {
      t0 <- visit_times[k0 + 1L]
      t1 <- visit_times[k0 + 2L]
      if (obs_end[i] <= t0 + 1e-12) next
      Lk <- L_mat[i, k0 + 1L]
      if (!is.finite(Lk)) next

      long[[pos]] <- data.table(
        id = i, k = k0, t_start = t0, t_end = t1, ell = t1 - t0,
        W1 = W1[i], W2 = W2[i], A0 = A0[i],
        L_k = Lk,
        L_next = if (k0 < K - 1L && obs_end[i] > t1 + 1e-12) L_mat[i, k0 + 2L] else NA_real_,
        T_time = T_time[i], R_time = R_time[i], C_time = C_time[i]
      )
      pos <- pos + 1L
    }
  }

  pi <- rbindlist(long[seq_len(pos - 1L)])
  list(
    subject = subj, person_interval = pi, visit_L = L_mat,
    visit_times = visit_times, tau = tau, params = params
  )
}

simulate_ittpp_truth_subject <- function(n, params, intervention, estimand = "PP",
                                          seed = NULL) {
  # Simulate under intervention d:
  #   ITT: fix A0=intervention, set C=Inf, leave R natural
  #   PP:  fix A0=intervention, set C=Inf, set R=Inf
  if (!is.null(seed)) set.seed(seed)
  p  <- params
  vt <- p$visit_times
  tau <- p$tau
  dt <- p$dt
  t_grid <- seq(0, tau, by = dt)
  visit_idx <- match(round(vt / dt) * dt, round(t_grid / dt) * dt)

  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  A0 <- rep(as.integer(intervention), n)
  Z <- W1
  L_current <- Z + rnorm(n, 0, p$sigma_L)

  U_T <- -log(runif(n))
  U_R <- -log(runif(n))
  H_T <- rep(0, n)
  H_R <- rep(0, n)
  T_time <- rep(Inf, n)
  R_time <- rep(Inf, n)
  alive <- rep(TRUE, n)            # at risk for failure (T not yet fired)

  # PP: never switch (A_curr fixed at A0). ITT: switch fires naturally and flips
  # the current treatment (crossover), but the subject is still followed for T.
  # Censoring is always removed (truth is failure-free survival under d).
  suppress_R <- identical(estimand, "PP")
  A_curr  <- A0
  on_prot <- rep(TRUE, n)

  next_visit <- 2L
  for (j in 2:length(t_grid)) {
    t_now <- t_grid[j]
    idx <- which(alive)
    if (!length(idx)) break

    drift <- (p$alpha_0 + p$alpha_A * A_curr[idx] +
                p$alpha_1 * W1[idx] + p$alpha_2 * W2[idx] -
                p$rho * Z[idx]) * dt
    Z[idx] <- Z[idx] + drift + p$sigma * sqrt(dt) * rnorm(length(idx))

    # Failure hazard uses CURRENT treatment
    lambda_T <- p$lambda0_T * exp(
      p$beta_A * A_curr[idx] + p$beta_L * L_current[idx] +
        p$beta_Z * Z[idx] + p$beta_W1 * W1[idx] + p$beta_W2 * W2[idx]
    )
    H_T[idx] <- H_T[idx] + lambda_T * dt

    # Switch (ITT only), among on-protocol subjects; flips current treatment
    if (!suppress_R && p$lambda0_R > 0) {
      op <- idx[on_prot[idx]]
      if (length(op)) {
        lambda_R <- p$lambda0_R * exp(
          p$eta_A * A0[op] + p$eta_L * L_current[op] +
            p$eta_Z * Z[op] + p$eta_W1 * W1[op] + p$eta_W2 * W2[op]
        )
        H_R[op] <- H_R[op] + lambda_R * dt
        sw <- op[H_R[op] >= U_R[op]]
        if (length(sw)) {
          R_time[sw] <- t_now
          A_curr[sw] <- 1L - A0[sw]
          on_prot[sw] <- FALSE
        }
      }
    }

    # Failure is the only terminal event (censoring removed for truth)
    cross_T <- H_T[idx] >= U_T[idx]
    if (any(cross_T)) {
      T_time[idx[cross_T]] <- t_now
      alive[idx[cross_T]] <- FALSE
    }

    if (next_visit <= length(vt) && j == visit_idx[next_visit]) {
      idx2 <- which(alive)
      if (length(idx2)) L_current[idx2] <- Z[idx2] + rnorm(length(idx2), 0, p$sigma_L)
      next_visit <- next_visit + 1L
    }
  }

  data.table(T_time = T_time, R_time = R_time)
}

simulate_ittpp_toy_constant <- function(n, visit_times, lambda_T, lambda_R = 0,
                                         lambda_C = 0, treatment = FALSE,
                                         lambda_T_by_A = NULL, seed = NULL) {
  # Toy DGM with constant hazards, no covariates
  if (!is.null(seed)) set.seed(seed)
  K <- length(visit_times) - 1L
  tau <- max(visit_times)
  lt_vec <- rep(lambda_T, length.out = K)
  lr_vec <- rep(lambda_R, length.out = K)
  lc_vec <- rep(lambda_C, length.out = K)

  W1 <- rep(0, n)
  W2 <- rep(0, n)
  A0 <- if (treatment) rbinom(n, 1, 0.5) else rep(0L, n)
  L_mat <- matrix(0, n, length(visit_times))
  T_time <- rep(Inf, n)
  R_time <- rep(Inf, n)
  C_time <- rep(Inf, n)

  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      if (is.finite(T_time[i]) || is.finite(R_time[i]) || is.finite(C_time[i])) break
      lt <- lt_vec[k]
      lr <- lr_vec[k]
      lc <- lc_vec[k]
      if (!is.null(lambda_T_by_A)) lt <- lambda_T_by_A[as.character(A0[i])]
      total_rate <- lt + lr + lc
      if (total_rate <= 0) next
      wait <- rexp(1, rate = total_rate)
      ell <- visit_times[k + 1L] - visit_times[k]
      if (wait <= ell) {
        u <- runif(1)
        if (u < lt / total_rate) {
          T_time[i] <- visit_times[k] + wait
        } else if (u < (lt + lr) / total_rate) {
          R_time[i] <- visit_times[k] + wait
        } else {
          C_time[i] <- visit_times[k] + wait
        }
      }
    }
  }

  make_ittpp_long(W1, W2, A0, L_mat, T_time, R_time, C_time, visit_times, tau)
}
