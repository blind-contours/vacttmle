###############################################################################
# Faithful discrete-time LTMLE comparator with an explicit within-interval
# ORDERING CONVENTION.
#
# The distinction vs. VA-CT-TMLE: a discrete LTMLE analyst must (1) SNAP every
# event/switch/censor time to a bin boundary of width Delta, and (2) when a
# switch and the failure (or a dropout and the failure) fall in the SAME bin,
# impose an arbitrary ordering to decide who "owns" the event. We implement the
# snapping + ordering as a data transform, then run the SAME validated discrete
# ICE-LTMLE recursion (ittpp_estimate_contrast on the snapped grid, which has no
# within-interval events, so it IS standard discrete LTMLE). The ordering
# convention is the only thing that changes between variants on identical data.
#
# Conventions (PP estimand, where switching is handled as artificial censoring):
#   "event_first" : L -> Y -> R -> C  in each bin. A switch and a death in the
#                   same bin -> the DEATH is counted (subject fails under the
#                   assigned regimen). ltmle default reading: outcome adjudicated
#                   before the treatment/censoring node it shares a bin with is
#                   NOT how ltmle orders, but it is the "event-before-switch"
#                   tie-break analysts use when the death is documented at the
#                   same scan that documents progression-then-switch.
#   "switch_first": L -> R -> C -> Y  in each bin (matches the ltmle/stremr
#                   ordering: censoring resolved before the outcome). A switch
#                   and a death in the same bin -> subject is CENSORED at switch,
#                   the death is NOT counted (standard "censor-at-switch").
#   "drop_first"  : like switch_first but the co-occurrence toggled is dropout
#                   vs. death (immortal-time contrast on the censoring side).
#
# Both event_first and switch_first are defensible, published conventions; they
# give different PP survival on the same data when co-occurrence is common.
###############################################################################

# Snap a scalar time to the Delta-grid bin BOUNDARY (event dated to the detecting
# assessment = interval end), capped at tau.
.snap_up <- function(t, delta, tau) {
  if (!is.finite(t)) return(Inf)
  pmin(ceiling(t / delta - 1e-9) * delta, tau)
}

#' Snap a continuous-time trial to a discrete grid under an ordering convention
#'
#' Converts a `va_switch_data` object to the discretized data a discrete-time
#' LTMLE analyst would construct: every failure/switch/censoring time is snapped
#' to its bin boundary, and same-bin co-occurrences are resolved by the chosen
#' within-interval ordering convention.
#'
#' @param obs_data A `va_switch_data` object (continuous-time; see
#'   [as_va_switch_data()]).
#' @param delta Bin width, on the time scale of the data.
#' @param ordering Within-bin ordering convention: `"switch_first"` (switch
#'   resolved before the outcome; a same-bin death is censored at switch),
#'   `"event_first"` (the death is adjudicated first and counted), or
#'   `"drop_first"` (dropout resolved before the outcome).
#' @param tau Analysis horizon.
#'
#' @return A discretized `va_switch_data`-style object on the `delta` grid.
#' @seealso [est_discrete_ltmle()]
snap_to_grid <- function(obs_data, delta, ordering = c("switch_first",
                                                       "event_first", "drop_first"),
                         tau = obs_data$tau) {
  ordering <- match.arg(ordering)
  subj <- copy(obs_data$subject)
  n <- nrow(subj)
  edges <- seq(0, tau, by = delta)
  if (abs(edges[length(edges)] - tau) > 1e-9) edges <- c(edges, tau)

  Ts <- vapply(subj$T_time, .snap_up, 0.0, delta = delta, tau = tau)
  Rs <- vapply(subj$R_time, .snap_up, 0.0, delta = delta, tau = tau)
  Cs <- vapply(subj$C_time, .snap_up, 0.0, delta = delta, tau = tau)

  # Resolve same-bin co-occurrence per the ordering convention by nudging the
  # snapped times to enforce the intended precedence within the shared bin.
  eps <- delta * 1e-3
  same_RT <- is.finite(Rs) & is.finite(Ts) & abs(Rs - Ts) < eps
  same_CT <- is.finite(Cs) & is.finite(Ts) & abs(Cs - Ts) < eps

  if (ordering == "event_first") {
    # death adjudicated before switch/dropout in the shared bin: keep T, push
    # R and C just AFTER T so the failure is observed (not censored).
    Rs[same_RT] <- Ts[same_RT] + eps
    Cs[same_CT] <- Ts[same_CT] + eps
  } else if (ordering == "switch_first") {
    # switch resolved before the outcome: push R just BEFORE T so the subject is
    # censored-at-switch and the death is not counted. Dropout keeps default.
    Rs[same_RT] <- pmax(Ts[same_RT] - eps, 0)
  } else if (ordering == "drop_first") {
    # dropout resolved before the outcome (immortal-time contrast on censoring):
    Cs[same_CT] <- pmax(Ts[same_CT] - eps, 0)
    Rs[same_RT] <- pmax(Ts[same_RT] - eps, 0)
  }

  subj[, T_time := Ts]
  subj[, R_time := Rs]
  subj[, C_time := Cs]
  subj[, first_event_time := pmin(T_time, C_time, tau)]
  subj[, failure_observed := as.integer(T_time <= pmin(C_time, tau))]
  subj[, switch_observed  := as.integer(R_time <= pmin(T_time, C_time, tau))]
  subj[, censor_observed  := as.integer(C_time < T_time & C_time <= tau)]

  # Rebuild person-interval on the Delta-grid with LOCF covariates from the
  # original visit measurements (covariates seen only at scan visits).
  vt_orig <- obs_data$visit_times
  L_orig  <- obs_data$visit_L                # n x length(vt_orig)
  K <- length(edges) - 1L
  # LOCF value at a query time = last visit covariate at or before that time.
  locf_L <- function(i, tq) {
    idx <- which(vt_orig <= tq + 1e-9)
    if (!length(idx)) return(L_orig[i, 1])
    v <- L_orig[i, max(idx)]
    if (!is.finite(v)) { av <- L_orig[i, idx]; av <- av[is.finite(av)]
      if (length(av)) av[length(av)] else L_orig[i, 1] } else v
  }

  L_grid <- matrix(NA_real_, n, length(edges))
  rows <- vector("list", n * K); pos <- 1L
  for (i in seq_len(n)) {
    oe <- subj$first_event_time[i]
    for (k0 in seq_len(K) - 1L) {
      t0 <- edges[k0 + 1L]; t1 <- edges[k0 + 2L]
      if (oe <= t0 + 1e-12) break
      Lk <- locf_L(i, t0)
      L_grid[i, k0 + 1L] <- Lk
      rows[[pos]] <- data.table(
        id = i, k = k0, t_start = t0, t_end = t1, ell = t1 - t0,
        W1 = subj$W1[i], W2 = subj$W2[i], A0 = subj$A0[i],
        L_k = Lk, L_next = NA_real_,
        T_time = subj$T_time[i], R_time = subj$R_time[i], C_time = subj$C_time[i])
      pos <- pos + 1L
    }
  }
  pi <- rbindlist(rows[seq_len(pos - 1L)])
  # fill L_next from the next in-grid row
  setorder(pi, id, k)
  pi[, L_next := shift(L_k, type = "lead"), by = id]

  list(subject = subj, person_interval = pi, visit_L = L_grid,
       visit_times = edges, tau = tau, params = obs_data$params)
}

#' Discrete-time LTMLE contrast under a fixed grid and ordering convention
#'
#' Runs a genuine discrete-time ICE-LTMLE on the output of [snap_to_grid()]:
#' event/switch/censoring times live on the bin boundaries, so this is the
#' estimator a discrete-time analyst would compute, including its dependence on
#' the arbitrary within-bin ordering convention. Used to reproduce the
#' ordering-sensitivity experiment of the accompanying paper.
#'
#' @param obs_data A `va_switch_data` object (continuous-time).
#' @param estimand `"ITT"` or `"PP"`.
#' @param tau Analysis horizon.
#' @param g0 Known randomization probability of arm 1.
#' @param delta Bin width passed to [snap_to_grid()].
#' @param ordering Within-bin ordering convention passed to [snap_to_grid()].
#' @param M Monte Carlo draws for visit-level integration.
#' @param B Sub-bins (kept at 1 for a genuine discrete-time fit).
#' @param seed Random seed.
#' @param true_rd Optional true risk difference, carried into the result.
#' @param scenario,replicate_id Optional bookkeeping labels.
#'
#' @return A one-row `data.table` with the estimate, SE, CI, and metadata.
#' @seealso [snap_to_grid()], [va_ct_switch()]
est_discrete_ltmle <- function(obs_data, estimand = "PP", tau = obs_data$tau,
                               g0 = 0.5, delta = 1, ordering = "switch_first",
                               M = 50, B = 1, seed = NULL, true_rd = NA_real_,
                               scenario = NA_character_, replicate_id = NA_integer_) {
  t0 <- proc.time()[["elapsed"]]
  snapped <- snap_to_grid(obs_data, delta = delta, ordering = ordering, tau = tau)
  # B = 1: no within-interval sub-binning (events are at boundaries by
  # construction), i.e. a genuine discrete-time ICE-LTMLE.
  res <- tryCatch(
    ittpp_estimate_contrast(snapped, estimator = "ct_tmle", estimand = estimand,
                            tau = tau, g0 = g0, M = M, B = 1L, max_iter = 5,
                            seed = seed, true_rd = true_rd, scenario = scenario,
                            replicate_id = replicate_id, update_method = "adaptive"),
    error = function(e) NULL)
  if (is.null(res)) {
    return(data.table(scenario = scenario, replicate_id = replicate_id,
      n = nrow(obs_data$subject), tau = tau,
      estimator = paste0("dltmle_", ordering), estimand = estimand,
      delta = delta, ordering = ordering,
      risk_difference = NA_real_, true_risk_difference = true_rd, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0))
  }
  res[, estimator := paste0("dltmle_", ordering)]
  res[, delta := delta]
  res[, ordering := ordering]
  res[, runtime_seconds := proc.time()[["elapsed"]] - t0]
  res[]
}
