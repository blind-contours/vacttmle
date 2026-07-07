# data.table non-standard-evaluation column names used inside [.data.table
# calls. Declared here so R CMD check does not flag them as undefined globals.
utils::globalVariables(c(
  "C_time", "Cplus", "G_cum", "G_int", "R_time", "T_time", "Xend",
  "censor_observed", "estimator", "event_C", "event_R", "event_T",
  "event_time", "failure_observed", "first_event_time", "id", "k", "k_f",
  "L_k", "L_next", "log_pt", "ordering", "runtime_seconds", "switch_observed",
  "t_end", "t_start", "time_at_risk", "treat_w", "A0", "M", "Q", "S", "W1",
  "W2", "Y_end", "delta", "ell"
))
