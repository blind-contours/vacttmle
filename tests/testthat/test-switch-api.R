test_that("va_ct_switch estimates both estimands on simulated data", {
  dat <- simulate_va_switch_trial(400, "A2", seed = 11)
  pp <- va_ct_switch(dat, "PP", M = 30, B = 5, max_iter = 3, seed = 11)
  itt <- va_ct_switch(dat, "ITT", M = 30, B = 5, max_iter = 3, seed = 11)
  expect_s3_class(pp, "vacttmle")
  expect_true(is.finite(pp$risk_difference))
  expect_true(is.finite(pp$se) && pp$se > 0)
  expect_lt(pp$ci["lower"], pp$ci["upper"])
  # crossover dilutes ITT relative to PP (protective effect)
  expect_lt(pp$risk_difference, 0)
  expect_lte(abs(itt$risk_difference), abs(pp$risk_difference) + 0.05)
})

test_that("weight truncation runs and caps weights", {
  dat <- simulate_va_switch_trial(400, "A7", seed = 12)
  fit <- va_ct_switch(dat, "PP", weight_trunc = "p95", M = 30, B = 5,
                      max_iter = 3, seed = 12)
  expect_true(is.finite(fit$risk_difference))
  mw <- suppressWarnings(max(unlist(fit$diagnostics[, c("max_weight_d0",
        "max_weight_d1")]), na.rm = TRUE))
  expect_lt(mw, 100)
})

test_that("as_va_switch_data accepts tidy input and round-trips", {
  sim <- simulate_va_switch_trial(200, "A2", seed = 13)
  s <- sim$subject
  df <- data.frame(id = s$id, A0 = s$A0, W1 = s$W1, W2 = s$W2,
                   T_time = s$T_time, R_time = s$R_time, C_time = s$C_time)
  dat <- as_va_switch_data(df, sim$visit_L, visit_times = c(0, 1, 3, 6, 12),
                           tau = 12)
  expect_s3_class(dat, "va_switch_data")
  expect_gt(nrow(dat$person_interval), 0)
})

test_that("as_va_switch_data rejects malformed input with clear messages", {
  sim <- simulate_va_switch_trial(50, "A2", seed = 14)
  s <- sim$subject
  ok <- data.frame(id = s$id, A0 = s$A0, W1 = s$W1, W2 = s$W2,
                   T_time = s$T_time, R_time = s$R_time, C_time = s$C_time)
  vt <- c(0, 1, 3, 6, 12)

  bad <- ok; bad$A0 <- bad$A0 + 1L            # arms coded 1/2
  expect_error(as_va_switch_data(bad, sim$visit_L, vt, 12), "coded 0/1")

  bad <- ok; bad$R_time[1] <- NA              # NA instead of Inf
  expect_error(as_va_switch_data(bad, sim$visit_L, vt, 12), "Inf")

  bad <- ok; bad$T_time[2] <- -1              # negative time
  expect_error(as_va_switch_data(bad, sim$visit_L, vt, 12), "negative")

  bad <- ok; bad$id[2] <- bad$id[1]           # duplicate ids
  expect_error(as_va_switch_data(bad, sim$visit_L, vt, 12), "duplicates")

  expect_error(as_va_switch_data(ok, sim$visit_L, c(1, 3, 6, 12), 12),
               "start at 0")                  # grid not anchored at 0
  expect_error(as_va_switch_data(ok, sim$visit_L, vt, 24),
               "tau")                         # horizon beyond grid
  expect_error(as_va_switch_data(ok, sim$visit_L[, 1:3], vt, 12),
               "one column per visit")        # wrong visit_L shape

  badL <- sim$visit_L; badL[1, 1] <- NA       # missing baseline L
  expect_error(as_va_switch_data(ok, badL, vt, 12), "baseline")
})

test_that("edge cases: zero switchers and zero censoring still run", {
  dat <- simulate_va_switch_trial(300, "A0", seed = 15)   # no R, no C
  fit <- va_ct_switch(dat, "PP", M = 30, B = 5, max_iter = 3, seed = 15)
  expect_true(is.finite(fit$risk_difference))
})

test_that("discrete-LTMLE comparator runs and orderings differ under stress", {
  dat <- simulate_va_switch_trial(400, "A7", seed = 16)
  a <- est_discrete_ltmle(dat, "PP", delta = 3, ordering = "switch_first",
                          M = 20, seed = 16)
  b <- est_discrete_ltmle(dat, "PP", delta = 3, ordering = "event_first",
                          M = 20, seed = 16)
  expect_true(is.finite(a$risk_difference) && is.finite(b$risk_difference))
  expect_false(isTRUE(all.equal(a$risk_difference, b$risk_difference)))
})
