test_that("public API returns bounded composite estimates", {
  dat <- simulate_va_trial(n = 200, scenario = "S1", seed = 101)
  fit <- va_ct_gcomp(dat, M = 5, seed = 102)
  expect_s3_class(fit, "vacttmle")
  if (!is.finite(fit$result$risk_difference)) fail("risk difference is not finite")
  if (!(fit$result$risk_d0 >= 0 && fit$result$risk_d0 <= 1)) fail("risk_d0 out of bounds")
  if (!(fit$result$risk_d1 >= 0 && fit$result$risk_d1 <= 1)) fail("risk_d1 out of bounds")
})

test_that("validation enforces composite event coding", {
  dat <- simulate_va_trial(n = 50, scenario = "S1", seed = 103)
  dat$person_interval$event_comp[1] <- 1L - dat$person_interval$event_comp[1]
  expect_error(validate_va_data(dat), "event_comp")
})
