test_that("package toy checks pass at small n", {
  checks <- run_vacttmle_toy_checks(n = 3000, seed = 7001)
  expect_error({
    if (!all(checks$passed)) stop("one or more package toy checks failed")
  }, NA)
})
