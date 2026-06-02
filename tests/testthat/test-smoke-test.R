test_that("installed smoke test runs", {
  smoke_env <- new.env(parent = globalenv())
  expect_error({
    source(system.file("examples", "trialist-smoke-test.R", package = "vacttmle"),
           local = smoke_env)
    if (!exists("smoke_summary", envir = smoke_env)) stop("missing smoke_summary")
    smoke_summary <- get("smoke_summary", envir = smoke_env)
    if (!identical(smoke_summary$status[1], "ok")) stop("smoke status was not ok")
  }, NA)
})
