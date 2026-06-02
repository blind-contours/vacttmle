library(vacttmle)
library(data.table)

cat("Running vacttmle trialist smoke test\n\n")

toy_checks <- run_vacttmle_toy_checks(n = 3000, seed = 7001)
print(toy_checks)
if (!all(toy_checks$passed)) {
  stop("Toy checks failed. Do not proceed to trial data until installation is fixed.")
}

dat <- simulate_va_trial(n = 500, scenario = "S1", seed = 1)

cat("\nComposite event counts by baseline arm\n")
print(dat$subject[, .N, by = .(A0, composite_event)][order(A0, composite_event)])

fit <- va_ct_tmle(dat, M = 20, B = 10, max_iter = 5, seed = 2)

cat("\nVA-CT-TMLE fit\n")
print(fit)

cat("\nDiagnostics\n")
print(fit$diagnostics)

smoke_summary <- data.table::data.table(
  analysis = "simulated_S1",
  status = if (isTRUE(fit$result$converged) && is.finite(fit$result$risk_difference)) "ok" else "check",
  risk_d0 = fit$result$risk_d0,
  risk_d1 = fit$result$risk_d1,
  risk_difference = fit$result$risk_difference,
  se = fit$result$se,
  ci_lower = fit$result$ci_lower,
  ci_upper = fit$result$ci_upper,
  converged = fit$result$converged,
  max_weight = max(fit$result$max_weight_d0, fit$result$max_weight_d1),
  min_ess = min(fit$result$ess_d0, fit$result$ess_d1)
)

cat("\nSmoke-test summary\n")
print(smoke_summary)

invisible(smoke_summary)
