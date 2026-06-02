# vacttmle

`vacttmle` is a standalone R package for visit-aligned continuous-time TMLE in
scheduled-visit trials with a composite time-to-event endpoint.

The current package targets

```text
P(T^d > tau, D^d > tau)
```

and the composite risk difference

```text
{1 - P(T^d1 > tau, D^d1 > tau)} -
{1 - P(T^d0 > tau, D^d0 > tau)}.
```

`D` is part of the composite endpoint. It is not treated as censoring in this
package version.

## Minimal example

```r
library(vacttmle)

dat <- simulate_va_trial(n = 1000, scenario = "S1", seed = 1)

fit <- va_ct_tmle(dat, M = 100, B = 10, max_iter = 10, seed = 2)
fit

fit$result
fit$diagnostics
```

## Public functions

- `va_ct_tmle()` estimates the composite risk difference with VA-CT-TMLE.
- `va_visit_tmle()` runs the visit-node-only ablation.
- `va_ct_gcomp()` runs visit-aligned continuous-time g-computation without targeting.
- `as_va_data()` adapts subject and person-interval tables to the package format.
- `validate_va_data()` checks endpoint coding and required fields.
- `simulate_va_trial()` generates example scheduled-visit trial data.
- `run_vacttmle_toy_checks()` runs lightweight analytic checks.

## Current scope

This is an early standalone package cut. The primary estimand is the composite
endpoint only. Hypothetical no-switch/no-intercurrent-event estimands require a
different estimator and are intentionally not exposed here.
