# Convergence diagnostics

This article explains the diagnostics returned by `vacttmle`.

## What convergence means

The targeting step tries to solve the empirical mean of each efficient
influence-function component:

- failure-intensity component, `D^T`
- intercurrent-event intensity component, `D^D`
- visit-level iterated conditional expectation component, `D^Q`

The default stopping threshold is

``` text
sigma(EIF) / (sqrt(n) * log(n)).
```

The fit is marked converged when all three maximum absolute component
means are below this threshold.

## Inspect component means

``` r
fit <- va_ct_tmle(obs, M = 100, B = 10, max_iter = 10, seed = 2026)
fit$diagnostics
```

Example:

| converged | n_outer_iterations | max_abs_eif_T | max_abs_eif_D | max_abs_eif_Q | max_abs_total_eif |
|-----------|-------------------:|--------------:|--------------:|--------------:|------------------:|
| TRUE      |                  1 |      0.000370 |      0.000424 |      0.001349 |          0.001094 |

Interpretation:

- `max_abs_eif_T`: largest failure-hazard EIF component mean across
  visits
- `max_abs_eif_D`: largest intercurrent-event-hazard EIF component mean
- `max_abs_eif_Q`: largest visit-level EIF component mean
- `max_abs_total_eif`: empirical mean of the total risk-difference EIF

Do not rely only on `max_abs_total_eif`. Component means can fail while
the total mean cancels numerically.

## Inspect weights

The diagnostic output includes cumulative treatment-history weight
summaries:

| Column                           | Meaning                                             |
|----------------------------------|-----------------------------------------------------|
| `max_weight_d0`, `max_weight_d1` | Maximum cumulative observed-regime weight           |
| `ess_d0`, `ess_d1`               | Effective sample size implied by cumulative weights |

For baseline randomized treatment with deterministic continuation,
weights should usually be `0` or `1 / g0`. With `g0 = 0.5`, the maximum
should be `2`.

For stochastic post-baseline treatment/deviation processes, high maximum
weights or low effective sample size indicate positivity stress. Treat
such analyses as diagnostic until a prespecified stabilization strategy
is chosen.

## If convergence fails

Recommended debugging ladder:

1.  Confirm the endpoint coding: `event_comp = event_T OR event_D`.
2.  Confirm that `D` is intended to be part of the composite endpoint.
3.  Run
    [`run_vacttmle_toy_checks()`](https://blind-contours.github.io/vacttmle/reference/run_vacttmle_toy_checks.md)
    to verify the installed estimator.
4.  Start with
    [`va_ct_gcomp()`](https://rdrr.io/pkg/vacttmle/man/va_ct_gcomp.html)
    to check initial nuisance behavior.
5.  Use `update_method = "adaptive"` for the TMLE fit.
6.  Increase `max_iter` only if component means are improving.
7.  Increase `M` to reduce Monte Carlo integration noise.
8.  Inspect cumulative weights and effective sample sizes.

Example:

``` r
fit <- va_ct_tmle(
  obs,
  M = 200,
  B = 10,
  max_iter = 20,
  update_method = "adaptive",
  seed = 2026
)

fit$diagnostics
```

## Common patterns

| Pattern              | Likely meaning                                 | Next step                                |
|----------------------|------------------------------------------------|------------------------------------------|
| High `max_abs_eif_T` | Failure-hazard targeting not solved            | Check event counts and hazard model fit  |
| High `max_abs_eif_D` | Intercurrent-event hazard targeting not solved | Check D incidence and endpoint coding    |
| High `max_abs_eif_Q` | Visit-level ICE targeting not solved           | Check `L_next` and survivor rows         |
| High max weights     | Positivity stress                              | Check treatment-history probabilities    |
| Low ESS              | Few subjects follow the regime                 | Treat sustained-regime result cautiously |

## Reporting template

``` r
library(data.table)

list(
  package_version = as.character(packageVersion("vacttmle")),
  tau = obs$tau,
  visit_times = obs$visit_times,
  endpoint = "composite",
  estimate = fit$result,
  diagnostics = fit$diagnostics,
  event_counts = obs$subject[, .N, by = .(A0, composite_event)]
)
```
