# Estimate ITT or per-protocol failure-free survival with VA-CT-TMLE

Estimate ITT or per-protocol failure-free survival with VA-CT-TMLE

## Usage

``` r
va_ct_switch(
  obs_data,
  estimand = c("ITT", "PP"),
  tau = obs_data$tau,
  g0 = 0.5,
  weight_trunc = NULL,
  M = 100,
  B = 10,
  max_iter = 10,
  seed = 1L
)
```

## Arguments

- obs_data:

  A `va_switch_data` object from
  [`as_va_switch_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_switch_data.md)
  or
  [`simulate_va_switch_trial()`](https://blind-contours.github.io/vacttmle/reference/simulate_va_switch_trial.md).

- estimand:

  `"ITT"` (treatment-policy: natural switching, censoring removed) or
  `"PP"` (hypothetical no-switching: switching and censoring removed).

- tau:

  Analysis horizon.

- g0:

  Known randomization probability of arm 1.

- weight_trunc:

  Optional cumulative-weight control under positivity stress: a numeric
  absolute cap, or a string like `"p95"` for the 95th percentile of the
  positive weights (recommended when switching is heavy). `NULL`
  (default) applies no truncation.

- M:

  Monte Carlo draws for visit-level integration.

- B:

  Sub-bins for within-interval hazard targeting.

- max_iter:

  Maximum targeting iterations.

- seed:

  Random seed.

## Value

A `vacttmle` object; the reported risk difference is risk under arm 1
minus risk under arm 0 (failure by `tau`).
