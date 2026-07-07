# Estimate a composite risk difference with VA-CT-TMLE

Estimate a composite risk difference with VA-CT-TMLE

## Usage

``` r
va_ct_tmle(
  obs_data,
  d0 = 0L,
  d1 = 1L,
  tau = obs_data$tau,
  g0 = 0.5,
  M = 100,
  B = 10,
  max_iter = 10,
  seed = 1L,
  endpoint = "composite",
  update_method = c("adaptive", "standard"),
  step_size = 1,
  min_step_size = 1e-04,
  max_line_iter = 10
)
```

## Arguments

- obs_data:

  Visit-aligned analysis object.

- d0, d1:

  Static binary regimes to compare. The reported risk difference is risk
  under `d1` minus risk under `d0`.

- tau:

  Analysis horizon.

- g0:

  Baseline randomization probability for treatment 1.

- M:

  Monte Carlo draws for visit-level integration.

- B:

  Number of sub-bins for hazard targeting.

- max_iter:

  Maximum targeting iterations.

- seed:

  Random seed.

- endpoint:

  Endpoint type. Only `"composite"` is supported.

- update_method:

  Targeting update method.

- step_size, min_step_size, max_line_iter:

  Adaptive update controls.

## Value

A `vacttmle` object.
