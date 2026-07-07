# Build a visit-aligned analysis object

`as_va_data()` converts subject-level and person-interval tables into
the analysis object used by the composite-endpoint VA-CT-TMLE
estimators.

## Usage

``` r
as_va_data(
  subject_data,
  interval_data,
  visit_times,
  tau,
  visit_L = NULL,
  visit_A = NULL,
  visit_g = NULL,
  params = list()
)
```

## Arguments

- subject_data:

  Subject-level data. Required columns are `id`, `W1`, `W2`, `A0`,
  `T_time`, `D_time`, `composite_time`, and `composite_event`.

- interval_data:

  Person-interval data. Required columns are `id`, `k`, `t_start`,
  `t_end`, `ell`, `W1`, `W2`, `A0`, `A_k`, `L_k`, `time_at_risk`,
  `event_T`, `event_D`, `event_comp`, and `Y_end`.

- visit_times:

  Numeric vector of scheduled visit times.

- tau:

  Analysis horizon.

- visit_L:

  Optional subject-by-visit matrix of visit covariates.

- visit_A:

  Optional subject-by-interval matrix of visit treatment values.

- visit_g:

  Optional subject-by-interval matrix of known treatment probabilities
  for the observed treatment path.

- params:

  Optional list of DGM or metadata parameters.

## Value

A list with class `va_data`.
