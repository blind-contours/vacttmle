# Build a visit-aligned analysis object for the switching estimands

Converts a one-row-per-subject table of exact event times plus a
subject-by-visit matrix of the time-varying covariate into the object
consumed by
[`va_ct_switch()`](https://blind-contours.github.io/vacttmle/reference/va_ct_switch.md).

## Usage

``` r
as_va_switch_data(subject_data, visit_L, visit_times, tau = max(visit_times))
```

## Arguments

- subject_data:

  Data frame with one row per subject and columns `id`, `A0` (randomized
  arm, 0/1), `W1`, `W2` (baseline covariates), `T_time` (failure time),
  `R_time` (treatment-switch time; `Inf` if never), `C_time`
  (censoring/dropout time; `Inf` if never). Times are on the same scale
  as `visit_times` and `tau`; administrative censoring beyond `tau`
  should be coded as `C_time = Inf`.

- visit_L:

  Numeric matrix, one row per subject (row order matching
  `subject_data`) and one column per entry of `visit_times`, giving the
  visit covariate `L` measured at each scheduled visit. Values after a
  subject leaves observation may be `NA`.

- visit_times:

  Numeric vector of scheduled visit times, starting at 0.

- tau:

  Analysis horizon.

## Value

A `va_switch_data` object.

## Details

The current version uses two baseline covariates (`W1`, `W2`) and a
single scheduled-visit covariate `L` (measured at each `visit_times`),
matching the estimator in McCoy (2026). Map your trial's covariates onto
these columns.
