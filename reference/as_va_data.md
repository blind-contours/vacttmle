# Build a visit-aligned analysis object

Converts subject-level and person-interval tables into the analysis
object used by the package estimators.

## Usage

``` r
as_va_data(subject_data, interval_data, visit_times, tau, visit_L = NULL,
  visit_A = NULL, visit_g = NULL, params = list())
```

## Arguments

- subject_data:

  Subject-level data.

- interval_data:

  Person-interval data.

- visit_times:

  Scheduled visit times.

- tau:

  Analysis horizon.

- visit_L, visit_A, visit_g:

  Optional visit matrices.

- params:

  Optional metadata list.

## Value

A `va_data` object.
