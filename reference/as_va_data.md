# Build a visit-aligned analysis object

Converts subject-level and person-interval tables into the analysis
object used by the package estimators. The person-interval table must
include start-of-visit and next-visit covariates (\`L_k\`, \`L_next\`)
plus exact elapsed event time within the interval (\`event_time\`) so
the continuous-time hazard targeting step can be computed.

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
