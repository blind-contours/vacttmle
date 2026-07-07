# Simulate a visit-aligned trial

This helper is intended for examples, tests, and teaching. It returns
the same `va_data` object accepted by the estimators.

## Usage

``` r
simulate_va_trial(n = 1000, scenario = "S1", intervention = NULL, seed = NULL)
```

## Arguments

- n:

  Number of subjects.

- scenario:

  Scenario name. Currently `"S1"`, `"S3"`, `"S4c"`, and `"S6"` are
  available.

- intervention:

  Optional static intervention, `0` or `1`. Leave `NULL` to sample
  treatment from the observed treatment mechanism.

- seed:

  Optional random seed.

## Value

A `va_data` object.
