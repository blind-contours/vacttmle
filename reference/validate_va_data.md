# Validate visit-aligned data

Checks required components, endpoint coding, and basic person-interval
consistency for the composite endpoint.

## Usage

``` r
validate_va_data(obs_data, endpoint = "composite")
```

## Arguments

- obs_data:

  Object produced by
  [`as_va_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_data.md)
  or
  [`simulate_va_trial()`](https://blind-contours.github.io/vacttmle/reference/simulate_va_trial.md).

- endpoint:

  Endpoint type. Only `"composite"` is supported.

## Value

The validated object, invisibly.
