# Simulate a visit-aligned trial

Generates example scheduled-visit trial data for package examples and
tests.

## Usage

``` r
simulate_va_trial(n = 1000, scenario = "S1", intervention = NULL, seed = NULL)
```

## Arguments

- n:

  Number of subjects.

- scenario:

  Scenario name.

- intervention:

  Optional static intervention, 0 or 1. Leave `NULL` for the observed
  treatment mechanism.

- seed:

  Optional random seed.

## Value

A `va_data` object.
