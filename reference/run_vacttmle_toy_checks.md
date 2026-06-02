# Run package toy checks

Runs lightweight analytic validation checks for the composite survival
target.

## Usage

``` r
run_vacttmle_toy_checks(n = 5000, seed = 7001)
```

## Arguments

- n:

  Analytic toy-simulation sample size.

- seed:

  Random seed.

## Value

A data.table with estimates, analytic truths, errors, and pass flags.
