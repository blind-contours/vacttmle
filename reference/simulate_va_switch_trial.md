# Simulate a scheduled-visit trial with switching and censoring

Convenience generator for examples, tests, and teaching; returns the
same `va_switch_data` object accepted by
[`va_ct_switch()`](https://blind-contours.github.io/vacttmle/reference/va_ct_switch.md).

## Usage

``` r
simulate_va_switch_trial(n = 1000, scenario = "A2", seed = NULL)
```

## Arguments

- n:

  Number of subjects.

- scenario:

  Scenario name: `"A0"`–`"A7"` from McCoy (2026). `"A2"` (crossover
  only) is a good default; `"A7"` is the oncology-crossover
  co-occurrence stress scenario.

- seed:

  Optional random seed.

## Value

A `va_switch_data` object.
