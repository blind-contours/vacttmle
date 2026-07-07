# Discrete-time LTMLE contrast under a fixed grid and ordering convention

Runs a genuine discrete-time ICE-LTMLE on the output of
[`snap_to_grid()`](https://blind-contours.github.io/vacttmle/reference/snap_to_grid.md):
event/switch/censoring times live on the bin boundaries, so this is the
estimator a discrete-time analyst would compute, including its
dependence on the arbitrary within-bin ordering convention. Used to
reproduce the ordering-sensitivity experiment of the accompanying paper.

## Usage

``` r
est_discrete_ltmle(
  obs_data,
  estimand = "PP",
  tau = obs_data$tau,
  g0 = 0.5,
  delta = 1,
  ordering = "switch_first",
  M = 50,
  B = 1,
  seed = NULL,
  true_rd = NA_real_,
  scenario = NA_character_,
  replicate_id = NA_integer_
)
```

## Arguments

- obs_data:

  A `va_switch_data` object (continuous-time).

- estimand:

  `"ITT"` or `"PP"`.

- tau:

  Analysis horizon.

- g0:

  Known randomization probability of arm 1.

- delta:

  Bin width passed to
  [`snap_to_grid()`](https://blind-contours.github.io/vacttmle/reference/snap_to_grid.md).

- ordering:

  Within-bin ordering convention passed to
  [`snap_to_grid()`](https://blind-contours.github.io/vacttmle/reference/snap_to_grid.md).

- M:

  Monte Carlo draws for visit-level integration.

- B:

  Sub-bins (kept at 1 for a genuine discrete-time fit).

- seed:

  Random seed.

- true_rd:

  Optional true risk difference, carried into the result.

- scenario, replicate_id:

  Optional bookkeeping labels.

## Value

A one-row `data.table` with the estimate, SE, CI, and metadata.

## See also

[`snap_to_grid()`](https://blind-contours.github.io/vacttmle/reference/snap_to_grid.md),
[`va_ct_switch()`](https://blind-contours.github.io/vacttmle/reference/va_ct_switch.md)
