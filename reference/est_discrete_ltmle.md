# Discrete-LTMLE PP/ITT contrast under a fixed ordering convention and grid.

Discrete-LTMLE PP/ITT contrast under a fixed ordering convention and
grid.

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
