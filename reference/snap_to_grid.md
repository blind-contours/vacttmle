# Snap a continuous-time trial to a discrete grid under an ordering convention

Converts a `va_switch_data` object to the discretized data a
discrete-time LTMLE analyst would construct: every
failure/switch/censoring time is snapped to its bin boundary, and
same-bin co-occurrences are resolved by the chosen within-interval
ordering convention.

## Usage

``` r
snap_to_grid(
  obs_data,
  delta,
  ordering = c("switch_first", "event_first", "drop_first"),
  tau = obs_data$tau
)
```

## Arguments

- obs_data:

  A `va_switch_data` object (continuous-time; see
  [`as_va_switch_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_switch_data.md)).

- delta:

  Bin width, on the time scale of the data.

- ordering:

  Within-bin ordering convention: `"switch_first"` (switch resolved
  before the outcome; a same-bin death is censored at switch),
  `"event_first"` (the death is adjudicated first and counted), or
  `"drop_first"` (dropout resolved before the outcome).

- tau:

  Analysis horizon.

## Value

A discretized `va_switch_data`-style object on the `delta` grid.

## See also

[`est_discrete_ltmle()`](https://blind-contours.github.io/vacttmle/reference/est_discrete_ltmle.md)
