# Build a snapped, discretized obs_data under a given ordering convention.

Build a snapped, discretized obs_data under a given ordering convention.

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

  original continuous-time obs_data (subject + person_interval)

- delta:

  bin width (months)

- ordering:

  "event_first" \| "switch_first" \| "drop_first"

- tau:

  horizon
