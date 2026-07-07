# Package index

## Package

- [`vacttmle-package`](https://blind-contours.github.io/vacttmle/reference/vacttmle-package.md)
  [`vacttmle`](https://blind-contours.github.io/vacttmle/reference/vacttmle-package.md)
  : Visit-aligned continuous-time TMLE

## Treatment-switching estimands (primary)

- [`va_ct_switch()`](https://blind-contours.github.io/vacttmle/reference/va_ct_switch.md)
  : Estimate ITT or per-protocol failure-free survival with VA-CT-TMLE
- [`as_va_switch_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_switch_data.md)
  : Build a visit-aligned analysis object for the switching estimands
- [`simulate_va_switch_trial()`](https://blind-contours.github.io/vacttmle/reference/simulate_va_switch_trial.md)
  : Simulate a scheduled-visit trial with switching and censoring

## Discrete-LTMLE comparator

- [`snap_to_grid()`](https://blind-contours.github.io/vacttmle/reference/snap_to_grid.md)
  : Build a snapped, discretized obs_data under a given ordering
  convention.
- [`est_discrete_ltmle()`](https://blind-contours.github.io/vacttmle/reference/est_discrete_ltmle.md)
  : Discrete-LTMLE PP/ITT contrast under a fixed ordering convention and
  grid.

## Composite-endpoint estimators

- [`va_ct_tmle()`](https://blind-contours.github.io/vacttmle/reference/va_ct_tmle.md)
  [`va_visit_tmle()`](https://blind-contours.github.io/vacttmle/reference/va_ct_tmle.md)
  [`va_ct_gcomp()`](https://blind-contours.github.io/vacttmle/reference/va_ct_tmle.md)
  : Estimate a composite risk difference
- [`va_visit_tmle()`](https://blind-contours.github.io/vacttmle/reference/va_visit_tmle.md)
  : Estimate a composite risk difference with visit-node TMLE
- [`va_ct_gcomp()`](https://blind-contours.github.io/vacttmle/reference/va_ct_gcomp.md)
  : Estimate a composite risk difference with VA-CT g-computation

## Data preparation

- [`as_va_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_data.md)
  : Build a visit-aligned analysis object
- [`validate_va_data()`](https://blind-contours.github.io/vacttmle/reference/validate_va_data.md)
  : Validate visit-aligned data

## Examples and checks

- [`simulate_va_trial()`](https://blind-contours.github.io/vacttmle/reference/simulate_va_trial.md)
  : Simulate a visit-aligned trial
- [`run_vacttmle_toy_checks()`](https://blind-contours.github.io/vacttmle/reference/run_vacttmle_toy_checks.md)
  : Run package toy checks
