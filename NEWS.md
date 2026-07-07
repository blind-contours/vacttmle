# vacttmle 0.1.0

* New primary functionality: **treatment-switching estimands**. `va_ct_switch()`
  estimates failure-free survival risk differences under the ICH E9(R1)
  treatment-policy (ITT) and hypothetical no-switching (per-protocol)
  strategies, using exact switch and dropout dates with visit-measured
  time-varying confounders.
* `as_va_switch_data()` builds the analysis object from a tidy one-row-per-
  subject table plus a subject-by-visit covariate matrix;
  `simulate_va_switch_trial()` generates example trials (scenarios A0–A7).
* Percentile weight control for heavy switching via
  `va_ct_switch(..., weight_trunc = "p95")`, with per-visit positivity
  diagnostics in the result object.
* Corrected efficient influence function: the intervened switching and
  censoring intensities enter only through the cumulative inverse weight and
  carry no score terms; SEs are calibrated (an earlier internal formulation
  including coarsening "augmentation" scores passed mean-zero checks but
  overstated variances — see the paper's gradient-pairing gate).
* Discrete-time LTMLE comparator (`snap_to_grid()`, `est_discrete_ltmle()`)
  with explicit within-interval ordering conventions, reproducing the paper's
  ordering-sensitivity experiment.
* Composite-endpoint estimators retained from 0.0.1. **Note:** the composite
  path has not yet been re-audited for the influence-function variance
  correction above; treat its SEs as potentially conservative.

# vacttmle 0.0.1

* Initial standalone package: visit-aligned continuous-time TMLE for composite
  time-to-event endpoints.
