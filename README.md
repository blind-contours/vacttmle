# vacttmle

<!-- badges: start -->
<!-- badges: end -->

**Visit-aligned continuous-time TMLE** for randomized trials in which treatment
switching (crossover) and loss to follow-up happen in *continuous time* — dated
to the day — while covariates are recorded only at *scheduled visits*.

`vacttmle` estimates failure-free survival and risk differences under the two
ICH E9(R1) strategies for treatment switching:

- **ITT / treatment-policy** — switching left as it naturally occurs, censoring
  removed: `P(T^{a, no-C} > tau)`.
- **Per-protocol / hypothetical no-switching** — switching *and* censoring
  removed: `P(T^{a, no-R, no-C} > tau)`.

It uses the exact switch and dropout dates (no discretization, no arbitrary
within-interval ordering of events), is doubly robust between the
switching-and-censoring mechanism and the outcome-regression sequence, and
provides influence-function inference, per-visit positivity diagnostics, and
percentile weight control for heavy switching. Methods: McCoy (2026).

## Installation

```r
# install.packages("remotes")
remotes::install_github("blind-contours/vacttmle")
```

Depends only on `data.table`, `survival`, and base R.

## Quick start (treatment switching)

```r
library(vacttmle)

# Simulated scheduled-visit trial with crossover + informative censoring
dat <- simulate_va_switch_trial(n = 1000, scenario = "A2", seed = 1)

# Per-protocol (hypothetical no-switching) 12-month risk difference
va_ct_switch(dat, estimand = "PP")

# Treatment-policy (ITT) from the same data
va_ct_switch(dat, estimand = "ITT")
```

### Using your own trial data

Provide one row per subject with the exact event times, plus a
subject-by-visit matrix of the visit covariate:

```r
subject_df <- data.frame(
  id, A0,             # randomized arm (0/1)
  W1, W2,             # baseline covariates
  T_time,             # failure time
  R_time,             # switch time  (Inf if never switched)
  C_time              # dropout time (Inf if never censored)
)
visit_L <- matrix(...) # n x length(visit_times): covariate L at each visit

dat <- as_va_switch_data(subject_df, visit_L,
                         visit_times = c(0, 1, 3, 6, 12), tau = 12)

# Heavy switching? Control the cumulative weight with percentile truncation:
va_ct_switch(dat, estimand = "PP", weight_trunc = "p95")
```

## Public functions

**Switching estimands (primary):**
- `va_ct_switch()` — ITT or per-protocol failure-free survival risk difference.
- `as_va_switch_data()` — build the analysis object from tidy trial data.
- `simulate_va_switch_trial()` — example scheduled-visit trials (scenarios A0–A7).
- `snap_to_grid()`, `est_discrete_ltmle()` — the discrete-time LTMLE comparator
  (snapping + explicit within-interval ordering convention) used to reproduce
  the paper's ordering-sensitivity experiment.

**Composite endpoint (also provided; SEs pending re-audit):**
- `va_ct_tmle()`, `va_visit_tmle()`, `va_ct_gcomp()` — composite-event-free
  survival `P(T^d > tau, D^d > tau)`; `as_va_data()`, `validate_va_data()`,
  `simulate_va_trial()`, `run_vacttmle_toy_checks()`. Note: the composite path
  has not yet been re-audited for the influence-function variance correction
  applied to the switching estimators (see `NEWS.md`); treat its standard
  errors as potentially conservative until then.

## Scope of this release (0.1.0)

- Baseline randomization with known `g0`; two baseline covariates (`W1`, `W2`)
  and a single scheduled-visit covariate `L`. Map your covariates onto these.
- Static regimes (assigned arm, with/without switching). The general
  switching-intensity intervention class of the paper (dynamic visit-gated and
  stochastic-scaled regimes) is prototype-validated but not yet in the public
  API.
- Piecewise-exponential working models; plug in flexible learners for
  nonparametric efficiency.

## Citation

McCoy, D. (2026). *Doubly robust, exact-time adjustment for treatment switching
in randomized trials: a visit-aligned continuous-time TMLE.* Preprint.
