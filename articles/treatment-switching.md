# ITT and per-protocol estimands under treatment switching

``` r
library(vacttmle)
```

## The problem

In many randomized trials the treatment assigned at baseline is not the
treatment received throughout follow-up: control-arm patients cross over
to the experimental treatment, patients switch off toxic regimens, and
others are lost to follow-up. These events happen in continuous time and
are recorded to the day, while the covariates that drive them (disease
status, biomarkers) are measured only at scheduled visits.

Two estimands answer different questions, both mapped to ICH E9(R1)
strategies:

- **Treatment-policy (ITT):** contrast the randomized arms leaving
  switching as it naturally occurs (censoring removed).
- **Hypothetical no-switching (per-protocol, PP):** contrast “assigned
  to arm *a* and never switch” (switching and censoring removed).

`vacttmle` estimates both from the exact switch and dropout dates,
adjusting for the visit-measured time-varying confounders that drive
switching and dropout.

## A simulated trial

``` r
dat <- simulate_va_switch_trial(n = 1000, scenario = "A2", seed = 1)
```

Scenario `"A2"` has crossover driven by the visit covariate and no
dropout; `"A1"` adds informative censoring, `"A7"` is a heavy-switching
oncology-style co-occurrence scenario.

## Estimating both estimands

``` r
pp  <- va_ct_switch(dat, estimand = "PP",  seed = 1)
itt <- va_ct_switch(dat, estimand = "ITT", seed = 1)

pp
#> Visit-aligned continuous-time estimator
#> Endpoint: switching
#> Estimand: per-protocol / hypothetical no-switching
#> Estimator: va_ct_tmle
#> Risk(d0): 0.3321
#> Risk(d1): 0.2073
#> Risk difference: -0.1248
#> SE: 0.0279
#> 95% CI: [-0.1796, -0.0701]
#> Converged: yes
itt
#> Visit-aligned continuous-time estimator
#> Endpoint: switching
#> Estimand: ITT / treatment-policy (natural switching)
#> Estimator: va_ct_tmle
#> Risk(d0): 0.3084
#> Risk(d1): 0.2396
#> Risk difference: -0.0689
#> SE: 0.0257
#> 95% CI: [-0.1193, -0.0185]
#> Converged: yes
```

The per-protocol effect is typically larger than the intention-to-treat
effect: removing switching un-dilutes the contrast, since control-arm
crossover to the more effective regimen narrows the observed ITT gap.

``` r
c(ITT = itt$risk_difference, PP = pp$risk_difference)
#>         ITT          PP 
#> -0.06889338 -0.12484571
```

Inference is influence-function based:

``` r
pp$ci
#>       lower       upper 
#> -0.17955336 -0.07013807
pp$diagnostics
#>    converged n_outer_iterations max_abs_total_eif max_weight_d0 max_weight_d1
#>       <lgcl>              <int>             <num>         <num>         <num>
#> 1:      TRUE                  1       0.005759134      2.562163      2.562163
#>      ess_d0   ess_d1
#>       <num>    <num>
#> 1: 1847.766 1789.098
```

## Your own data

Supply one row per subject with exact times (`Inf` where an event never
happens) and a subject-by-visit covariate matrix:

``` r
dat <- as_va_switch_data(
  subject_data = my_subjects,  # id, A0, W1, W2, T_time, R_time, C_time
  visit_L      = my_visit_L,   # n x length(visit_times)
  visit_times  = c(0, 1, 3, 6, 12),
  tau          = 12
)
va_ct_switch(dat, estimand = "PP")
```

## Heavy switching and positivity

When switching is frequent and strongly covariate-driven, the cumulative
inverse-probability weight can grow large. Percentile truncation on the
correctly specified mechanism controls it while retaining the covariate
dependence:

``` r
heavy <- simulate_va_switch_trial(n = 1000, scenario = "A7", seed = 2)
va_ct_switch(heavy, estimand = "PP", weight_trunc = "p95", seed = 2)
#> Visit-aligned continuous-time estimator
#> Endpoint: switching
#> Estimand: per-protocol / hypothetical no-switching
#> Estimator: va_ct_tmle
#> Risk(d0): 0.5129
#> Risk(d1): 0.1786
#> Risk difference: -0.3344
#> SE: 0.0459
#> 95% CI: [-0.4243, -0.2445]
#> Converged: yes
```

The per-visit weight maxima and effective sample sizes are in
`fit$diagnostics`, supporting the positivity checks that
switching-adjustment guidance recommends.
