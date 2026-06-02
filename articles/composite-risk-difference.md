# Composite Risk Difference with VA-CT-TMLE

``` r
library(vacttmle)
```

This package estimates a visit-aligned continuous-time TMLE for the
composite survival target

``` text
P(T^d > tau, D^d > tau).
```

The reported risk difference is

``` text
{1 - P(T^d1 > tau, D^d1 > tau)} -
{1 - P(T^d0 > tau, D^d0 > tau)}.
```

The intercurrent event `D` is part of the composite endpoint in this
package version. It is not a censoring process.

``` r
dat <- simulate_va_trial(n = 500, scenario = "S1", seed = 1)
fit <- va_ct_tmle(dat, M = 20, B = 10, max_iter = 5, seed = 2)
fit
```

    ## Visit-aligned continuous-time estimator
    ## Endpoint: composite
    ## Estimator: ct_tmle
    ## Risk(d0): 0.2875
    ## Risk(d1): 0.1556
    ## Risk difference: -0.1319
    ## SE: 0.0347
    ## 95% CI: [-0.1999, -0.0638]
    ## Converged: yes

Diagnostics include component-wise EIF means and cumulative-weight
summaries.

``` r
fit$diagnostics
```

    ##    converged n_outer_iterations max_abs_eif_T max_abs_eif_D max_abs_eif_Q
    ##       <lgcl>              <int>         <num>         <num>         <num>
    ## 1:      TRUE                  1    0.00037049  0.0004238737   0.001348733
    ##    max_abs_total_eif max_weight_d0 max_weight_d1 ess_d0 ess_d1
    ##                <num>         <num>         <num>  <num>  <num>
    ## 1:       0.001094326             2             2    952    971

For trial data, use
[`as_va_data()`](https://blind-contours.github.io/vacttmle/reference/as_va_data.md)
to adapt subject-level and person-interval tables, then call
[`validate_va_data()`](https://blind-contours.github.io/vacttmle/reference/validate_va_data.md)
before estimation.

See
[`vignette("trialist-quickstart", package = "vacttmle")`](https://blind-contours.github.io/vacttmle/articles/trialist-quickstart.md)
for the full trialist workflow and data-layout checklist.
