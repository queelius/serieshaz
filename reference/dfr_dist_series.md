# Series System Distribution from DFR Components

Composes m `dfr_dist` component distributions into a series system
distribution. A series system fails when any component fails, so the
system hazard is the sum of component hazards: \\h\_{sys}(t) = \sum_j
h_j(t)\\.

## Usage

``` r
dfr_dist_series(components, par = NULL, n_par = NULL)
```

## Arguments

- components:

  A list of `dfr_dist` objects representing system components.

- par:

  Optional concatenated parameter vector \\\theta = (\theta_1, \ldots,
  \theta_m)\\. If `NULL`, parameters are concatenated from component
  objects.

- n_par:

  Optional integer vector giving the number of parameters per component.
  Inferred from component `par` if not supplied; required when any
  component has `NULL` parameters.

## Value

A `dfr_dist_series` object (inherits `dfr_dist`). Extra fields:
`$components`, `$layout`, `$m`, `$n_par`.

## Details

The resulting object inherits from `dfr_dist`, so all existing methods
(hazard, survival, CDF, density, sampling, log-likelihood, MLE fitting)
work automatically.

**Parameter layout**: Parameters are stored as a single concatenated
vector. The `$layout` field maps global indices to component indices.
For example, if component 1 has 2 parameters and component 2 has 1, then
`layout = list(1:2, 3)`.

**Analytical cumulative hazard**: If *all* components provide
`cum_haz_rate`, the series system gets an analytical \\H\_{sys}(t) =
\sum_j H_j(t)\\. Otherwise, falls back to numerical integration.

**Score and Hessian**: Fall back to
[`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) and
[`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html) on
the (correct) composed log-likelihood.

**Identifiability**: Exponential series systems are *not* identifiable
from system-level data alone — only the sum of rates is identifiable.
When fitting to data, check `sum(coef(result))` rather than individual
rate parameters. Mixed-type series systems (e.g., Weibull + Gompertz)
are generally identifiable because the components have different hazard
shapes.

**Nested series**: A `dfr_dist_series` is itself a `dfr_dist`, so it can
be used as a component in another series system. The resulting nested
system's hazard is the sum of all leaf-component hazards.

**Class hierarchy**: `dfr_dist_series` inherits from `dfr_dist` -\>
`likelihood_model` -\> `univariate_dist` -\> `dist`. All methods from
these parent classes work automatically.

## See also

[`is_dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/is_dfr_dist_series.md)
for the type predicate,
[`ncomponents`](https://queelius.github.io/dfr.dist.series/reference/ncomponents.md)
and
[`component`](https://queelius.github.io/dfr.dist.series/reference/component.md)
for introspection,
[`param_layout`](https://queelius.github.io/dfr.dist.series/reference/param_layout.md)
for parameter index mapping,
[`component_hazard`](https://queelius.github.io/dfr.dist.series/reference/component_hazard.md)
for per-component hazard closures,
[`sample_components`](https://queelius.github.io/dfr.dist.series/reference/sample_components.md)
for sampling component lifetimes,
[`dfr_dist`](https://queelius.github.io/dfr.dist/reference/dfr_dist.html)
for the parent class constructor,
[`hazard`](https://queelius.github.io/algebraic.dist/reference/hazard.html)
for distribution generics

Other series system:
[`assumptions.dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/assumptions.dfr_dist_series.md),
[`is_dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/is_dfr_dist_series.md),
[`print.dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/print.dfr_dist_series.md)

## Examples

``` r
# \donttest{
library(dfr.dist)

# --- Basic exponential series ---
# Three exponential components -> equivalent to single exponential
sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_exponential(0.2),
    dfr_exponential(0.3)
))
# System hazard = 0.6 (constant)
h <- hazard(sys)
h(10)  # 0.6
#> [1] 0.6

# System survival at t = 5
S <- surv(sys)
S(5)   # exp(-0.6 * 5)
#> [1] 0.04978707

# --- Mixed Weibull + Gompertz series ---
sys2 <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_gompertz(a = 0.01, b = 0.1)
))
h2 <- hazard(sys2)
h2(50)  # sum of Weibull and Gompertz hazards at t=50
#> [1] 1.494132

# --- Nested series ---
subsystem <- dfr_dist_series(list(
    dfr_exponential(0.05),
    dfr_exponential(0.10)
))
full_system <- dfr_dist_series(list(
    subsystem,
    dfr_weibull(shape = 2, scale = 200)
))

# --- Fitting workflow ---
solver <- fit(sys)
# result <- solver(df, par = c(0.1, 0.2, 0.3))
# coef(result)   # fitted parameters
# vcov(result)   # variance-covariance matrix
# logLik(result) # maximized log-likelihood
# }
```
