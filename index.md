# dfr.dist.series

Series System Distributions from Dynamic Failure Rate Components

**dfr.dist.series** composes multiple `dfr_dist` objects into a series
system distribution. A series system fails when *any* component fails,
so the system hazard is the sum of component hazards:

$$h_{sys}(t) = \sum\limits_{j = 1}^{m}h_{j}\left( t,\theta_{j} \right)$$

The resulting object inherits from `dfr_dist`, so all existing methods —
hazard, survival, CDF, density, sampling, log-likelihood, and MLE
fitting — work automatically.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("queelius/dfr.dist.series")
```

## Quick Start

``` r
library(dfr.dist.series)

# Three-component server with different failure modes
server <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 500),     # disk wear-out
    dfr_exponential(0.001),            # random memory failure
    dfr_gompertz(a = 0.0001, b = 0.02)       # PSU degradation
))

# Evaluate system hazard and survival
h <- hazard(server)
S <- surv(server)

h(100)   # system hazard at t = 100
#> [1] 0.002538906
S(100)   # probability of surviving past t = 100
#> [1] 0.8420252
```

``` r
# Sample system lifetimes
set.seed(42)
samp <- sampler(server)
times <- samp(5)
times
#> [1] 294.9884 302.0998 150.3888 274.2224 236.8369
```

``` r
# Introspect: which component contributes most at t = 200?
for (j in 1:ncomponents(server)) {
    hj <- component_hazard(server, j)
    cat(sprintf("Component %d hazard at t=200: %.6f\n", j, hj(200)))
}
#> Component 1 hazard at t=200: 0.001600
#> Component 2 hazard at t=200: 0.001000
#> Component 3 hazard at t=200: 0.005460
```

## Key Features

- **Composition**: Combine any `dfr_dist` objects (Weibull, exponential,
  Gompertz, log-logistic, custom) into series systems
- **Full interface**: All distribution methods (hazard, survival, CDF,
  density, quantile, sampling) work out of the box
- **MLE fitting**: Fit series system parameters to observed failure data
  with [`fit()`](https://generics.r-lib.org/reference/fit.html)
- **Introspection**:
  [`ncomponents()`](https://queelius.github.io/dfr.dist.series/reference/ncomponents.md),
  [`component()`](https://queelius.github.io/dfr.dist.series/reference/component.md),
  [`param_layout()`](https://queelius.github.io/dfr.dist.series/reference/param_layout.md),
  [`component_hazard()`](https://queelius.github.io/dfr.dist.series/reference/component_hazard.md),
  [`sample_components()`](https://queelius.github.io/dfr.dist.series/reference/sample_components.md)
- **Nesting**: Series systems can be nested as components of larger
  series systems
- **Analytical cumulative hazard**: When all components provide
  closed-form cumulative hazard, the series system does too

## Ecosystem

dfr.dist.series builds on:

- [algebraic.dist](https://github.com/queelius/algebraic.dist) — Base
  distribution generics
- [likelihood.model](https://github.com/queelius/likelihood.model) —
  Statistical inference generics
- [dfr.dist](https://github.com/queelius/dfr.dist) — Dynamic failure
  rate distributions

## Documentation

- [`vignette("series-overview")`](https://queelius.github.io/dfr.dist.series/articles/series-overview.md)
  — Package overview and quick start
- [`vignette("series-math")`](https://queelius.github.io/dfr.dist.series/articles/series-math.md)
  — Mathematical foundations
- [`vignette("series-fitting")`](https://queelius.github.io/dfr.dist.series/articles/series-fitting.md)
  — MLE fitting and inference
- [`vignette("series-advanced")`](https://queelius.github.io/dfr.dist.series/articles/series-advanced.md)
  — Advanced composition patterns
