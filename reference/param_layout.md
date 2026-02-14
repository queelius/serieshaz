# Get the parameter layout for a system

Returns the mapping from global (flat) parameter indices to
per-component parameter indices, enabling the series system to
distribute a single parameter vector across its components.

## Usage

``` r
param_layout(x, ...)

# S3 method for class 'dfr_dist_series'
param_layout(x, ...)
```

## Arguments

- x:

  A system object (e.g.,
  [`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)).

- ...:

  Additional arguments passed to methods.

## Value

A list of integer vectors, one per component, containing global
parameter indices.

## Details

Parameters across all components are stored as a single concatenated
vector \\\theta = (\theta_1, \ldots, \theta_m)\\. The layout maps global
indices back to each component. For example, with:

- Component 1: Weibull (shape, scale) — 2 parameters

- Component 2: Exponential (rate) — 1 parameter

- Component 3: Gompertz (a, b) — 2 parameters

the layout is `list(1:2, 3, 4:5)`, so the global parameter vector
`c(shape1, scale1, rate2, a3, b3)` gets sliced as `par[1:2]` for
component 1, `par[3]` for component 2, and `par[4:5]` for component 3.

This design enables standard optimizers to work on a flat vector while
the series system internally distributes parameters to the correct
components.

## Methods (by class)

- `param_layout(dfr_dist_series)`: Parameter index mapping for a series
  system.

## See also

[`component`](https://queelius.github.io/serieshaz/reference/component.md)
to extract a component with its parameters,
[`params`](https://queelius.github.io/algebraic.dist/reference/params.html)
to get the full parameter vector,
[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor

Other system introspection:
[`component()`](https://queelius.github.io/serieshaz/reference/component.md),
[`component_hazard()`](https://queelius.github.io/serieshaz/reference/component_hazard.md),
[`ncomponents()`](https://queelius.github.io/serieshaz/reference/ncomponents.md),
[`sample_components()`](https://queelius.github.io/serieshaz/reference/sample_components.md)

## Examples

``` r
# \donttest{
library(flexhaz)

sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05),
    dfr_gompertz(a = 0.01, b = 0.1)
))
param_layout(sys)
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3
#> 
#> [[3]]
#> [1] 4 5
#> 
# list(1:2, 3, 4:5)
# }
```
