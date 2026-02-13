# Get the number of components in a system

Returns the number of components \\m\\ in a series system, corresponding
to the number of terms in \\h\_{sys}(t) = \sum\_{j=1}^{m} h_j(t)\\.

## Usage

``` r
ncomponents(x, ...)

# S3 method for class 'dfr_dist_series'
ncomponents(x, ...)
```

## Arguments

- x:

  A system object (e.g.,
  [`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)).

- ...:

  Additional arguments passed to methods.

## Value

Integer, the number of components.

## Details

For a
[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)
object created from a list of `m` components, this simply returns `m`.
This is useful for programmatically iterating over components, e.g., for
plotting hazard decompositions or computing failure attribution.

## Methods (by class)

- `ncomponents(dfr_dist_series)`: Number of components in a series
  system.

## See also

[`component`](https://queelius.github.io/dfr.dist.series/reference/component.md)
to extract individual components,
[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)
for the constructor

Other system introspection:
[`component()`](https://queelius.github.io/dfr.dist.series/reference/component.md),
[`component_hazard()`](https://queelius.github.io/dfr.dist.series/reference/component_hazard.md),
[`param_layout()`](https://queelius.github.io/dfr.dist.series/reference/param_layout.md),
[`sample_components()`](https://queelius.github.io/dfr.dist.series/reference/sample_components.md)

## Examples

``` r
# \donttest{
library(dfr.dist)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_weibull(shape = 2, scale = 100),
    dfr_gompertz(a = 0.01, b = 0.05)
))
ncomponents(sys)  # 3
#> [1] 3
# }
```
