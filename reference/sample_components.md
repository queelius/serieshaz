# Sample component lifetimes from a system

Generates an \\n \times m\\ matrix where column \\j\\ contains
independent samples from component \\j\\'s lifetime distribution. The
system lifetime is the row-wise minimum.

## Usage

``` r
sample_components(x, n, ...)

# S3 method for class 'dfr_dist_series'
sample_components(x, n, par = NULL, ...)
```

## Arguments

- x:

  A system object (e.g.,
  [`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)).

- n:

  Number of samples (rows).

- ...:

  Additional arguments passed to methods.

- par:

  Optional parameter vector override.

## Value

An \\n \times m\\ numeric matrix of component lifetimes, with columns
named `comp1`, `comp2`, etc.

## Details

Each column is sampled independently using the component's own sampler.
Since the series system fails when *any* component fails, the system
lifetime for each observation is:

    t_sys <- apply(mat, 1, min)

The failing component for each observation can be identified via:

    failing <- apply(mat, 1, which.min)

This enables failure attribution analysis: what proportion of system
failures are caused by each component?

## Methods (by class)

- `sample_components(dfr_dist_series)`: Sample component lifetimes from
  a series system. Returns an n x m matrix where column j holds samples
  from component j. The system lifetime is `apply(mat, 1, min)`.

## See also

[`sampler`](https://queelius.github.io/algebraic.dist/reference/sampler.html)
for system-level sampling,
[`component`](https://queelius.github.io/serieshaz/reference/component.md)
to extract individual component objects,
[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor

Other system introspection:
[`component()`](https://queelius.github.io/serieshaz/reference/component.md),
[`component_hazard()`](https://queelius.github.io/serieshaz/reference/component_hazard.md),
[`ncomponents()`](https://queelius.github.io/serieshaz/reference/ncomponents.md),
[`param_layout()`](https://queelius.github.io/serieshaz/reference/param_layout.md)

## Examples

``` r
# \donttest{
library(flexhaz)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_exponential(0.2),
    dfr_exponential(0.3)
))

set.seed(42)
mat <- sample_components(sys, n = 1000)
dim(mat)  # 1000 x 3
#> [1] 1000    3

# System lifetimes
t_sys <- apply(mat, 1, min)

# Which component caused each failure?
failing <- apply(mat, 1, which.min)
table(failing) / 1000
#> failing
#>     1     2     3 
#> 0.173 0.345 0.482 
# Proportions ~= c(1/6, 2/6, 3/6) for rates (0.1, 0.2, 0.3)
# }
```
