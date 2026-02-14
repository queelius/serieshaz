# Print method for series system distributions

Displays a human-readable summary of a series system distribution,
including the number of components, per-component parameter counts and
values, and the system hazard/survival formulas.

## Usage

``` r
# S3 method for class 'dfr_dist_series'
print(x, ...)
```

## Arguments

- x:

  A `dfr_dist_series` object.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

The output includes:

- Header with the number of components

- One line per component showing its parameter count and current
  parameter values (or "unknown" if parameters are `NULL`)

- The system hazard formula: \\h\_{sys}(t) = \sum_j h_j(t, \theta_j)\\

- The system survival formula: \\S\_{sys}(t) = \prod_j S_j(t,
  \theta_j)\\

## See also

[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor

Other series system:
[`assumptions.dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/assumptions.dfr_dist_series.md),
[`dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md),
[`is_dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/is_dfr_dist_series.md)

## Examples

``` r
# \donttest{
library(flexhaz)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_weibull(shape = 2, scale = 100)
))
print(sys)
#> Series system distribution with 2 components
#>   Component 1: 1 param(s) [0.1]
#>   Component 2: 2 param(s) [  2, 100]
#> System hazard: h_sys(t) = sum_j h_j(t, theta_j)
#> Survival: S_sys(t) = exp(-H_sys(t)) = prod_j S_j(t, theta_j)
# Series system distribution with 2 components
#   Component 1: 1 param(s) [0.1]
#   Component 2: 2 param(s) [2, 100]
# System hazard: h_sys(t) = sum_j h_j(t, theta_j)
# Survival: S_sys(t) = exp(-H_sys(t)) = prod_j S_j(t, theta_j)
# }
```
