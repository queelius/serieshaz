# Get the hazard function for a specific component

Returns a closure that computes the hazard rate for component `j` of a
series system. Useful for plotting hazard decompositions and
understanding each component's contribution to system risk.

## Usage

``` r
component_hazard(x, j, ...)

# S3 method for class 'dfr_dist_series'
component_hazard(x, j, ...)
```

## Arguments

- x:

  A system object (e.g.,
  [`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)).

- j:

  Component index (integer, `1 <= j <= ncomponents(x)`).

- ...:

  Additional arguments passed to methods.

## Value

A closure `function(t, par = NULL, ...)` that evaluates component `j`'s
hazard rate. If `par` is `NULL`, the component's default parameters
(from the system) are used.

## Details

The returned closure evaluates \\h_j(t, \theta_j)\\ for component `j`.
The `par` argument accepts *component-local* parameters (not the full
system parameter vector). This is useful for:

- Plotting individual hazard contributions

- Verifying that \\\sum_j h_j(t) = h\_{sys}(t)\\

- Sensitivity analysis on a single component

## Methods (by class)

- `component_hazard(dfr_dist_series)`: Hazard closure for component j of
  a series system. Returns `function(t, par_j = NULL, ...)` where
  `par_j` are the component-local parameters.

## See also

[`component`](https://queelius.github.io/serieshaz/reference/component.md)
to extract the full component object,
[`hazard`](https://queelius.github.io/algebraic.dist/reference/hazard.html)
for the system-level hazard,
[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor

Other system introspection:
[`component()`](https://queelius.github.io/serieshaz/reference/component.md),
[`ncomponents()`](https://queelius.github.io/serieshaz/reference/ncomponents.md),
[`param_layout()`](https://queelius.github.io/serieshaz/reference/param_layout.md),
[`sample_components()`](https://queelius.github.io/serieshaz/reference/sample_components.md)

## Examples

``` r
# \donttest{
library(flexhaz)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_exponential(0.2)
))

h1 <- component_hazard(sys, 1)
h2 <- component_hazard(sys, 2)
h_sys <- hazard(sys)

# Verify hazard sum property
t <- 10
h1(t) + h2(t)  # 0.3
#> [1] 0.3
h_sys(t)        # 0.3 (same!)
#> [1] 0.3
# }
```
