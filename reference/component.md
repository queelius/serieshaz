# Extract a component from a system

Extracts component `j` from a series system as a standalone
[`dfr_dist`](https://queelius.github.io/dfr.dist/reference/dfr_dist.html)
object, with its parameters set to the current values from the system's
parameter vector (via the layout).

## Usage

``` r
component(x, j, ...)

# S3 method for class 'dfr_dist_series'
component(x, j, ...)
```

## Arguments

- x:

  A system object (e.g.,
  [`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)).

- j:

  Component index (integer, `1 <= j <= ncomponents(x)`).

- ...:

  Additional arguments passed to methods.

## Value

A
[`dfr_dist`](https://queelius.github.io/dfr.dist/reference/dfr_dist.html)
object for component `j`.

## Details

The returned component object is a copy of the original component with
its `par` field updated to reflect the current system-level parameter
vector. This means you can evaluate the extracted component's hazard,
survival, etc. directly:

    comp1 <- component(sys, 1)
    h1 <- hazard(comp1)
    h1(10)  # evaluates using parameters from the system

Changes to the extracted component do *not* propagate back to the
original series system.

## Methods (by class)

- `component(dfr_dist_series)`: Extract component j as a standalone
  `dfr_dist` with its current parameters from the series system's
  parameter vector.

## See also

[`ncomponents`](https://queelius.github.io/dfr.dist.series/reference/ncomponents.md)
for the component count,
[`component_hazard`](https://queelius.github.io/dfr.dist.series/reference/component_hazard.md)
for getting just the hazard closure,
[`param_layout`](https://queelius.github.io/dfr.dist.series/reference/param_layout.md)
for parameter index mapping,
[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)
for the constructor

Other system introspection:
[`component_hazard()`](https://queelius.github.io/dfr.dist.series/reference/component_hazard.md),
[`ncomponents()`](https://queelius.github.io/dfr.dist.series/reference/ncomponents.md),
[`param_layout()`](https://queelius.github.io/dfr.dist.series/reference/param_layout.md),
[`sample_components()`](https://queelius.github.io/dfr.dist.series/reference/sample_components.md)

## Examples

``` r
# \donttest{
library(dfr.dist)

sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
))

# Extract the Weibull component
wb <- component(sys, 1)
params(wb)  # c(2, 100)
#> [1]   2 100

# Evaluate its hazard independently
h_wb <- hazard(wb)
h_wb(50)
#> [1] 0.01
# }
```
