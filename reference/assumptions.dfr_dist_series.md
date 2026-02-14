# Assumptions for series system distributions

Returns the statistical and structural assumptions underlying a series
system model, which are important for the validity of MLE-based
inference.

## Usage

``` r
# S3 method for class 'dfr_dist_series'
assumptions(model, ...)
```

## Arguments

- model:

  A `dfr_dist_series` object.

- ...:

  Additional arguments (unused).

## Value

Character vector of model assumptions.

## Details

The assumptions returned are:

- **Series structure**: The system fails when any component fails
  (weakest-link model)

- **Component independence**: Component lifetimes are statistically
  independent

- **Non-negative hazard**: Each component hazard satisfies \\h_j(t) \geq
  0\\ for all \\t \> 0\\

- **Proper distribution**: The cumulative hazard diverges, ensuring
  \\S\_{sys}(t) \to 0\\ as \\t \to \infty\\

- **Positive support**: The time domain is \\(0, \infty)\\

- **Independent observations**: The observed lifetimes are independent

- **Censoring convention**: `delta = 1` for exact, `0` for
  right-censored, `-1` for left-censored

- **Non-informative censoring**: The censoring mechanism carries no
  information about the failure process

These assumptions are required for the MLE fitting procedure
([`fit`](https://generics.r-lib.org/reference/fit.html)) to produce
valid estimates. Violation of component independence, in particular,
invalidates the hazard-sum property that defines series systems.

## See also

[`assumptions`](https://queelius.github.io/likelihood.model/reference/assumptions.html)
for the generic,
[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor,
[`vignette("series-fitting")`](https://queelius.github.io/serieshaz/articles/series-fitting.md)
for how assumptions affect inference

Other series system:
[`dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md),
[`is_dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/is_dfr_dist_series.md),
[`print.dfr_dist_series()`](https://queelius.github.io/serieshaz/reference/print.dfr_dist_series.md)

## Examples

``` r
# \donttest{
library(flexhaz)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_weibull(shape = 2, scale = 100)
))
assumptions(sys)
#> [1] "Series system: system fails when any component fails"            
#> [2] "Component independence: component lifetimes are independent"     
#> [3] "Non-negative hazard: h_j(t) >= 0 for all j, t > 0"               
#> [4] "Cumulative hazard diverges: lim(t->Inf) H_sys(t) = Inf"          
#> [5] "Support is positive reals: t in (0, Inf)"                        
#> [6] "Observations are independent"                                    
#> [7] "Censoring indicator: 1=exact, 0=right-censored, -1=left-censored"
#> [8] "Non-informative censoring"                                       
# }
```
