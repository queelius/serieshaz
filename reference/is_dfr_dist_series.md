# Test whether an object is a dfr_dist_series

Returns `TRUE` if `x` inherits from `"dfr_dist_series"`, `FALSE`
otherwise.

## Usage

``` r
is_dfr_dist_series(x)
```

## Arguments

- x:

  Object to test.

## Value

Logical scalar.

## Details

Since `dfr_dist_series` inherits from `dfr_dist`, an object that passes
`is_dfr_dist_series()` will also pass
[`is_dfr_dist()`](https://queelius.github.io/dfr.dist/reference/is_dfr_dist.html).
Use this function when you need to distinguish series systems from
ordinary `dfr_dist` objects.

## See also

[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md)
for the constructor,
[`is_dfr_dist`](https://queelius.github.io/dfr.dist/reference/is_dfr_dist.html)
for the parent class predicate

Other series system:
[`assumptions.dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/assumptions.dfr_dist_series.md),
[`dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.md),
[`print.dfr_dist_series()`](https://queelius.github.io/dfr.dist.series/reference/print.dfr_dist_series.md)

## Examples

``` r
# \donttest{
library(dfr.dist)

sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_exponential(0.2)
))
is_dfr_dist_series(sys)  # TRUE
#> [1] TRUE
is_dfr_dist(sys)         # also TRUE (inherits dfr_dist)
#> [1] TRUE

single <- dfr_exponential(0.5)
is_dfr_dist_series(single)  # FALSE
#> [1] FALSE
is_dfr_dist(single)         # TRUE
#> [1] TRUE

is_dfr_dist_series(42)  # FALSE
#> [1] FALSE
# }
```
