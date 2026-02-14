# serieshaz: Series System Distributions from Dynamic Failure Rate Components

Compose multiple dynamic failure rate (DFR) distributions into series
system distributions. A series system fails when any component fails,
giving a system hazard equal to the sum of component hazards. The
resulting object inherits from 'dfr_dist' so all existing methods
(hazard, survival, CDF, density, sampling, MLE fitting) work
automatically.

## Details

The serieshaz package composes multiple dynamic failure rate
([`dfr_dist`](https://queelius.github.io/flexhaz/reference/dfr_dist.html))
distributions into a series system distribution. A series system fails
when *any* component fails, so the system hazard is the sum of component
hazards: \$\$h\_{sys}(t) = \sum\_{j=1}^{m} h_j(t, \theta_j)\$\$ and the
system survival is the product of component survivals: \$\$S\_{sys}(t) =
\prod\_{j=1}^{m} S_j(t, \theta_j)\$\$

The series system object inherits from `dfr_dist`, which in turn
inherits from `likelihood_model`, `univariate_dist`, and `dist`. This
means all existing methods — hazard, survival, CDF, density, quantile
function, sampling, log-likelihood, score, Hessian, and MLE fitting —
work automatically on series systems without reimplementation.

Parameters across all components are stored as a single flat vector,
with a *layout* that maps global indices to per-component indices. This
design enables standard optimizers (e.g.,
[`optim`](https://rdrr.io/r/stats/optim.html)) to work directly on the
concatenated parameter vector.

## Package functions

- [`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md):

  Constructor: compose components into a series system

- [`is_dfr_dist_series`](https://queelius.github.io/serieshaz/reference/is_dfr_dist_series.md):

  Type predicate

- [`ncomponents`](https://queelius.github.io/serieshaz/reference/ncomponents.md):

  Number of components

- [`component`](https://queelius.github.io/serieshaz/reference/component.md):

  Extract a single component

- [`param_layout`](https://queelius.github.io/serieshaz/reference/param_layout.md):

  Parameter index mapping

- [`component_hazard`](https://queelius.github.io/serieshaz/reference/component_hazard.md):

  Component-level hazard closure

- [`sample_components`](https://queelius.github.io/serieshaz/reference/sample_components.md):

  Sample component lifetimes

## See also

[`dfr_dist_series`](https://queelius.github.io/serieshaz/reference/dfr_dist_series.md)
for the constructor,
[`dfr_dist`](https://queelius.github.io/flexhaz/reference/dfr_dist.html)
for the parent class,
[`hazard`](https://queelius.github.io/algebraic.dist/reference/hazard.html)
for distribution generics,
[`loglik`](https://queelius.github.io/likelihood.model/reference/loglik.html)
for statistical inference generics

[`vignette("series-overview")`](https://queelius.github.io/serieshaz/articles/series-overview.md)
for a quick-start guide,
[`vignette("series-math")`](https://queelius.github.io/serieshaz/articles/series-math.md)
for mathematical foundations,
[`vignette("series-fitting")`](https://queelius.github.io/serieshaz/articles/series-fitting.md)
for MLE fitting workflows

## Author

**Maintainer**: Alexander Towell <lex@metafunctor.com>
