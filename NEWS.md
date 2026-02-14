# serieshaz 0.1.0

Initial release.

## Features

* `dfr_dist_series()` constructor composes multiple `dfr_dist` components into a
  series system distribution where the system hazard is the sum of component
  hazards.
* Full inheritance from `dfr_dist`: hazard, survival, CDF, density, quantile,
  sampling, log-likelihood, score, Hessian, and MLE fitting all work
  automatically.
* Analytical cumulative hazard when all components provide `cum_haz_rate`;
  numerical integration fallback otherwise.
* Series-specific introspection: `ncomponents()`, `component()`,
  `param_layout()`, `component_hazard()`, `sample_components()`.
* Nested series system support (a `dfr_dist_series` can be a component in
  another series system).
* Four vignettes: overview, mathematical foundations, MLE fitting workflows,
  and advanced composition.
