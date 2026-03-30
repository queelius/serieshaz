# serieshaz 0.1.1

## Improvements

* Decomposed system score function: instead of `numDeriv::grad` on the full
  system log-likelihood, the score is now computed per-component. When
  components provide analytical `score_fn`, the cumulative hazard derivative
  is extracted analytically (via the all-censored trick); the hazard derivative
  uses `numDeriv::jacobian` per-component (lower-dimensional). This yields
  roughly m-fold speedup for systems with m components.

* Input validation for `component()` and `component_hazard()`: the `j`
  parameter now rejects non-integer, NA, and vector inputs (previously `j = 1.5`
  silently truncated to `j = 1`).

* `component_hazard()` closure now reads the system parameter vector lazily,
  so it reflects updates to `sys$par` after construction (e.g., after fitting).

* Guard against zero-parameter components: the constructor now rejects
  `n_par` entries less than 1.

* Aligned `sys_cum_haz` initialization to `numeric(length(t))` for explicit
  vectorization consistency with `sys_rate`.

* Added `numDeriv` to Imports (previously an undeclared transitive dependency
  through `flexhaz`).

* Removed `Remotes` field from DESCRIPTION (CRAN submission blocker).

## Bug Fixes

* Fixed stale parameter capture in `component_hazard()` closure: default
  parameters were eagerly captured at closure-creation time, causing the
  closure to use outdated values if `sys$par` was updated post-construction.

## Tests

* Added tests for left-censored data (`delta = -1`), boundary conditions
  (`t = 0`), `j` parameter validation, `sample_components` error path,
  and decomposed score accuracy.

## Known Limitations

* Exponential series systems are not identifiable from system-level data alone
  (only the sum of rates is identifiable).
* The Hessian still falls back to `numDeriv::hessian` on the full system
  log-likelihood. A decomposed Hessian using the block structure
  (within-component + cross-component terms) is a future optimization.

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
