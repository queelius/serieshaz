# serieshaz 0.2.0

## New features

* `dfr_dist_series` objects now also inherit from `dist_structure` (from
  the `dist.structure` package). Series-system topology generics work
  out of the box: `phi`, `min_paths`, `min_cuts`, `system_signature`,
  plus the defaults from `dist.structure` (`structural_importance`,
  `reliability`, `dual`, `is_coherent`, `critical_states`,
  `system_lifetime`, `system_censoring`, `birnbaum_importance`,
  `criticality_importance`, `vesely_fussell_importance`,
  `substitute_component`, `compose_systems`).

* The local `ncomponents` and `component` generics are removed and
  re-exported from `dist.structure` instead. This eliminates a latent
  namespace collision when both `serieshaz` and `dist.structure` are
  loaded. Existing user code that calls `ncomponents(sys)` or
  `component(sys, j)` is unaffected.

## Internal

* Class chain on `dfr_dist_series` is now
  `c("dfr_dist_series", "dfr_dist", "dist_structure", "univariate_dist", "dist")`.
  The `"dist_structure"` insertion is positioned after `"dfr_dist"` so
  the optimized `flexhaz` survival, hazard, sampler, and CDF methods
  continue to dispatch first; `"dist_structure"` provides defaults only
  for generics that `flexhaz`/`serieshaz` do not already implement.

# serieshaz 0.1.1

## Improvements

* Decomposed system score and Hessian: instead of `numDeriv` on the full
  system log-likelihood, both are computed per-component. When components
  provide analytical `score_fn`/`hess_fn`, the cumulative hazard derivatives
  are extracted analytically (via the all-censored trick); the hazard
  derivatives use `numDeriv::jacobian`/`numDeriv::hessian` per-component
  (lower-dimensional). The Hessian exploits block structure: cross-component
  blocks reuse rate Jacobians from the score. Roughly m-fold speedup for
  systems with m components.

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
