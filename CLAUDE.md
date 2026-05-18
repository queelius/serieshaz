# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Where this fits in the rlang monorepo

This package is one node in the `~/github/rlang/` reliability ecosystem. The parent
`../CLAUDE.md` documents cross-package architecture (dependency graph,
`dist.structure` protocol design, class-hierarchy invariants, `:F` vs `:G`
conventions in `kofn`, the masked-data protocol, release sequencing) and is the
authoritative reference for anything that spans packages. This file is
intentionally **package-local**: file layout, dispatch ordering inside the
class chain, derivative-decomposition tricks, identifiability edges, and test
fixture conventions specific to serieshaz.

CRAN status snapshot: 0.1.1 is on CRAN; the local tree is 0.2.0. With
`dist.structure` 0.5.0 now on CRAN (landed 2026-05-07), v0.2.0 is no
longer blocked and is being submitted in the next wave alongside
`flexhaz` 0.5.2 and `maskedcauses` 0.10.0.

## Build and Test Commands

```bash
# Install dependencies (algebraic.dist / likelihood.model / flexhaz /
# dist.structure all live in github::queelius/...; r-universe also serves them)
Rscript -e 'remotes::install_deps(dependencies = TRUE)'

# Build and check
R CMD build .
R CMD check serieshaz_*.tar.gz

# Tests
Rscript -e 'testthat::test_dir("tests/testthat")'                              # all
Rscript -e 'testthat::test_file("tests/testthat/test-constructor.R")'          # one file
Rscript -e 'testthat::test_dir("tests/testthat", filter = "dist-structure")'   # pattern
Rscript -e 'testthat::test_dir("tests/testthat", filter = "score-decomposed")' # pattern

# Regenerate NAMESPACE and man/ pages (after roxygen edits)
Rscript -e 'roxygen2::roxygenise()'

# Coverage
Rscript -e 'covr::package_coverage()'
```

When working interactively, `cd serieshaz/` then start R; running
`devtools::load_all()` from `rlang/` fails because rlang is not itself an
R package.

## Mathematical model

A series system fails when *any* component fails, so for components
`{1, ..., m}` with hazards `h_j(t, theta_j)`:

```
h_sys(t)   = sum_j h_j(t, theta_j)
S_sys(t)   = prod_j S_j(t, theta_j)
H_sys(t)   = sum_j H_j(t, theta_j)   (only when all components have closed-form cum_haz_rate)
T_sys      = min_j T_j               (component lifetimes are independent)
```

Everything in this package is a direct mechanical consequence of those four
identities.

## Class chain and dispatch ordering

The class vector on a `dfr_dist_series` object is:

```
c("dfr_dist_series", "dfr_dist", "dist_structure", "univariate_dist", "dist")
```

The ordering encodes a dispatch contract:

1. **`dfr_dist_series` (head)**: package-specific specializations (print,
   assumptions, the four `dist.structure` topology methods, the introspection
   generics).
2. **`dfr_dist`**: parent class from `flexhaz`. Its `surv`, `cdf`, `sampler`,
   `hazard`, default `loglik` / `score` / `hess_loglik` methods are reused as-is.
3. **`dist_structure`**: inserted *after* `dfr_dist`, not at the head. This is
   deliberate. It means `dist.structure` defaults provide topology generics
   that nobody else implements (`reliability`, `dual`, `structural_importance`,
   `critical_states`, Birnbaum / criticality / Vesely-Fussell importance,
   `substitute_component`, `compose_systems`, `system_lifetime`,
   `system_censoring`, `is_coherent`) but `flexhaz`'s optimized
   distribution methods win when both define a method (`surv`, `cdf`, `hazard`,
   `sampler`). See `R/dfr_dist_series.R:425-427` for the exact insertion.
4. **`univariate_dist` -> `dist`**: from `algebraic.dist`, untouched.

**If you reorder this chain, distribution methods will silently regress to
`dist_structure` defaults that compose component dists generically through the
topology.** That is correct mathematically but throws away the analytical
sum-of-hazards shortcut that motivates this package.

## Three protocols this package participates in

### `flexhaz::dfr_dist` (parent class, full inheritance)

Constructing via `dfr_dist(rate = ..., par = ..., cum_haz_rate = ..., score_fn = ..., hess_fn = ...)`
yields hazard / survival / CDF / density / quantile / sampling / log-likelihood
/ MLE fitting for free. `dfr_dist_series()` calls this constructor and then
*adds* its own metadata (`$components`, `$layout`, `$m`, `$n_par`) and
re-classes the result.

### `likelihood.model` (inherited via `dfr_dist`)

`loglik`, `score`, `hess_loglik`, `assumptions`, and `fit` are the
relevant generics. `assumptions.dfr_dist_series` is the only one this package
implements directly; the rest dispatch into the analytical `score_fn` /
`hess_fn` closures attached at construction time (see Decomposed
score/Hessian below).

### `dist.structure` (v0.2.0, topology and reliability operations)

Implemented directly on `dfr_dist_series` (specialized for the trivial series
topology):

| Generic | Method | Result |
|---|---|---|
| `phi(x, state)` | series specialization | `as.integer(all(state == 1))` |
| `min_paths(x)` | series specialization | `list(seq_len(m))` (one path = all components) |
| `min_cuts(x)` | series specialization | `lapply(seq_len(m), function(j) j)` (each component is a singleton cut) |
| `system_signature(x)` | series specialization | `c(1, rep(0, m - 1))` |

Inherited from `dist.structure` defaults (no method registered here):
`reliability`, `dual`, `is_coherent`, `critical_states`,
`structural_importance`, `birnbaum_importance`, `criticality_importance`,
`vesely_fussell_importance`, `substitute_component`, `compose_systems`,
`system_lifetime`, `system_censoring`. These all derive from `phi` / `min_paths`,
so the four explicit specializations above are sufficient to make the entire
`dist.structure` surface work on series systems.

`ncomponents` and `component` are *re-exported* from `dist.structure` (not
locally defined as of v0.2.0). The methods registered here (`ncomponents.dfr_dist_series`,
`component.dfr_dist_series`) implement those re-exported generics. The previous
local-generic version caused a latent namespace collision when both
packages were loaded; v0.2.0 fixed that.

## Parameter layout system

Parameters across all components are stored as a single concatenated vector
`theta = (theta_1, ..., theta_m)`. The `$layout` field maps global indices
back to per-component slices. With Weibull(shape, scale) + Exponential(rate) +
Gompertz(a, b):

```r
layout    = list(1:2, 3, 4:5)
par       = c(shape1, scale1, rate2, a3, b3)
par[layout[[j]]]  -> parameters for component j
```

This is the mechanism that lets a flat-vector optimizer (`stats::optim`,
`nlm`, etc.) drive a structured multi-component model. The system-level
closures internally loop `for (j in 1:m)` and call
`components[[j]]$rate(t, par[layout[[j]]], ...)`.

## Closure-based design

All distribution functions return *closures* `function(t, par = NULL, ...)`,
not values. The system-level `sys_rate` closure (the `rate` slot passed to
`dfr_dist()`) sums component rates. If **all** components provide
`cum_haz_rate`, the system gets an analytical `sys_cum_haz`; otherwise
`cum_haz_rate` is `NULL` and the parent `dfr_dist` machinery falls back to
numerical integration via `stats::integrate`.

`component_hazard()` returns a closure that *lazily* reads `x$par` so updates
to system parameters after construction (e.g., after `fit()`) propagate to
previously-extracted hazard closures. This is intentional: the v0.1.1 NEWS
calls out a bug where eager parameter capture caused stale-default behavior.

## Decomposed score and Hessian (the non-obvious part)

The constructor attaches custom `score_fn` and `hess_fn` to the underlying
`dfr_dist` object. **Do not replace these with `numDeriv` on the full system
log-likelihood.** They exploit four structural facts:

1. **All-censored trick**: For a single component `j`, calling its analytical
   `score_fn` on data with `delta` overwritten to all-zeros gives
   `-sum_i dH_j(t_i)/dtheta_j` exactly (no contribution from exact-failure
   hazard terms). This is how `cum_haz_deriv` is extracted analytically when
   the component provides `score_fn` (e.g., for the closed-form exponential
   and Weibull `dfr_dist`s).

2. **Per-component dimensions**: When the component has no analytical
   `score_fn`, fall back to `numDeriv::grad` on a per-component
   cumulative-hazard sum (parameter dimension `n_par[j]`), not the full system
   (dimension `sum(n_par)`). For an `m`-component system, this is roughly an
   `m`-fold speedup.

3. **Hazard derivative via Jacobian**: For exact observations, the hazard
   contribution to the score is
   `sum_{i: exact} dh_j(t_i)/dtheta_j / h_sys(t_i)`. We compute the per-component
   Jacobian once via `numDeriv::jacobian` and reuse it.

4. **Hessian block structure**:
   - Diagonal blocks `(k, k)`:
     `sum_exact[d2h_k/dtheta_k^2 / h_sys - (dh_k)(dh_k)^T / h_sys^2] - sum_i d2H_k/dtheta_k^2`.
   - Off-diagonal blocks `(k, l)`, `k != l`:
     `-sum_exact (dh_k)(dh_l)^T / h_sys^2`. **These cross blocks reuse the per-component
     `dh_list` already computed for the score.** No additional numerical-Jacobian
     calls are needed for cross terms.

This decomposition is the reason `fit()` on a 5- or 10-component mixed-family
series is tractable. Reverting to a naive `numDeriv::grad(neg_loglik, par)`
would multiply runtime by roughly the parameter count.

## Source file map

| File | Contents |
|---|---|
| `R/dfr_dist_series.R` | `dfr_dist_series()` constructor (also defines `sys_rate`, `sys_cum_haz`, `sys_score_fn`, `sys_hess_fn` inline as closures), `is_dfr_dist_series()`, `print.dfr_dist_series`, `assumptions.dfr_dist_series`, and the four `dist.structure` topology specializations (`phi`, `min_paths`, `min_cuts`, `system_signature` are actually in `R/methods.R`; everything else in this file). |
| `R/methods.R` | S3 method bodies: `ncomponents`, `component`, `param_layout`, `component_hazard`, `sample_components`, and the four `dist.structure` topology specializations (`phi`, `min_paths`, `min_cuts`, `system_signature`). |
| `R/generic_functions.R` | `UseMethod` shells and roxygen for the series-specific generics (`param_layout`, `component_hazard`, `sample_components`). Note: `ncomponents` and `component` are *not* defined here, they come from `dist.structure`. |
| `R/reexports.R` | Re-exports the public surface from dependency packages so `library(serieshaz)` alone gives users `hazard`, `surv`, `cdf`, `inv_cdf`, `sampler`, `params`, `loglik`, `score`, `hess_loglik`, `assumptions`, `fit`, `cum_haz`, `dfr_dist`, `is_dfr_dist`, `dfr_exponential`, `dfr_weibull`, `dfr_gompertz`, `dfr_loglogistic`, `ncomponents`, `component`, `phi`, `min_paths`, `min_cuts`, `system_signature`. |

## Test layout

Test fixtures live in `tests/testthat/helper-fixtures.R` (`make_exp_series`,
`make_weibull_series`, `make_mixed_series`, and friends for censored
data frames). Each test file targets a specific concern:

| File | Verifies |
|---|---|
| `test-constructor.R` | Component-list validation, layout construction, `par` resolution, `n_par` inference. |
| `test-composition.R` | Nested series, mixed-family series, the hazard-sum identity, `S_sys = prod S_j`. |
| `test-methods.R` | Introspection generics and their input validation (`j` must be a single integer in range). |
| `test-sampling.R` | `sample_components` and the `min_j T_j = T_sys` identity via Monte Carlo. |
| `test-fit.R` | MLE fitting, including the identifiability workaround for exponential-only series. |
| `test-derivatives.R` | Analytical score/Hessian agree with `numDeriv` on the full log-likelihood. |
| `test-score-decomposed.R` / `test-hessian-decomposed.R` | The all-censored trick produces the correct cumulative-hazard derivative; cross-component Hessian blocks match the outer-product formula. |
| `test-edge-cases.R` | Left-censored data, `t = 0` boundaries, `j` validation paths. |
| `test-dist-structure-interop.R` | 45 tests covering every `dist.structure` generic against a series system: `phi`, `min_paths`, `min_cuts`, `system_signature`, `reliability`, `dual`, `is_coherent`, `critical_states`, `structural_importance`, `birnbaum_importance`, `criticality_importance`, `vesely_fussell_importance`, `substitute_component`, `compose_systems`, `system_lifetime`, `system_censoring`. |

Gold-standard verification strategies recurring across files:

- **Sum-of-rates collapse**: 3 exponential components with rates
  `c(r1, r2, r3)` must behave identically to a single `dfr_exponential(r1 + r2 + r3)`.
- **Product-of-survivals**: `S_sys(t)` must equal `prod_j S_j(t)` for any
  family mix.
- **Minimum of components**: empirical `min(sample_components(sys, n))` must
  match `sampler(sys)(n)` in distribution.

## Identifiability constraint

Exponential series systems are **not identifiable** from system-level data
alone. Only `sum(rates)` is identifiable. Fitting tests check
`sum(coef(result))` rather than individual parameters when the model is pure
exponential. Mixed-family series (e.g., Weibull + Gompertz, or Weibull +
Exponential) *are* identifiable because the components have qualitatively
different hazard shapes that pin down each parameter through temporal signal.

This constraint shows up only at the *statistical* level. All distributional
operations (`hazard`, `surv`, `sampler`, ...) are well-defined for any
combination of component parameters.

## Conventions

- **No em-dashes** anywhere (project-wide convention, enforced by a PostToolUse
  hook). Use commas, colons, periods, or parentheses instead.
- **Survival data columns**: `t` for observation time, `delta` for censoring
  (1 exact, 0 right-censored, -1 left-censored).
- **Component prototype pattern**: throughout the rlang ecosystem
  (`flexhaz`, `kofn`, `serieshaz`, `dist.structure`), components are specified
  by passing a `dfr_dist` *object* (e.g., `dfr_exponential()`,
  `dfr_weibull(shape = 2, scale = 100)`), not a family string. The S3 class
  drives dispatch; the parameter values on the prototype are ignored when
  the parameters will be estimated by `fit()`.
