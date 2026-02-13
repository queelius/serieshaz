# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
# Install dependencies (from GitHub remotes)
Rscript -e 'remotes::install_deps(dependencies = TRUE)'

# Build and check
R CMD build .
R CMD check dfr.dist.series_*.tar.gz

# Run all tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-constructor.R")'

# Run tests matching a pattern
Rscript -e 'testthat::test_dir("tests/testthat", filter = "composition")'

# Regenerate NAMESPACE and man/ pages (after roxygen changes)
Rscript -e 'roxygen2::roxygenise()'

# Test coverage
Rscript -e 'covr::package_coverage()'
```

## Architecture

This package composes multiple `dfr_dist` objects (from the `dfr.dist` package) into a **series system distribution**. A series system fails when *any* component fails, so the system hazard is the sum of component hazards: `h_sys(t) = sum_j h_j(t)`.

### Class Hierarchy

`dfr_dist_series` inherits from `dfr_dist` which inherits from `likelihood_model` -> `univariate_dist` -> `dist`. This inheritance chain means all existing methods (hazard, survival, CDF, density, sampling, log-likelihood, MLE fitting) work automatically on series systems without reimplementation.

### Three Dependency Packages (all from `github::queelius/`)

- **`algebraic.dist`**: Base distribution generics (`hazard`, `surv`, `cdf`, `inv_cdf`, `sampler`, `params`)
- **`likelihood.model`**: Statistical inference generics (`loglik`, `score`, `hess_loglik`, `assumptions`, `fit`)
- **`dfr.dist`**: Dynamic failure rate distributions (`dfr_dist` constructor, named distributions like `dfr_exponential`, `dfr_weibull`, `dfr_gompertz`, `dfr_loglogistic`)

### Parameter Layout System

Parameters across all components are stored as a single concatenated vector. The `$layout` field maps global indices to per-component indices. For example, with Weibull(shape, scale) + Exponential(rate): `layout = list(1:2, 3)`, so the global parameter vector `c(shape1, scale1, rate2)` gets sliced as `par[1:2]` for component 1 and `par[3]` for component 2. This design enables standard optimizers to work on a flat parameter vector while the series system distributes parameters to the correct components.

### Closure-Based Design

All distribution functions (`hazard`, `surv`, `cdf`, etc.) return **closures** — `function(t, par = NULL, ...)` — not values. The system-level closures internally loop over components, slice the parameter vector via `layout`, and sum hazards. If all components provide `cum_haz_rate`, the series system gets an analytical cumulative hazard; otherwise it falls back to numerical integration. Score and Hessian always use `numDeriv` fallback.

### Source File Organization

- `R/dfr_dist_series.R`: Constructor (`dfr_dist_series`), type predicate, print method, assumptions
- `R/methods.R`: S3 method implementations (`ncomponents`, `component`, `param_layout`, `component_hazard`, `sample_components`)
- `R/generic_functions.R`: Generic function definitions for series-specific generics
- `R/reexports.R`: Re-exports from dependency packages

### Testing Patterns

Test fixtures are in `tests/testthat/helper-fixtures.R`. The gold-standard verification strategy is: 3 exponential components with rates (r1, r2, r3) must behave identically to a single exponential with rate (r1 + r2 + r3). Tests also verify S_sys(t) = product of S_j(t) for Weibull series, and that min(component samples) follows the system distribution.

### Key Identifiability Constraint

Exponential series systems are **not identifiable** from system-level data alone — only the sum of rates is identifiable. The fitting tests account for this by checking `sum(coef(result))` rather than individual parameters.
