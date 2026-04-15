---
title: 'serieshaz: Series System Distributions from Composed Hazard Functions in R'
tags:
  - R
  - reliability engineering
  - survival analysis
  - series system
  - hazard function
  - maximum likelihood estimation
authors:
  - name: Alexander Towell
    orcid: 0000-0001-6443-9897
    affiliation: 1
affiliations:
  - name: Southern Illinois University Edwardsville
    index: 1
date: "29 March 2026"
bibliography: paper.bib
---

# Summary

`serieshaz` is an R package for composing multiple dynamic failure rate
distributions into series system distributions. A series system fails when
*any* component fails, so the system hazard is the sum of component hazards:
$h_{\text{sys}}(t) = \sum_{j=1}^{m} h_j(t, \theta_j)$. The package accepts
any combination of `dfr_dist` objects from the `flexhaz` package
[@flexhaz]---Weibull, exponential, Gompertz, log-logistic, or user-defined
hazard functions---and produces a single distribution object with the full
inferential machinery: hazard, survival, CDF, density, quantile, sampling,
log-likelihood, score, Hessian, and maximum likelihood estimation via `fit()`.

```r
library(serieshaz)
server <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 500),     # disk wear-out
    dfr_exponential(0.001),                   # random memory failure
    dfr_gompertz(a = 0.0001, b = 0.02)       # PSU degradation
))
S <- surv(server)
S(100)   # system survival probability at t = 100
```

The resulting `dfr_dist_series` object inherits from `dfr_dist`, so all
methods defined at any level of the class hierarchy (`dfr_dist`,
`likelihood_model`, `univariate_dist`, `dist`) work automatically without
reimplementation. Series-specific introspection functions---`ncomponents()`,
`component()`, `component_hazard()`, `param_layout()`, and
`sample_components()`---enable hazard decomposition and failure attribution
analysis.

# Statement of Need

Reliability engineers routinely analyze multi-component systems where each
component exhibits a distinct failure mechanism: mechanical wear-out (Weibull),
random shocks (exponential), or accelerating degradation (Gompertz). When
these components are arranged in series---the system fails at the first
component failure---the analyst needs a system-level distribution that
combines heterogeneous hazard functions and supports maximum likelihood
estimation from system-level failure data. Researchers working with series
system models in reliability analysis and survival analysis need a tool that
composes arbitrary parametric hazard functions into a system distribution,
provides the full distributional interface, and delivers efficient MLE with
decomposed derivatives that exploit the additive structure of the series
hazard.

# State of the Field

Several R packages address survival modeling and system reliability, but none
compose arbitrary user-defined hazard functions into series system
distributions with decomposed inference.

The `survival` package [@survival] provides Kaplan-Meier estimation, Cox
proportional hazards, and parametric accelerated failure time models, but does
not model multi-component system structure or compose component hazards.
`flexsurv` [@flexsurv] supports flexible parametric survival models with
user-defined distributions, but requires the user to specify the density or
cumulative hazard alongside the hazard, and targets regression modeling
rather than system composition. The `ReliabilityTheory` package
[@ReliabilityTheory] computes system signatures and survival signatures for
structural reliability analysis, including series and parallel configurations,
but operates on discrete system topology rather than composing parametric
hazard functions for MLE. The `FaultTree` package [@FaultTree] constructs
fault trees for probabilistic risk assessment but addresses system-level
failure probability calculation, not parametric distribution composition.
The Python `reliability` library [@reliability_py] provides a
`Competing_Risks_Model` that multiplies survival functions from built-in
distributions, but restricts components to its fixed catalog and does not
provide decomposed score and Hessian computation. The `fitdistrplus` package
[@fitdistrplus] fits univariate parametric distributions to censored data but
has no notion of multi-component system structure.

`serieshaz` fills this gap by accepting any `dfr_dist` component---including
user-defined hazard functions with no closed-form cumulative hazard---and
producing a series system distribution with efficient decomposed derivatives.

# Software Design

The `dfr_dist_series` constructor accepts a list of `dfr_dist` objects and
builds a composite `dfr_dist` whose hazard rate sums the component hazards.
Parameters across all $m$ components are stored as a single concatenated
vector $\boldsymbol{\theta} = (\theta_1, \ldots, \theta_m)$, with a layout
mapping global indices to per-component slices. This flat-vector design
enables standard optimizers to work directly on the system parameter space.

The score and Hessian exploit the series sum structure rather than applying
`numDeriv` [@numDeriv] to the full system log-likelihood. Because
$H_{\text{sys}} = \sum_j H_j$ and $h_{\text{sys}} = \sum_j h_j$, the
gradient with respect to component $k$'s parameters decomposes into a hazard
derivative term and a cumulative hazard derivative term, each involving only
component $k$'s lower-dimensional parameter space. When a component provides
an analytical score function, the cumulative hazard derivative is extracted
analytically via the "all-censored trick": evaluating the component's score
with all observations marked right-censored isolates the
$-\partial H_k / \partial \theta_k$ term. The Hessian has a natural block
structure---cross-component blocks
$\partial^2 \ell / \partial \theta_k \partial \theta_l^\top$ ($k \neq l$)
are free, computed entirely from the rate Jacobians already obtained for the
score. This decomposition yields roughly an $m$-fold speedup over naive
numerical differentiation of the full system log-likelihood.

When all components provide analytical cumulative hazard functions, the series
system computes $H_{\text{sys}}(t) = \sum_j H_j(t)$ exactly; otherwise it
falls back to numerical integration. Series systems can be nested---a
`dfr_dist_series` is itself a valid `dfr_dist` component---enabling
hierarchical system modeling.

# Research Impact Statement

`serieshaz` is part of an ecosystem of R packages for reliability analysis
with masked failure data. It serves as the compositional layer between
`flexhaz`, which defines individual component distributions via arbitrary
hazard functions, and `maskedhaz`, which constructs likelihood models for
masked series system data where the failing component is not observed. A
companion package, `maskedcauses` [@maskedcauses], provides closed-form MLE
for exponential and Weibull masked series systems; `maskedhaz` cross-validates
its numerical results against these analytical solutions. The ecosystem is
available on r-universe at <https://queelius.r-universe.dev>.

# AI Usage Disclosure

Claude Code (Anthropic) was used to assist with code refactoring, test
generation, and drafting of documentation including this manuscript. All
generated content was reviewed, edited, and validated by the author. The
package design, mathematical formulations, and research decisions are
entirely the author's own work.

# Acknowledgements

The author thanks the R Core Team for the R statistical computing environment
[@R] and the developers of `numDeriv` for numerical derivative infrastructure.

# References
