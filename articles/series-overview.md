# Series System Distributions: Overview

## The Problem

Consider a server with three independent failure modes:

- **Disk** — mechanical wear follows a Weibull distribution
- **Memory** — random bit-flip failures follow an exponential
  distribution
- **Power supply** — aging degradation follows a Gompertz distribution

The server fails when *any* component fails. This is a **series system**
— like a chain that breaks at its weakest link. The `serieshaz` package
lets you compose these failure modes into a single distribution that
captures the full system behavior.

## What Is a Series System?

Given $m$ independent components with hazard functions
$h_{1}(t),\ldots,h_{m}(t)$, the series system hazard is simply the
**sum**:

$$h_{sys}(t) = \sum\limits_{j = 1}^{m}h_{j}\left( t,\theta_{j} \right)$$

Equivalently, the system survival is the **product** of component
survivals:

$$S_{sys}(t) = \prod\limits_{j = 1}^{m}S_{j}\left( t,\theta_{j} \right)$$

This follows directly from component independence: the probability that
the system survives to time $t$ is the probability that *all* components
survive to time $t$.

## Quick Start

``` r
library(serieshaz)

# Define three failure mode components
disk    <- dfr_weibull(shape = 2, scale = 500)
memory  <- dfr_exponential(0.001)
psu     <- dfr_gompertz(a = 0.0001, b = 0.02)

# Compose into a series system
server <- dfr_dist_series(list(disk, memory, psu))
print(server)
#> Series system distribution with 3 components
#>   Component 1: 2 param(s) [  2, 500]
#>   Component 2: 1 param(s) [0.001]
#>   Component 3: 2 param(s) [1e-04, 2e-02]
#> System hazard: h_sys(t) = sum_j h_j(t, theta_j)
#> Survival: S_sys(t) = exp(-H_sys(t)) = prod_j S_j(t, theta_j)
```

The system object has 5 parameters total: 2 from Weibull + 1 from
exponential + 2 from Gompertz. You can see how they map to components:

``` r
param_layout(server)
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3
#> 
#> [[3]]
#> [1] 4 5
```

## Everything from flexhaz Works

Because `dfr_dist_series` inherits from `dfr_dist`, all standard
distribution methods work out of the box — no special code needed.

### Hazard and Survival

``` r
h <- hazard(server)
S <- surv(server)

# Evaluate at t = 100 hours
cat(sprintf("System hazard at t=100:  %.6f\n", h(100)))
#> System hazard at t=100:  0.002539
cat(sprintf("System survival at t=100: %.4f\n", S(100)))
#> System survival at t=100: 0.8420
```

The hazard is low ($\approx 0.0025$) because the Weibull scale is 500
and the exponential rate is only 0.001. The survival probability of
$\approx 84\%$ at $t = 100$ tells us most servers are still running at
100 hours.

### CDF and Density

``` r
F_sys <- cdf(server)
f_sys <- density(server)

cat(sprintf("P(T <= 100) = %.4f\n", F_sys(100)))
#> P(T <= 100) = 0.1580
cat(sprintf("f(100)      = %.6f\n", f_sys(100)))
#> f(100)      = 0.002138
```

The CDF and survival sum to 1: $F(100) + S(100) = 1$. The density
$f(t) = h(t) \cdot S(t)$ gives the instantaneous failure rate weighted
by the probability of surviving to that point.

### Sampling

``` r
samp <- sampler(server)
set.seed(42)
times <- samp(5)
times
#> [1] 294.9884 302.0998 150.3888 274.2224 236.8369
```

### Log-Likelihood

``` r
ll <- loglik(server)
df <- data.frame(t = c(100, 200, 150, 300, 250), delta = c(1, 1, 0, 1, 0))
ll(df)
#> [1] -18.97089
```

### Fitting (MLE)

Fitting a 5-parameter model (Weibull + exponential + Gompertz) from
system-level data alone is challenging — the components’ hazard
contributions overlap, making individual parameters hard to identify.
For a cleaner demonstration, we fit a simpler 2-component model:

``` r
# Simpler 2-component model for a clean fit demo
fit_model <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01)
))
true_par <- c(2, 100, 0.01)

set.seed(42)
fit_samp <- sampler(fit_model)
train <- data.frame(t = fit_samp(500), delta = rep(1L, 500))

solver <- fit(fit_model)
result <- suppressWarnings(solver(train, par = c(1.5, 80, 0.02)))

cat("True parameters:  ", true_par, "\n")
#> True parameters:   2 100 0.01
cat("Fitted parameters:", round(coef(result), 3), "\n")
#> Fitted parameters: 2.494 112.753 0.013
```

With 500 observations and a mixed-type model (Weibull + exponential),
the MLE recovers the parameters reasonably well. The Weibull’s
increasing hazard shape is distinguishable from the exponential’s
constant hazard, making this system identifiable. See
[`vignette("series-fitting")`](https://queelius.github.io/serieshaz/articles/series-fitting.md)
for a thorough treatment of identifiability, censoring, and model
selection.

## System Introspection

The package provides functions specifically for understanding series
systems.

### Component Count and Extraction

``` r
ncomponents(server)  # 3
#> [1] 3

# Extract the Weibull (disk) component
disk_comp <- component(server, 1)
params(disk_comp)
#> [1]   2 500
```

### Component Hazard Decomposition

``` r
# Get per-component hazard closures
h1 <- component_hazard(server, 1)  # disk
h2 <- component_hazard(server, 2)  # memory
h3 <- component_hazard(server, 3)  # PSU

# At t = 200, which component contributes most?
t0 <- 200
hazards <- c(Disk = h1(t0), Memory = h2(t0), PSU = h3(t0))
cat(sprintf("Disk:   %.6f (%.1f%%)\n", hazards[1], 100 * hazards[1] / sum(hazards)))
#> Disk:   0.001600 (19.9%)
cat(sprintf("Memory: %.6f (%.1f%%)\n", hazards[2], 100 * hazards[2] / sum(hazards)))
#> Memory: 0.001000 (12.4%)
cat(sprintf("PSU:    %.6f (%.1f%%)\n", hazards[3], 100 * hazards[3] / sum(hazards)))
#> PSU:    0.005460 (67.7%)
cat(sprintf("System: %.6f\n", h(t0)))
#> System: 0.008060
cat(sprintf("Sum:    %.6f (matches system)\n", sum(hazards)))
#> Sum:    0.008060 (matches system)
```

At $t = 200$, the PSU’s Gompertz degradation dominates system risk
(~68%), with the disk’s Weibull wear-out contributing ~20% and the
memory’s constant exponential hazard only ~12%. This decomposition is
the key advantage of series system modeling: it reveals *which* failure
mode drives system risk at any given age. Early on, the constant memory
hazard is the largest contributor, but the accelerating Gompertz term
overtakes it by around $t = 150$.

### Failure Attribution via Sampling

``` r
set.seed(42)
mat <- sample_components(server, n = 10000)

# System lifetime = min across components
t_sys <- apply(mat, 1, min)

# Which component caused each failure?
failing <- apply(mat, 1, which.min)
cat("Failure proportions:\n")
#> Failure proportions:
round(table(failing) / length(failing), 3)
#> failing
#>     1     2     3 
#> 0.182 0.200 0.618
```

## The Ecosystem

`serieshaz` builds on three packages:

- **[algebraic.dist](https://github.com/queelius/algebraic.dist)** —
  Base distribution generics: `hazard`, `surv`, `cdf`, `inv_cdf`,
  `sampler`, `params`
- **[likelihood.model](https://github.com/queelius/likelihood.model)** —
  Statistical inference: `loglik`, `score`, `hess_loglik`, `fit`,
  `assumptions`
- **[flexhaz](https://github.com/queelius/flexhaz)** — DFR
  distributions: `dfr_dist`, `dfr_exponential`, `dfr_weibull`,
  `dfr_gompertz`, `dfr_loglogistic`

The inheritance chain is: `dfr_dist_series` → `dfr_dist` →
`likelihood_model` → `univariate_dist` → `dist`. Every method defined at
any level in this chain works automatically on series systems.

## Where to Go Next

- [`vignette("series-math")`](https://queelius.github.io/serieshaz/articles/series-math.md)
  — Mathematical foundations and derivations
- [`vignette("series-fitting")`](https://queelius.github.io/serieshaz/articles/series-fitting.md)
  — MLE fitting workflows, identifiability, and diagnostics
- [`vignette("series-advanced")`](https://queelius.github.io/serieshaz/articles/series-advanced.md)
  — Nested systems, custom components, and failure attribution
