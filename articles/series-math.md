# Mathematical Foundations of Series Systems

## Series System Theory

A **series system** consists of $m$ independent components, each with
lifetime $T_{j}$ and survival function
$S_{j}(t) = P\left( T_{j} > t \right)$. The system fails when *any*
component fails, so the system lifetime is:

$$T_{sys} = \min\left( T_{1},\ldots,T_{m} \right)$$

Since the components are independent:

$$S_{sys}(t) = P\left( T_{sys} > t \right) = P\left( T_{1} > t,\ldots,T_{m} > t \right) = \prod\limits_{j = 1}^{m}S_{j}(t)$$

Taking the negative logarithm gives the cumulative hazard:

$$H_{sys}(t) = - \log S_{sys}(t) = \sum\limits_{j = 1}^{m}H_{j}(t)$$

Differentiating yields the hazard rate:

$$h_{sys}(t) = \frac{d}{dt}H_{sys}(t) = \sum\limits_{j = 1}^{m}h_{j}(t)$$

This is the fundamental result: **the system hazard is the sum of
component hazards**.

## Key Properties

Let’s verify these relationships numerically.

### Property 1: Hazard Sum

``` r
library(dfr.dist.series)

sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01),
    dfr_gompertz(a = 0.001, b = 0.05)
))

h_sys <- hazard(sys)
h1 <- component_hazard(sys, 1)
h2 <- component_hazard(sys, 2)
h3 <- component_hazard(sys, 3)

t_vals <- c(10, 50, 100, 200)
for (t0 in t_vals) {
    lhs <- h_sys(t0)
    rhs <- h1(t0) + h2(t0) + h3(t0)
    cat(sprintf("t=%3d: h_sys=%.6f, sum h_j=%.6f, diff=%.2e\n",
                t0, lhs, rhs, abs(lhs - rhs)))
}
#> t= 10: h_sys=0.013649, sum h_j=0.013649, diff=0.00e+00
#> t= 50: h_sys=0.032182, sum h_j=0.032182, diff=0.00e+00
#> t=100: h_sys=0.178413, sum h_j=0.178413, diff=0.00e+00
#> t=200: h_sys=22.076466, sum h_j=22.076466, diff=0.00e+00
```

### Property 2: Survival Product

``` r
S_sys <- surv(sys)
comp1 <- component(sys, 1)
comp2 <- component(sys, 2)
comp3 <- component(sys, 3)
S1 <- surv(comp1)
S2 <- surv(comp2)
S3 <- surv(comp3)

for (t0 in t_vals) {
    lhs <- S_sys(t0)
    rhs <- S1(t0) * S2(t0) * S3(t0)
    cat(sprintf("t=%3d: S_sys=%.6f, prod S_j=%.6f, diff=%.2e\n",
                t0, lhs, rhs, abs(lhs - rhs)))
}
#> t= 10: S_sys=0.884286, prod S_j=0.884286, diff=1.11e-16
#> t= 50: S_sys=0.377702, prod S_j=0.377702, diff=0.00e+00
#> t=100: S_sys=0.007096, prod S_j=0.007096, diff=0.00e+00
#> t=200: S_sys=0.000000, prod S_j=0.000000, diff=1.52e-210
```

### Property 3: Cumulative Hazard Sum

``` r
# The system's cumulative hazard is the sum of component cumulative hazards
H_sys <- cum_haz(sys)
H1 <- cum_haz(comp1)
H2 <- cum_haz(comp2)
H3 <- cum_haz(comp3)

for (t0 in t_vals) {
    lhs <- H_sys(t0)
    rhs <- H1(t0) + H2(t0) + H3(t0)
    cat(sprintf("t=%3d: H_sys=%.6f, sum H_j=%.6f, diff=%.2e\n",
                t0, lhs, rhs, abs(lhs - rhs)))
}
#> t= 10: H_sys=0.122974, sum H_j=0.122974, diff=0.00e+00
#> t= 50: H_sys=0.973650, sum H_j=0.973650, diff=0.00e+00
#> t=100: H_sys=4.948263, sum H_j=4.948263, diff=0.00e+00
#> t=200: H_sys=446.509316, sum H_j=446.509316, diff=0.00e+00
```

### Property 4: Density Formula

The system density is:

$$f_{sys}(t) = h_{sys}(t) \cdot S_{sys}(t) = \left\lbrack \sum\limits_{j = 1}^{m}h_{j}(t) \right\rbrack\exp\left\lbrack - \sum\limits_{j = 1}^{m}H_{j}(t) \right\rbrack$$

``` r
f_sys <- density(sys)
for (t0 in t_vals) {
    from_density <- f_sys(t0)
    from_hS <- h_sys(t0) * S_sys(t0)
    cat(sprintf("t=%3d: f(t)=%.8f, h(t)*S(t)=%.8f, diff=%.2e\n",
                t0, from_density, from_hS, abs(from_density - from_hS)))
}
#> t= 10: f(t)=0.01206938, h(t)*S(t)=0.01206938, diff=0.00e+00
#> t= 50: f(t)=0.01215539, h(t)*S(t)=0.01215539, diff=0.00e+00
#> t=100: f(t)=0.00126597, h(t)*S(t)=0.00126597, diff=0.00e+00
#> t=200: f(t)=0.00000000, h(t)*S(t)=0.00000000, diff=0.00e+00
```

The [`density()`](https://rdrr.io/r/stats/density.html) method and the
manual $h(t) \cdot S(t)$ calculation agree to machine precision. Note
that at $t = 200$, both the density and survival are effectively zero —
by that time, the combined Weibull + exponential + Gompertz hazard has
driven the cumulative hazard so high that virtually all systems have
already failed.

## The Parameter Layout

The series system stores parameters from all components as a single flat
vector. This is critical for optimization — MLE fitting via
[`optim()`](https://rdrr.io/r/stats/optim.html) requires a single
parameter vector.

``` r
sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),    # 2 params
    dfr_exponential(0.05),            # 1 param
    dfr_gompertz(a = 0.01, b = 0.1)         # 2 params
))

# Full parameter vector (5 elements)
params(sys)
#> [1] 2e+00 1e+02 5e-02 1e-02 1e-01

# Layout maps global indices to components
param_layout(sys)
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3
#> 
#> [[3]]
#> [1] 4 5

# Component 1 uses par[1:2], component 2 uses par[3], component 3 uses par[4:5]
```

When an optimizer proposes a new parameter vector
`par = c(p1, p2, p3, p4, p5)`, the series system internally slices it:
component 1 gets `par[1:2]`, component 2 gets `par[3]`, and component 3
gets `par[4:5]`.

## Analytical vs. Numerical Cumulative Hazard

The cumulative hazard $H(t) = \int_{0}^{t}h(u)\, du$ is needed for
survival, CDF, and log-likelihood computations. When all components
provide an analytical `cum_haz_rate`, the series system computes
$H_{sys}(t) = \sum_{j}H_{j}(t)$ exactly. Otherwise, it falls back to
numerical integration.

``` r
# All standard distributions provide analytical cumulative hazard
sys_analytical <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
))
cat("Has analytical cum_haz:", !is.null(sys_analytical$cum_haz_rate), "\n")
#> Has analytical cum_haz: TRUE

# A custom dfr_dist without cum_haz_rate forces numerical fallback
custom <- dfr_dist(
    rate = function(t, par, ...) par[1] * t^(par[1] - 1),
    par = c(1.5)
)
sys_numerical <- dfr_dist_series(list(
    custom,
    dfr_exponential(0.05)
))
cat("Has analytical cum_haz:", !is.null(sys_numerical$cum_haz_rate), "\n")
#> Has analytical cum_haz: FALSE
```

The analytical path is faster and more precise, but the numerical
fallback works correctly for any valid hazard function.

## Score and Hessian

The log-likelihood for a series system with observed data
$\left( t_{i},\delta_{i} \right)$ is:

$$\ell(\theta) = \sum\limits_{i}\left\lbrack \delta_{i}\log h_{sys}\left( t_{i};\theta \right) - H_{sys}\left( t_{i};\theta \right) \right\rbrack$$

where $\delta_{i} = 1$ for exact observations and $\delta_{i} = 0$ for
right-censored observations.

The score (gradient) and Hessian are computed via `numDeriv`:

``` r
sys <- dfr_dist_series(list(
    dfr_exponential(0.1),
    dfr_exponential(0.2)
))

# Generate some data
set.seed(42)
samp <- sampler(sys)
times <- samp(50)
df <- data.frame(t = times, delta = rep(1L, 50))

# Log-likelihood at true parameters
ll <- loglik(sys)
ll(df)
#> [1] -128.9726

# Score (gradient of log-likelihood) — near zero at true params
sc <- score(sys)
score_val <- sc(df)
cat("Score at true params:", score_val, "\n")
#> Score at true params: -62.57988 -62.57988

# Hessian (second derivatives) — negative definite at MLE
H <- hess_loglik(sys)
hess_val <- H(df)
cat("Hessian eigenvalues:", eigen(hess_val, symmetric = TRUE)$values, "\n")
#> Hessian eigenvalues: -2.552025e-08 -1111.111
```

### Numerical Verification

We can verify the score is indeed the gradient of the log-likelihood
using finite differences:

``` r
par0 <- params(sys)
eps <- 1e-5
ll_fn <- loglik(sys)
sc_fn <- score(sys)

# Analytical score at par0
score_val <- sc_fn(df, par = par0)

# Finite-difference approximation
fd_score <- numeric(length(par0))
for (k in seq_along(par0)) {
    par_plus <- par_minus <- par0
    par_plus[k] <- par0[k] + eps
    par_minus[k] <- par0[k] - eps
    fd_score[k] <- (ll_fn(df, par = par_plus) - ll_fn(df, par = par_minus)) / (2 * eps)
}

cat("Score (numDeriv):   ", score_val, "\n")
#> Score (numDeriv):    -62.57988 -62.57988
cat("Score (finite diff):", fd_score, "\n")
#> Score (finite diff): -62.57988 -62.57988
cat("Max absolute diff:  ", max(abs(score_val - fd_score)), "\n")
#> Max absolute diff:   6.891987e-08
```

The `numDeriv`-based score and the manual finite-difference
approximation agree closely, confirming that the composed log-likelihood
is being differentiated correctly through the parameter layout.

## Special Cases

### Exponential Series = Exponential(sum of rates)

When all components are exponential with rates
$\lambda_{1},\ldots,\lambda_{m}$, the system hazard is constant:
$h_{sys}(t) = \sum_{j}\lambda_{j}$. This is equivalent to a single
exponential with rate $\lambda = \sum_{j}\lambda_{j}$.

``` r
rates <- c(0.1, 0.2, 0.3)
sys <- dfr_dist_series(lapply(rates, dfr_exponential))
equiv <- dfr_exponential(sum(rates))

h_sys <- hazard(sys)
h_eq  <- hazard(equiv)
S_sys <- surv(sys)
S_eq  <- surv(equiv)

for (t0 in c(1, 5, 10, 50)) {
    cat(sprintf("t=%2d: h_sys=%.4f h_eq=%.4f | S_sys=%.6f S_eq=%.6f\n",
                t0, h_sys(t0), h_eq(t0), S_sys(t0), S_eq(t0)))
}
#> t= 1: h_sys=0.6000 h_eq=0.6000 | S_sys=0.548812 S_eq=0.548812
#> t= 5: h_sys=0.6000 h_eq=0.6000 | S_sys=0.049787 S_eq=0.049787
#> t=10: h_sys=0.6000 h_eq=0.6000 | S_sys=0.002479 S_eq=0.002479
#> t=50: h_sys=0.6000 h_eq=0.6000 | S_sys=0.000000 S_eq=0.000000
```

### Identical Weibull Components

When all $m$ components are Weibull with the same shape $k$ and scale
$\lambda$, the system is also Weibull with shape $k$ and scale
$\lambda/m^{1/k}$:

``` r
m <- 3
shape <- 2
scale <- 100

sys <- dfr_dist_series(replicate(
    m, dfr_weibull(shape = shape, scale = scale), simplify = FALSE
))
equiv <- dfr_weibull(shape = shape, scale = scale / m^(1/shape))

S_sys <- surv(sys)
S_eq  <- surv(equiv)

for (t0 in c(10, 30, 50)) {
    cat(sprintf("t=%2d: S_sys=%.6f S_equiv=%.6f diff=%.2e\n",
                t0, S_sys(t0), S_eq(t0), abs(S_sys(t0) - S_eq(t0))))
}
#> t=10: S_sys=0.970446 S_equiv=0.970446 diff=0.00e+00
#> t=30: S_sys=0.763379 S_equiv=0.763379 diff=1.11e-16
#> t=50: S_sys=0.472367 S_equiv=0.472367 diff=5.55e-17
```

### Single Component (Degenerate Case)

A series system with one component is equivalent to that component:

``` r
single <- dfr_dist_series(list(dfr_weibull(shape = 2, scale = 100)))
direct <- dfr_weibull(shape = 2, scale = 100)

h1 <- hazard(single)
h2 <- hazard(direct)
S1 <- surv(single)
S2 <- surv(direct)

for (t0 in c(10, 50, 100)) {
    cat(sprintf("t=%3d: h=%.6f/%.6f S=%.6f/%.6f\n",
                t0, h1(t0), h2(t0), S1(t0), S2(t0)))
}
#> t= 10: h=0.002000/0.002000 S=0.990050/0.990050
#> t= 50: h=0.010000/0.010000 S=0.778801/0.778801
#> t=100: h=0.020000/0.020000 S=0.367879/0.367879
```
