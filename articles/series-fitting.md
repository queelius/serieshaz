# Fitting Series Systems to Data

## Introduction

In practice, we observe system-level failure data: the time at which a
system failed (or was censored), but often *not* which component caused
the failure. The `dfr.dist.series` package fits series system models to
such data via maximum likelihood estimation (MLE).

## Data Format

The fitting functions expect a data frame with two columns:

| Column  | Type    | Meaning             |
|---------|---------|---------------------|
| `t`     | numeric | Observed time       |
| `delta` | integer | Censoring indicator |

The censoring indicator values are:

- `1` — exact observation (the system failed at time `t`)
- `0` — right-censored (the system was still running at time `t`)
- `-1` — left-censored (the system had already failed by time `t`)

## Basic Fitting Workflow

### Step 1: Define the System Model

``` r
library(dfr.dist.series)

# Hypothesize: system has two failure modes
#   - Wear-out (Weibull)
#   - Random shock (exponential)
model <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01)
))
```

### Step 2: Simulate Training Data

In practice you’d have real field data. Here we simulate from known
parameters to demonstrate the workflow.

``` r
set.seed(42)
true_par <- c(2, 100, 0.01)  # shape, scale, rate

samp <- sampler(model)
n <- 300
times <- samp(n, par = true_par)
train <- data.frame(t = times, delta = rep(1L, n))

summary(train$t)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#>   0.02386  21.45824  46.56538  53.75536  77.11423 193.31375
```

### Step 3: Fit the Model

``` r
solver <- fit(model)

# Provide initial parameter guesses
init_par <- c(1.5, 80, 0.02)
result <- solver(train, par = init_par)
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
```

### Step 4: Extract Results

``` r
# Fitted parameters
cat("True parameters:", true_par, "\n")
#> True parameters: 2 100 0.01
cat("Fitted parameters:", coef(result), "\n")
#> Fitted parameters: 2.47432 109.1372 0.01211101

# Log-likelihood at MLE
cat("Log-likelihood:", logLik(result), "\n")
#> Log-likelihood: -1475.163

# Variance-covariance matrix
vcov(result)
#>              [,1]         [,2]         [,3]
#> [1,] 0.2503925591   4.37634657 7.702069e-04
#> [2,] 4.3763465699 139.44446132 1.998348e-02
#> [3,] 0.0007702069   0.01998348 4.079105e-06
```

## Initial Parameter Selection

The choice of initial parameters can affect convergence. Strategies
include:

1.  **Domain knowledge**: Use engineering estimates or handbook values
2.  **Simpler model first**: Fit a single exponential to get a rough
    rate, then use that as a starting point
3.  **Multiple starts**: Try several initial values and keep the best
    fit

``` r
# Strategy: try a few initial guesses, keep the best
init_guesses <- list(
    c(1.5, 80,  0.02),
    c(2.5, 120, 0.005),
    c(1.0, 50,  0.05)
)

best_ll <- -Inf
best_result <- NULL
for (init in init_guesses) {
    res <- tryCatch(solver(train, par = init), error = function(e) NULL)
    if (!is.null(res) && logLik(res) > best_ll) {
        best_ll <- logLik(res)
        best_result <- res
    }
}
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
cat("Best log-likelihood:", best_ll, "\n")
#> Best log-likelihood: -1475.163
cat("Best parameters:", coef(best_result), "\n")
#> Best parameters: 2.484501 109.2152 0.01213407
```

## Identifiability

### Exponential Series: Non-Identifiable

When all components are exponential, only the **sum of rates** is
identifiable from system-level data. The system hazard
$h_{sys}(t) = \sum_{j}\lambda_{j}$ is constant, so any partition of
rates that sums to the same total produces the same likelihood.

``` r
# True rates: 0.1, 0.2, 0.3 (sum = 0.6)
true_rates <- c(0.1, 0.2, 0.3)
exp_sys <- dfr_dist_series(lapply(true_rates, dfr_exponential))

set.seed(123)
exp_samp <- sampler(exp_sys)
exp_data <- data.frame(t = exp_samp(500), delta = rep(1L, 500))

exp_solver <- fit(exp_sys)
exp_result <- exp_solver(exp_data, par = c(0.15, 0.25, 0.25))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced

cat("True sum of rates:", sum(true_rates), "\n")
#> True sum of rates: 0.6
cat("Fitted sum of rates:", sum(coef(exp_result)), "\n")
#> Fitted sum of rates: 0.6033517
cat("Individual rates (unreliable):", coef(exp_result), "\n")
#> Individual rates (unreliable): 0.1344505 0.2344506 0.2344506
```

The sum is accurately estimated, but individual rates are not
meaningful.

### Mixed Types: Identifiable

When components have different distributional forms, the system is
generally identifiable because each component shapes the hazard
differently.

``` r
# Weibull (increasing hazard) + Exponential (constant hazard)
mixed <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01)
))

set.seed(42)
mixed_samp <- sampler(mixed)
mixed_data <- data.frame(t = mixed_samp(500), delta = rep(1L, 500))

mixed_solver <- fit(mixed)
mixed_result <- mixed_solver(mixed_data, par = c(1.5, 80, 0.02))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced

cat("True:   shape=2.00, scale=100.0, rate=0.010\n")
#> True:   shape=2.00, scale=100.0, rate=0.010
cat(sprintf("Fitted: shape=%.2f, scale=%.1f, rate=%.3f\n",
            coef(mixed_result)[1], coef(mixed_result)[2], coef(mixed_result)[3]))
#> Fitted: shape=2.49, scale=112.8, rate=0.013
```

## Handling Censored Data

In reliability studies, some units are removed from the test before
failing (right-censoring). The MLE correctly accounts for this.

``` r
set.seed(42)
model <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01)
))
true_par <- c(2, 100, 0.01)

# Generate true failure times
n <- 500
true_times <- sampler(model)(n, par = true_par)

# Right-censor at a fixed time horizon
censor_time <- 80
observed_t <- pmin(true_times, censor_time)
delta <- as.integer(true_times <= censor_time)
cens_data <- data.frame(t = observed_t, delta = delta)

cat("Proportion censored:", round(1 - mean(delta), 3), "\n")
#> Proportion censored: 0.232

# Fit with censored data
solver <- fit(model)
cens_result <- solver(cens_data, par = c(1.5, 80, 0.02))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced

cat("True parameters:    ", true_par, "\n")
#> True parameters:     2 100 0.01
cat("Fitted (censored):  ", round(coef(cens_result), 4), "\n")
#> Fitted (censored):   1.9106 105.7704 0.011
cat("Fitted (uncensored):", round(coef(result), 4), "\n")
#> Fitted (uncensored): 2.4743 109.1372 0.0121
```

With censoring, we lose information about what would have happened after
$t = 80$, which primarily affects the Weibull shape estimate (shape
governs late-life behavior). Still, the MLE correctly accounts for the
censoring mechanism and recovers parameters reasonably well. Larger
sample sizes or lower censoring proportions improve precision.

## Model Diagnostics

### Log-Likelihood Comparison

Compare competing models using log-likelihood, AIC, or BIC:

``` r
# Model 1: Weibull + Exponential (3 parameters)
m1 <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.01)
))

# Model 2: Two Weibulls (4 parameters)
m2 <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_weibull(shape = 1.5, scale = 200)
))

# Generate data from Model 1
set.seed(42)
data_m1 <- data.frame(
    t = sampler(m1)(300, par = c(2, 100, 0.01)),
    delta = rep(1L, 300)
)

r1 <- fit(m1)(data_m1, par = c(1.5, 80, 0.02))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
r2 <- fit(m2)(data_m1, par = c(1.5, 80, 1.2, 150))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced

k1 <- length(coef(r1))
k2 <- length(coef(r2))
ll1 <- logLik(r1)
ll2 <- logLik(r2)

aic1 <- -2 * ll1 + 2 * k1
aic2 <- -2 * ll2 + 2 * k2

cat(sprintf("Model 1 (Weibull+Exp):    LL=%.2f, k=%d, AIC=%.2f\n", ll1, k1, aic1))
#> Model 1 (Weibull+Exp):    LL=-1475.16, k=3, AIC=2956.33
cat(sprintf("Model 2 (Weibull+Weibull): LL=%.2f, k=%d, AIC=%.2f\n", ll2, k2, aic2))
#> Model 2 (Weibull+Weibull): LL=-1474.34, k=4, AIC=2956.68
cat("Preferred model (lower AIC):", ifelse(aic1 < aic2, "Model 1", "Model 2"), "\n")
#> Preferred model (lower AIC): Model 1
```

### Model Assumptions

Review what assumptions your model makes:

``` r
assumptions(m1)
#> [1] "Series system: system fails when any component fails"            
#> [2] "Component independence: component lifetimes are independent"     
#> [3] "Non-negative hazard: h_j(t) >= 0 for all j, t > 0"               
#> [4] "Cumulative hazard diverges: lim(t->Inf) H_sys(t) = Inf"          
#> [5] "Support is positive reals: t in (0, Inf)"                        
#> [6] "Observations are independent"                                    
#> [7] "Censoring indicator: 1=exact, 0=right-censored, -1=left-censored"
#> [8] "Non-informative censoring"
```

## Confidence Intervals

### Wald Confidence Intervals

Use the variance-covariance matrix from
[`vcov()`](https://rdrr.io/r/stats/vcov.html) to construct Wald-type
confidence intervals:

``` r
result <- fit(m1)(data_m1, par = c(1.5, 80, 0.02))
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
#> Warning in log(h_exact): NaNs produced
est <- coef(result)
se <- sqrt(diag(vcov(result)))

alpha <- 0.05
z <- qnorm(1 - alpha / 2)

ci_lower <- est - z * se
ci_upper <- est + z * se

ci_table <- data.frame(
    Parameter = c("shape", "scale", "rate"),
    Estimate = est,
    SE = se,
    Lower = ci_lower,
    Upper = ci_upper
)
print(ci_table, digits = 4)
#>   Parameter  Estimate       SE     Lower     Upper
#> 1     shape   2.47432  0.50039  1.493569   3.45507
#> 2     scale 109.13719 11.80866 85.992646 132.28174
#> 3      rate   0.01211  0.00202  0.008153   0.01607
```

### Delta Method for MTBF

The mean time between failures (MTBF) is a function of the parameters.
The delta method provides approximate confidence intervals for functions
of parameters.

``` r
# For this system, MTBF is not available in closed form, but we can
# approximate it numerically and use the delta method
par_hat <- coef(result)
V <- vcov(result)

# Numerical MTBF estimate via integration of survival function
S_fn <- surv(m1)

# Wrap S_fn for integrate() which may pass vectorized t values
compute_mtbf <- function(p) {
    integrate(function(t) sapply(t, function(ti) S_fn(ti, par = p)),
              0, Inf, subdivisions = 1000L)$value
}

mtbf <- compute_mtbf(par_hat)
cat("Estimated MTBF:", round(mtbf, 2), "\n")
#> Estimated MTBF: 53.75

# Delta method: gradient of MTBF w.r.t. parameters
eps <- 1e-5
grad_mtbf <- numeric(length(par_hat))
for (k in seq_along(par_hat)) {
    par_plus <- par_minus <- par_hat
    par_plus[k] <- par_hat[k] + eps
    par_minus[k] <- par_hat[k] - eps
    grad_mtbf[k] <- (compute_mtbf(par_plus) - compute_mtbf(par_minus)) / (2 * eps)
}

se_mtbf <- sqrt(t(grad_mtbf) %*% V %*% grad_mtbf)
cat(sprintf("MTBF = %.2f +/- %.2f (95%% CI: [%.2f, %.2f])\n",
            mtbf, 1.96 * se_mtbf, mtbf - 1.96 * se_mtbf, mtbf + 1.96 * se_mtbf))
#> MTBF = 53.75 +/- 4.39 (95% CI: [49.35, 58.14])
```
