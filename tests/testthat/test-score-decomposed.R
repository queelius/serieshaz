# Tests for the decomposed system score function.
# Verifies that the per-component score decomposition produces identical
# results to finite differences of the log-likelihood.

test_that("decomposed score matches finite-difference for exponential series", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_exact_data(c(1, 2, 3, 0.5, 1.5))
    par <- c(0.1, 0.2, 0.3)

    s <- score(sys)
    score_val <- s(df, par = par)

    # Manual finite difference of log-likelihood
    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-4)
})

test_that("decomposed score matches finite-difference for Weibull series", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))
    df <- make_exact_data(c(50, 80, 120, 30, 90))
    par <- c(2, 100, 1.5, 200)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(4)
    for (k in 1:4) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-3)
})

test_that("decomposed score matches finite-difference for mixed types", {
    sys <- dfr_dist_series(list(
        dfr_weibull(shape = 2, scale = 100),
        dfr_exponential(0.01)
    ))
    df <- make_exact_data(c(30, 50, 80, 40, 60))
    par <- c(2, 100, 0.01)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-3)
})

test_that("decomposed score works with right-censored data", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_mixed_data(c(1, 2, 3), c(5, 10))
    par <- c(0.1, 0.2, 0.3)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-4)
})

test_that("decomposed score works with all-censored data", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_censored_data(c(1, 2, 3, 5))
    par <- c(0.1, 0.2, 0.3)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-4)
})

test_that("decomposed score at MLE is approximately zero", {
    set.seed(42)
    sys <- make_exp_series(c(0.5, 0.5))

    samp <- sampler(sys)
    df <- make_exact_data(samp(500))

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(0.6, 0.6)))
    mle_par <- coef(result)

    s <- score(sys)
    score_at_mle <- s(df, par = mle_par)
    expect_true(all(abs(score_at_mle) < 0.5))
})

test_that("score with Gompertz components (has score_fn, no hess_fn)", {
    sys <- make_mixed_series()  # Weibull + Gompertz
    df <- make_exact_data(c(10, 20, 30, 50, 80))
    par <- c(2, 100, 0.01, 0.05)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(4)
    for (k in 1:4) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-3)
})

test_that("score with custom component (no score_fn, has cum_haz_rate)", {
    # Custom declining hazard -- no score_fn but has cum_haz_rate
    custom <- dfr_dist(
        rate = function(t, par, ...) par[1] * exp(-par[2] * t),
        par = c(0.5, 0.1),
        cum_haz_rate = function(t, par, ...) (par[1] / par[2]) * (1 - exp(-par[2] * t))
    )
    sys <- dfr_dist_series(list(dfr_exponential(0.1), custom))
    df <- make_exact_data(c(1, 3, 5, 8, 12))
    par <- c(0.1, 0.5, 0.1)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-3)
})

test_that("score with custom component (no score_fn, no cum_haz_rate)", {
    # Custom component without score_fn or cum_haz_rate -- exercises the
    # numerical integration fallback in the score function
    custom <- dfr_dist(
        rate = function(t, par, ...) par[1] * exp(-par[2] * t),
        par = c(0.5, 0.1)
    )
    sys <- dfr_dist_series(list(dfr_exponential(0.1), custom))
    df <- make_exact_data(c(1, 3, 5))
    par <- c(0.1, 0.5, 0.1)

    s <- score(sys)
    score_val <- s(df, par = par)

    ll <- loglik(sys)
    eps <- 1e-6
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par; par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-2)
})
