# Tests for the decomposed system Hessian function.
# Verifies block-structured Hessian against finite differences of the score.

test_that("decomposed Hessian matches finite-difference of score for exponential series", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_exact_data(c(1, 2, 3, 0.5, 1.5))
    par <- c(0.1, 0.2, 0.3)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    # Finite-difference of score
    sc_fn <- score(sys)
    eps <- 1e-5
    fd_hess <- matrix(0, 3, 3)
    for (k in 1:3) {
        par_plus <- par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_hess[, k] <- (sc_fn(df, par = par_plus) - sc_fn(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(hess_val, fd_hess, tolerance = 1e-3)
})

test_that("decomposed Hessian is symmetric", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))
    df <- make_exact_data(c(50, 80, 120, 30, 90))
    par <- c(2, 100, 1.5, 200)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    expect_equal(hess_val, t(hess_val), tolerance = 1e-6)
    expect_equal(nrow(hess_val), 4)
})

test_that("decomposed Hessian matches finite-difference for mixed types", {
    sys <- dfr_dist_series(list(
        dfr_weibull(shape = 2, scale = 100),
        dfr_exponential(0.01)
    ))
    df <- make_exact_data(c(30, 50, 80, 40, 60))
    par <- c(2, 100, 0.01)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    sc_fn <- score(sys)
    eps <- 1e-5
    fd_hess <- matrix(0, 3, 3)
    for (k in 1:3) {
        par_plus <- par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_hess[, k] <- (sc_fn(df, par = par_plus) - sc_fn(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(hess_val, fd_hess, tolerance = 1e-2)
})

test_that("cross-component Hessian blocks are correct", {
    # For exponential series, the cross-component Hessian should be:
    #   d2ll/d(lambda_k)*d(lambda_l) = -n_exact / (sum_rates)^2
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_exact_data(c(1, 2, 3, 4, 5))
    par <- c(0.1, 0.2, 0.3)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    # For exponential: h_j(t) = lambda_j, constant
    # d(h_j)/d(lambda_j) = 1
    # Cross term: -sum(1 * 1 / h_sys^2) = -n_exact / 0.6^2
    expected_cross <- -5 / 0.6^2
    expect_equal(hess_val[1, 2], expected_cross, tolerance = 1e-3)
    expect_equal(hess_val[1, 3], expected_cross, tolerance = 1e-3)
    expect_equal(hess_val[2, 3], expected_cross, tolerance = 1e-3)
})

test_that("decomposed Hessian works with right-censored data", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_mixed_data(c(1, 2, 3), c(5, 10))
    par <- c(0.1, 0.2, 0.3)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    # Should be symmetric
    expect_equal(hess_val, t(hess_val), tolerance = 1e-6)

    # Verify against finite-difference of score
    sc_fn <- score(sys)
    eps <- 1e-5
    fd_hess <- matrix(0, 3, 3)
    for (k in 1:3) {
        par_plus <- par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_hess[, k] <- (sc_fn(df, par = par_plus) - sc_fn(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(hess_val, fd_hess, tolerance = 1e-3)
})

test_that("decomposed Hessian works with all-censored data", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_censored_data(c(1, 2, 3, 5))
    par <- c(0.1, 0.2, 0.3)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    # With all censored: no hazard term, only cumulative hazard
    # For exponential: d2H_j/d(lambda_j)^2 = 0 (H = lambda * t, linear)
    # So the Hessian should be all zeros
    expect_equal(hess_val, matrix(0, 3, 3), tolerance = 1e-6)
})

test_that("decomposed Hessian with Gompertz (has score_fn, no hess_fn)", {
    sys <- make_mixed_series()  # Weibull + Gompertz
    df <- make_exact_data(c(10, 20, 30, 50, 80))
    par <- c(2, 100, 0.01, 0.05)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    expect_equal(nrow(hess_val), 4)
    expect_equal(hess_val, t(hess_val), tolerance = 1e-6)

    # Verify against finite-difference of score
    sc_fn <- score(sys)
    eps <- 1e-5
    fd_hess <- matrix(0, 4, 4)
    for (k in 1:4) {
        par_plus <- par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_hess[, k] <- (sc_fn(df, par = par_plus) - sc_fn(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(hess_val, fd_hess, tolerance = 0.1)
})

test_that("decomposed Hessian with custom component (no score_fn, has cum_haz_rate)", {
    custom <- dfr_dist(
        rate = function(t, par, ...) par[1] * exp(-par[2] * t),
        par = c(0.5, 0.1),
        cum_haz_rate = function(t, par, ...) (par[1] / par[2]) * (1 - exp(-par[2] * t))
    )
    sys <- dfr_dist_series(list(dfr_exponential(0.1), custom))
    df <- make_exact_data(c(1, 3, 5))
    par <- c(0.1, 0.5, 0.1)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    expect_equal(nrow(hess_val), 3)
    expect_equal(hess_val, t(hess_val), tolerance = 1e-4)
})

test_that("decomposed Hessian with custom component (no score_fn, no cum_haz_rate)", {
    custom <- dfr_dist(
        rate = function(t, par, ...) par[1] * exp(-par[2] * t),
        par = c(0.5, 0.1)
    )
    sys <- dfr_dist_series(list(dfr_exponential(0.1), custom))
    df <- make_exact_data(c(1, 3, 5))
    par <- c(0.1, 0.5, 0.1)

    H_fn <- hess_loglik(sys)
    hess_val <- H_fn(df, par = par)

    expect_equal(nrow(hess_val), 3)
    expect_equal(hess_val, t(hess_val), tolerance = 1e-3)
})
