# Tests for edge cases, boundary conditions, and previously uncovered paths.

test_that("log-likelihood handles left-censored observations", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ll <- loglik(sys)
    df_left <- make_left_censored_data(c(1, 2, 3))
    val <- ll(df_left, par = c(0.1, 0.2, 0.3))
    expect_true(is.finite(val))
    expect_true(val < 0)
})

test_that("log-likelihood with mixed exact/right/left censoring", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ll <- loglik(sys)
    df <- make_full_mixed_data(c(1, 2), c(5, 10), c(0.5, 1.5))
    val <- ll(df, par = c(0.1, 0.2, 0.3))
    expect_true(is.finite(val))
    expect_true(val < 0)
})

test_that("left-censored log-likelihood matches reference for exponential", {
    # For left-censored: ll_i = log(F(t_i)) = log(1 - exp(-lambda * t_i))
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    ll_sys <- loglik(sys)
    ll_ref <- loglik(ref)

    df_left <- make_left_censored_data(c(1, 2, 3))
    expect_equal(ll_sys(df_left, par = c(0.1, 0.2, 0.3)),
                 ll_ref(df_left, par = 0.6), tolerance = 1e-10)
})

test_that("hazard and survival at t=0", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    h <- hazard(sys)
    S <- surv(sys)
    Fc <- cdf(sys)

    # Exponential hazard is constant, so h(0) = 0.6
    expect_equal(h(0), 0.6, tolerance = 1e-10)
    # S(0) = 1
    expect_equal(S(0), 1.0, tolerance = 1e-10)
    # F(0) = 0
    expect_equal(Fc(0), 0.0, tolerance = 1e-10)
})

test_that("Weibull series at t=0", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))
    h <- hazard(sys)
    S <- surv(sys)

    # Weibull(k>1) has h(0)=0, S(0)=1
    expect_equal(h(0), 0.0, tolerance = 1e-10)
    expect_equal(S(0), 1.0, tolerance = 1e-10)
})

test_that("sample_components error when par is NULL everywhere", {
    # Construct with explicit n_par and par, then NULL out par
    sys <- dfr_dist_series(
        list(dfr_exponential(0.1), dfr_exponential(0.2)),
        par = NULL,
        n_par = c(1L, 1L)
    )
    # The system has NULL par, but components have par
    # sample_components falls through to component par
    # To trigger the error, we need both sys$par and par argument to be NULL
    # AND component par to be NULL -- construct without par
    sys2 <- dfr_dist_series(
        list(dfr_exponential(), dfr_exponential()),
        par = c(0.1, 0.2),
        n_par = c(1L, 1L)
    )
    # Now NULL out the system par to trigger the error path
    sys2$par <- NULL
    expect_error(sample_components(sys2, n = 10), "Parameters required")
})

test_that("constructor rejects zero-parameter components", {
    expect_error(
        dfr_dist_series(
            list(dfr_exponential(0.1)),
            n_par = c(0L)
        ),
        "at least one parameter"
    )
})

test_that("vectorized cum_haz works with analytical path", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    H <- cum_haz(sys)

    # Vectorized call
    t_vec <- c(1, 5, 10, 20)
    result <- H(t_vec)
    expected <- 0.6 * t_vec
    expect_equal(result, expected, tolerance = 1e-10)
    expect_equal(length(result), length(t_vec))
})
