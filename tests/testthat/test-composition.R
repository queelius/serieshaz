# Test that series composition correctly produces system-level distributions.
# Gold standard: 3 exponentials with rates (r1, r2, r3) must behave
# identically to a single exponential with rate (r1 + r2 + r3).

test_that("exponential series hazard equals sum of component rates", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    h_sys <- hazard(sys)
    h_ref <- hazard(ref)

    for (t in c(0.1, 1, 5, 10, 100)) {
        expect_equal(h_sys(t), h_ref(t), tolerance = 1e-10)
    }
})

test_that("exponential series survival matches single exponential", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    S_sys <- surv(sys)
    S_ref <- surv(ref)

    for (t in c(0.5, 1, 5, 10, 50)) {
        expect_equal(S_sys(t), S_ref(t), tolerance = 1e-10)
    }
})

test_that("exponential series CDF matches single exponential", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    F_sys <- cdf(sys)
    F_ref <- cdf(ref)

    for (t in c(0.5, 1, 5, 10)) {
        expect_equal(F_sys(t), F_ref(t), tolerance = 1e-10)
    }
})

test_that("exponential series density matches single exponential", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    f_sys <- density(sys)
    f_ref <- density(ref)

    for (t in c(0.5, 1, 5, 10)) {
        expect_equal(f_sys(t), f_ref(t), tolerance = 1e-10)
    }
})

test_that("exponential series cumulative hazard matches single exponential", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    H_sys <- cum_haz(sys)
    H_ref <- cum_haz(ref)

    for (t in c(0.5, 1, 5, 10)) {
        expect_equal(H_sys(t), H_ref(t), tolerance = 1e-10)
    }
})

test_that("exponential series log-likelihood matches single exponential", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    ll_sys <- loglik(sys)
    ll_ref <- loglik(ref)

    df <- make_exact_data(c(1, 2, 3, 0.5))
    expect_equal(ll_sys(df, par = c(0.1, 0.2, 0.3)),
                 ll_ref(df, par = 0.6), tolerance = 1e-10)

    # With right-censored data
    df_mixed <- make_mixed_data(c(1, 2), c(5, 10))
    expect_equal(ll_sys(df_mixed, par = c(0.1, 0.2, 0.3)),
                 ll_ref(df_mixed, par = 0.6), tolerance = 1e-10)
})

test_that("Weibull series survival equals product of component survivals", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))

    S_sys <- surv(sys)
    S1 <- surv(dfr_weibull(shape = 2, scale = 100))
    S2 <- surv(dfr_weibull(shape = 1.5, scale = 200))

    for (t in c(10, 50, 100, 150)) {
        expected <- S1(t) * S2(t)
        expect_equal(S_sys(t), expected, tolerance = 1e-10)
    }
})

test_that("series system with single component equals that component", {
    comp <- dfr_weibull(shape = 2, scale = 100)
    sys <- dfr_dist_series(list(comp))

    h_sys <- hazard(sys)
    h_comp <- hazard(comp)

    S_sys <- surv(sys)
    S_comp <- surv(comp)

    for (t in c(10, 50, 100)) {
        expect_equal(h_sys(t), h_comp(t), tolerance = 1e-10)
        expect_equal(S_sys(t), S_comp(t), tolerance = 1e-10)
    }
})

test_that("parameter override works via par argument", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    h <- hazard(sys)

    # Default params -> 0.6
    expect_equal(h(1), 0.6, tolerance = 1e-10)

    # Override params -> 0.9
    expect_equal(h(1, par = c(0.3, 0.3, 0.3)), 0.9, tolerance = 1e-10)
})

test_that("analytical cum_haz_rate is used when all components provide it", {
    # All exponentials have analytical cum_haz_rate
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    expect_false(is.null(sys$cum_haz_rate))

    # Verify it's correct
    H <- cum_haz(sys)
    expect_equal(H(10), 6.0, tolerance = 1e-10)  # 0.6 * 10
})

test_that("vectorized t input works for hazard, surv, cdf, density", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    ref <- dfr_exponential(0.6)

    t_vec <- c(0.5, 1, 5, 10, 50)

    # Vectorized calls
    h_sys <- hazard(sys)
    S_sys <- surv(sys)
    F_sys <- cdf(sys)
    f_sys <- density(sys)

    h_ref <- hazard(ref)
    S_ref <- surv(ref)
    F_ref <- cdf(ref)
    f_ref <- density(ref)

    expect_equal(h_sys(t_vec), h_ref(t_vec), tolerance = 1e-10)
    expect_equal(S_sys(t_vec), S_ref(t_vec), tolerance = 1e-10)
    expect_equal(F_sys(t_vec), F_ref(t_vec), tolerance = 1e-10)
    expect_equal(f_sys(t_vec), f_ref(t_vec), tolerance = 1e-10)

    # Verify vectorized == element-wise scalar calls
    h_scalar <- vapply(t_vec, h_sys, numeric(1))
    expect_equal(h_sys(t_vec), h_scalar, tolerance = 1e-10)
})

test_that("mixed series with no cum_haz falls back to numerical", {
    # Create a component without analytical cum_haz_rate
    custom <- dfr_dist(
        rate = function(t, par, ...) par[1] * exp(-par[2] * t),
        par = c(0.5, 0.1)
    )
    sys <- dfr_dist_series(list(dfr_exponential(0.1), custom))

    # No analytical cum_haz for system
    expect_null(sys$cum_haz_rate)

    # But numerical integration still works
    H <- cum_haz(sys)
    # At t=0: H should be 0
    expect_equal(H(0), 0, tolerance = 1e-3)
    # At t > 0: H should be positive
    expect_true(H(5) > 0)
})
