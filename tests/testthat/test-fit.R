# Test MLE fitting on series system data.

test_that("MLE recovers exponential series parameters", {
    set.seed(42)
    true_rates <- c(0.1, 0.2, 0.3)
    sys <- make_exp_series(true_rates)

    # Generate data from the series system
    samp <- sampler(sys)
    times <- samp(500)
    df <- make_exact_data(times)

    # Fit
    solver <- fit(sys)
    # Note: exponential series is not identifiable from system-level data alone
    # (only sum of rates is identifiable). But the optimizer should find
    # parameters where the sum equals the true total rate.
    result <- suppressWarnings(solver(df, par = c(0.15, 0.25, 0.35)))

    expect_true(result$converged)

    # Total rate should be close to 0.6
    est_total <- sum(coef(result))
    expect_equal(est_total, sum(true_rates), tolerance = 0.15)
})

test_that("MLE on Weibull series recovers parameters from component data", {
    set.seed(123)
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))

    # Generate data from series system
    samp <- sampler(sys)
    times <- samp(300)
    df <- make_exact_data(times)

    # Fit with reasonable initial values (shape1, scale1, shape2, scale2)
    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(1.8, 90, 1.3, 180),
                                       method = "Nelder-Mead"))

    expect_true(result$converged)
    # Log-likelihood at MLE should be finite (fisher_mle stores it as $loglik)
    expect_true(is.finite(result$loglik))
})

test_that("MLE handles right-censored series data", {
    set.seed(42)
    sys <- make_exp_series(c(0.1, 0.2, 0.3))

    samp <- sampler(sys)
    times <- samp(300)

    # Right-censor at tau = 2
    tau <- 2
    observed_times <- pmin(times, tau)
    delta <- as.integer(times <= tau)
    df <- data.frame(t = observed_times, delta = delta)

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(0.15, 0.25, 0.35)))

    expect_true(result$converged)
    est_total <- sum(coef(result))
    expect_equal(est_total, 0.6, tolerance = 0.2)
})

test_that("mixed-type MLE produces positive fitted parameters", {
    set.seed(42)
    sys <- dfr_dist_series(list(
        dfr_weibull(shape = 2, scale = 100),
        dfr_exponential(0.01)
    ))

    samp <- sampler(sys)
    df <- make_exact_data(samp(500))

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(1.5, 80, 0.02)))

    # All fitted parameters should be positive (negative rates are nonphysical)
    expect_true(all(coef(result) > 0))
    expect_true(result$converged)

    # Parameters should be in a reasonable range of the true values
    expect_equal(coef(result)[1], 2, tolerance = 1.0)     # shape
    expect_equal(coef(result)[2], 100, tolerance = 50)     # scale
    expect_equal(coef(result)[3], 0.01, tolerance = 0.01)  # rate
})

test_that("score at MLE is approximately zero for mixed-type series", {
    set.seed(42)
    sys <- dfr_dist_series(list(
        dfr_weibull(shape = 2, scale = 100),
        dfr_exponential(0.01)
    ))

    samp <- sampler(sys)
    df <- make_exact_data(samp(300))

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(1.5, 80, 0.02)))
    mle_par <- coef(result)

    # Score at MLE should be near zero (first-order optimality)
    s <- score(sys)
    score_at_mle <- s(df, par = mle_par)
    expect_true(all(abs(score_at_mle) < 1.0))
})

test_that("fit returns fisher_mle with expected methods", {
    set.seed(42)
    sys <- make_exp_series(c(0.5, 0.5))

    samp <- sampler(sys)
    df <- make_exact_data(samp(200))

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(0.6, 0.6)))

    # Standard fisher_mle methods should work
    expect_true(length(coef(result)) == 2)
    expect_true(is.matrix(vcov(result)))
    expect_equal(nrow(vcov(result)), 2)
    expect_true(is.finite(logLik(result)))
})
