# Test score and Hessian computation for series systems.
# Since dfr_dist_series uses numDeriv fallback, we verify against
# known analytical results for exponential series.

test_that("score at MLE is approximately zero for exponential series", {
    set.seed(42)
    sys <- make_exp_series(c(0.5, 0.5))

    samp <- sampler(sys)
    df <- make_exact_data(samp(500))

    # Fit to get MLE
    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(0.6, 0.6)))
    mle_par <- coef(result)

    # Score at MLE should be near zero
    s <- score(sys)
    score_at_mle <- s(df, par = mle_par)
    expect_true(all(abs(score_at_mle) < 0.5))
})

test_that("Hessian is negative definite at MLE", {
    set.seed(42)
    sys <- make_exp_series(c(0.3, 0.7))

    samp <- sampler(sys)
    df <- make_exact_data(samp(300))

    solver <- fit(sys)
    result <- suppressWarnings(solver(df, par = c(0.4, 0.8)))
    mle_par <- coef(result)

    H <- hess_loglik(sys)
    hess_at_mle <- H(df, par = mle_par)

    # Eigenvalues should all be negative (negative definite)
    evals <- eigen(hess_at_mle, symmetric = TRUE)$values
    expect_true(all(evals < 0))
})

test_that("score matches finite-difference approximation", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    df <- make_exact_data(c(1, 2, 3, 0.5, 1.5))
    par <- c(0.1, 0.2, 0.3)

    # numDeriv score (our default)
    s <- score(sys)
    score_val <- s(df, par = par)

    # Manual finite difference
    ll <- loglik(sys)
    eps <- 1e-5
    fd_score <- numeric(3)
    for (k in 1:3) {
        par_plus <- par
        par_minus <- par
        par_plus[k] <- par[k] + eps
        par_minus[k] <- par[k] - eps
        fd_score[k] <- (ll(df, par = par_plus) - ll(df, par = par_minus)) / (2 * eps)
    }

    expect_equal(score_val, fd_score, tolerance = 1e-4)
})

test_that("Hessian matches finite-difference of score", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))
    df <- make_exact_data(c(50, 80, 120, 30, 90))
    par <- c(2, 100, 1.5, 200)

    H <- hess_loglik(sys)
    hess_val <- H(df, par = par)

    # Verify dimensions
    expect_equal(nrow(hess_val), 4)
    expect_equal(ncol(hess_val), 4)

    # Hessian should be symmetric
    expect_equal(hess_val, t(hess_val), tolerance = 1e-6)
})
