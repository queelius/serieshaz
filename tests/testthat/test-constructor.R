# Test constructor validation and edge cases.

test_that("constructor validates component types", {
    expect_error(
        dfr_dist_series(list(1, 2, 3)),
        "dfr_dist"
    )
    expect_error(
        dfr_dist_series(list(dfr_exponential(0.1), "not a dist")),
        "dfr_dist"
    )
})

test_that("constructor requires non-empty component list", {
    expect_error(dfr_dist_series(list()))
})

test_that("constructor validates parameter length", {
    comps <- list(dfr_exponential(0.1), dfr_exponential(0.2))
    expect_error(
        dfr_dist_series(comps, par = c(0.1, 0.2, 0.3)),
        "Expected 2 parameters but got 3"
    )
})

test_that("constructor infers params from components", {
    sys <- dfr_dist_series(list(
        dfr_exponential(0.1),
        dfr_exponential(0.2)
    ))
    expect_equal(params(sys), c(0.1, 0.2))
})

test_that("constructor allows explicit par to override", {
    sys <- dfr_dist_series(
        list(dfr_exponential(0.1), dfr_exponential(0.2)),
        par = c(0.5, 0.5)
    )
    expect_equal(params(sys), c(0.5, 0.5))
})

test_that("constructor requires n_par for NULL-param components", {
    expect_error(
        dfr_dist_series(list(dfr_exponential(), dfr_exponential())),
        "Cannot infer parameter count"
    )
})

test_that("constructor works with n_par for NULL-param components", {
    sys <- dfr_dist_series(
        list(dfr_exponential(), dfr_exponential()),
        par = c(0.1, 0.2),
        n_par = c(1L, 1L)
    )
    expect_equal(params(sys), c(0.1, 0.2))
    expect_equal(ncomponents(sys), 2L)
})

test_that("n_par length must match number of components", {
    expect_error(
        dfr_dist_series(
            list(dfr_exponential(0.1), dfr_exponential(0.2)),
            n_par = c(1L, 1L, 1L)
        ),
        "n_par length"
    )
})

test_that("class hierarchy is correct", {
    sys <- make_exp_series()
    expect_s3_class(sys, "dfr_dist_series")
    expect_s3_class(sys, "dfr_dist")
    expect_s3_class(sys, "likelihood_model")
    expect_s3_class(sys, "univariate_dist")
    expect_s3_class(sys, "dist")
})

test_that("nested series systems work", {
    # Inner series: Exp(0.1) + Exp(0.2)
    inner <- dfr_dist_series(list(
        dfr_exponential(0.1),
        dfr_exponential(0.2)
    ))
    # Outer series: inner + Exp(0.3)
    outer <- dfr_dist_series(list(inner, dfr_exponential(0.3)))

    # Should behave like Exp(0.6)
    h <- hazard(outer)
    expect_equal(h(5), 0.6, tolerance = 1e-10)

    # Layout: inner has 2 params, outer component has 1
    expect_equal(param_layout(outer), list(1:2, 3L))
})
