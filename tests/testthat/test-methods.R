# Tests for series-specific methods: component, ncomponents, param_layout,
# component_hazard, sample_components

test_that("ncomponents returns correct count", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    expect_equal(ncomponents(sys), 3L)

    sys2 <- make_weibull_series()
    expect_equal(ncomponents(sys2), 2L)
})

test_that("param_layout returns correct index mapping", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    lay <- param_layout(sys)

    expect_equal(length(lay), 3L)
    expect_equal(lay[[1]], 1L)
    expect_equal(lay[[2]], 2L)
    expect_equal(lay[[3]], 3L)

    # Weibull: 2 params each
    sys2 <- make_weibull_series()
    lay2 <- param_layout(sys2)
    expect_equal(lay2[[1]], 1:2)
    expect_equal(lay2[[2]], 3:4)
})

test_that("component extracts standalone dfr_dist with correct params", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))

    c1 <- component(sys, 1)
    expect_true(is_dfr_dist(c1))
    expect_equal(params(c1), 0.1)

    c2 <- component(sys, 2)
    expect_equal(params(c2), 0.2)

    c3 <- component(sys, 3)
    expect_equal(params(c3), 0.3)
})

test_that("component extracts Weibull with correct params", {
    sys <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))

    c1 <- component(sys, 1)
    expect_equal(params(c1), c(2, 100))

    c2 <- component(sys, 2)
    expect_equal(params(c2), c(1.5, 200))
})

test_that("component validates index bounds", {
    sys <- make_exp_series(c(0.1, 0.2))
    expect_error(component(sys, 0), "out of range")
    expect_error(component(sys, 3), "out of range")
})

test_that("component_hazard returns correct closure", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))

    h1 <- component_hazard(sys, 1)
    h2 <- component_hazard(sys, 2)
    h3 <- component_hazard(sys, 3)

    expect_equal(h1(5), 0.1, tolerance = 1e-10)
    expect_equal(h2(5), 0.2, tolerance = 1e-10)
    expect_equal(h3(5), 0.3, tolerance = 1e-10)
})

test_that("component_hazard allows parameter override", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    h1 <- component_hazard(sys, 1)

    expect_equal(h1(5, par = 0.5), 0.5, tolerance = 1e-10)
})

test_that("component_hazard validates index", {
    sys <- make_exp_series(c(0.1, 0.2))
    expect_error(component_hazard(sys, 3), "out of range")
})

test_that("is_dfr_dist_series identifies series objects", {
    sys <- make_exp_series(c(0.1, 0.2))
    expect_true(is_dfr_dist_series(sys))
    expect_true(is_dfr_dist(sys))  # also a dfr_dist

    comp <- dfr_exponential(0.1)
    expect_false(is_dfr_dist_series(comp))
})

test_that("print method works without error", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    expect_output(print(sys), "Series system distribution with 3 components")
    expect_output(print(sys), "Component 1")
})

test_that("assumptions returns series-specific assumptions", {
    sys <- make_exp_series(c(0.1, 0.2))
    a <- assumptions(sys)
    expect_true(any(grepl("series", a, ignore.case = TRUE)))
    expect_true(any(grepl("independence", a, ignore.case = TRUE)))
})

test_that("params returns concatenated parameter vector", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    expect_equal(params(sys), c(0.1, 0.2, 0.3))

    sys2 <- make_weibull_series(shapes = c(2, 1.5), scales = c(100, 200))
    # Params are concatenated per component: (shape1, scale1, shape2, scale2)
    expect_equal(params(sys2), c(2, 100, 1.5, 200))
})
