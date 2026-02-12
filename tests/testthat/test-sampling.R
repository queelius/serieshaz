# Test sampling from series systems.
# Key property: system lifetime = min(component lifetimes).

test_that("sample_components returns correct dimensions", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    mat <- sample_components(sys, n = 100)

    expect_equal(nrow(mat), 100)
    expect_equal(ncol(mat), 3)
    expect_equal(colnames(mat), c("comp1", "comp2", "comp3"))
})

test_that("sample_components allows parameter override", {
    sys <- make_exp_series(c(0.1, 0.2, 0.3))
    mat <- sample_components(sys, n = 50, par = c(1, 2, 3))

    expect_equal(nrow(mat), 50)
    expect_equal(ncol(mat), 3)
    # Higher rates -> shorter lifetimes on average
    expect_true(all(mat > 0))
})

test_that("sample_components requires parameters", {
    comps <- list(dfr_exponential(), dfr_exponential())
    # Can't construct without params
    expect_error(dfr_dist_series(comps))
})

test_that("system samples follow min(component lifetimes) distribution", {
    set.seed(42)
    sys <- make_exp_series(c(0.1, 0.2, 0.3))

    # Sample component lifetimes
    mat <- sample_components(sys, n = 5000)
    sys_times_from_components <- apply(mat, 1, min)

    # Sample directly from series system (via inverse CDF)
    samp <- sampler(sys)
    sys_times_direct <- samp(5000)

    # Both should follow Exp(0.6)
    # Compare means: E[T] = 1/0.6 ≈ 1.667
    expect_equal(mean(sys_times_from_components), 1 / 0.6,
                 tolerance = 0.1)
    expect_equal(mean(sys_times_direct), 1 / 0.6,
                 tolerance = 0.1)

    # KS test: both should be Exp(0.6)
    ks1 <- ks.test(sys_times_from_components, "pexp", rate = 0.6)
    ks2 <- ks.test(sys_times_direct, "pexp", rate = 0.6)
    expect_true(ks1$p.value > 0.01)
    expect_true(ks2$p.value > 0.01)
})

test_that("Weibull series samples have correct survival distribution", {
    set.seed(123)
    sys <- make_weibull_series(shapes = c(2, 2), scales = c(100, 100))

    # System survival: S_sys(t) = S1(t) * S2(t) = exp(-2*(t/100)^2)
    # This is Weibull(shape=2, scale=100/sqrt(2))
    samp <- sampler(sys)
    samples <- samp(5000)

    effective_scale <- 100 / sqrt(2)
    ks <- ks.test(samples, "pweibull", shape = 2, scale = effective_scale)
    expect_true(ks$p.value > 0.01)
})
