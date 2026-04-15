# ==========================================================================
# dist.structure interop: dfr_dist_series implements the dist_structure
# protocol so users can call phi, min_paths, min_cuts, system_signature,
# structural_importance, reliability, dual, critical_states, importance
# measures (birnbaum, criticality, vesely-fussell), substitute_component,
# and compose_systems on serieshaz objects.
# ==========================================================================


library(flexhaz)


test_that("dfr_dist_series inherits dist_structure", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3)
  ))
  expect_true(dist.structure::is_dist_structure(sys))
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dfr_dist_series")
})


test_that("phi for series returns 1 iff every component functions", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  expect_equal(dist.structure::phi(sys, c(1, 1, 1)), 1L)
  expect_equal(dist.structure::phi(sys, c(1, 0, 1)), 0L)
  expect_equal(dist.structure::phi(sys, c(0, 0, 0)), 0L)
})


test_that("min_paths for series is a single path containing all components", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  paths <- dist.structure::min_paths(sys)
  expect_equal(length(paths), 1L)
  expect_equal(sort(paths[[1]]), 1:3)
})


test_that("min_cuts for series are singleton sets", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  cuts <- dist.structure::min_cuts(sys)
  expect_equal(length(cuts), 3L)
  expect_setequal(unlist(cuts), 1:3)
})


test_that("system_signature for series is (1, 0, ..., 0)", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  sig <- dist.structure::system_signature(sys)
  expect_equal(sig, c(1, 0, 0))
})


test_that("structural_importance for series matches prod of other p's at iid", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  imp <- dist.structure::structural_importance(sys, 1L)
  # Birnbaum at p = 0.5 each: structural importance = 0.5^(m-1) = 0.25
  expect_equal(imp, 0.25)
})


test_that("reliability for series at iid p equals p^m", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3)
  ))
  expect_equal(dist.structure::reliability(sys, 0.7), 0.49,
               tolerance = 1e-12)
})


test_that("birnbaum_importance for series equals product of other p's", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  p <- c(0.9, 0.8, 0.7)
  for (j in 1:3) {
    expected <- prod(p[-j])
    expect_equal(dist.structure::birnbaum_importance(sys, j, p),
                 expected, tolerance = 1e-10)
  }
})


test_that("flexhaz/serieshaz surv method wins dispatch over dist.structure default", {
  # The serieshaz/flexhaz surv computation should win over the
  # dist.structure default; dispatch order in the class chain is
  # c("dfr_dist_series", "dfr_dist", "dist_structure", ...), so
  # surv.dfr_dist (from flexhaz) is selected before surv.dist_structure.
  # System: Exp(0.5) and Exp(0.3) in series = Exp(0.8); S(1) = exp(-0.8).
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3)
  ))
  S <- algebraic.dist::surv(sys)
  expect_equal(S(1), exp(-0.8), tolerance = 1e-12)
  expect_equal(S(2.5), exp(-2), tolerance = 1e-12)
})


test_that("system_lifetime via dist.structure default works on serieshaz objects", {
  sys <- dfr_dist_series(list(
    dfr_exponential(1), dfr_exponential(1), dfr_exponential(1)
  ))
  # System fails at the FIRST component failure (series).
  expect_equal(dist.structure::system_lifetime(sys, c(2, 5, 1)), 1)
  expect_equal(dist.structure::system_lifetime(sys, c(0.5, 1.0, 0.7)), 0.5)
})


test_that("system_censoring labels each component correctly", {
  sys <- dfr_dist_series(list(
    dfr_exponential(1), dfr_exponential(1), dfr_exponential(1)
  ))
  cens <- dist.structure::system_censoring(sys, c(0.5, 1.0, 1.5))
  expect_equal(cens$system_time, 0.5)
  expect_equal(cens$component_status, c("exact", "right", "right"))
})


test_that("ncomponents and component still work after adoption", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3)
  ))
  expect_equal(ncomponents(sys), 2L)
  c1 <- component(sys, 1L)
  expect_true(is_dfr_dist(c1))
})
