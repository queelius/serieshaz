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


test_that("dist.structure critical_states default dispatches via inheritance", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  # For a series system, component j is critical iff all other
  # components are functioning (sum_others = m - 1).
  crit <- dist.structure::critical_states(sys, j = 1L)
  expect_equal(nrow(crit), 1L)
  expect_equal(ncol(crit), 2L)
  expect_equal(as.integer(crit[1, ]), c(1L, 1L))
})


test_that("reliability equals p^m for series at iid p", {
  sys <- dfr_dist_series(list(
    dfr_exponential(1),
    dfr_exponential(1),
    dfr_exponential(1),
    dfr_exponential(1)
  ))
  for (p in c(0.1, 0.5, 0.9)) {
    expect_equal(dist.structure::reliability(sys, p), p^4,
                 tolerance = 1e-12)
  }
})


test_that("dual of series behaves like parallel via phi identity", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3),
    dfr_exponential(0.2)
  ))
  dsys <- dist.structure::dual(sys)
  # phi_dual(state) = 1 - phi(1 - state). Series phi is AND, so
  # dual phi is OR (parallel behavior).
  expect_equal(dist.structure::phi(dsys, c(0, 0, 0)), 0L)
  expect_equal(dist.structure::phi(dsys, c(1, 0, 0)), 1L)
  expect_equal(dist.structure::phi(dsys, c(1, 1, 1)), 1L)
})


test_that("is_coherent confirms series-specific topology", {
  sys <- dfr_dist_series(list(
    dfr_exponential(0.5),
    dfr_exponential(0.3)
  ))
  expect_true(dist.structure::is_coherent(sys))
})


test_that("importance measures work on dfr_dist_series objects", {
  sys <- dfr_dist_series(list(
    dfr_exponential(1),
    dfr_exponential(1),
    dfr_exponential(1)
  ))
  # Birnbaum: for series, I_B(j; p) = prod(p[-j]) (= p^(m-1) at iid).
  for (p in c(0.5, 0.8)) {
    expect_equal(dist.structure::birnbaum_importance(sys, 1L, p),
                 p^2, tolerance = 1e-12)
  }
  # Criticality at time t for series with iid Exp(1):
  # S(t) = exp(-t); I_B * F_j / F_sys = p^(m-1) * (1-p) / (1 - p^m).
  t0 <- 0.5
  p <- exp(-t0)
  expected <- p^2 * (1 - p) / (1 - p^3)
  expect_equal(dist.structure::criticality_importance(sys, 1L, t0),
               expected, tolerance = 1e-8)
})


test_that("substitute_component returns a coherent_dist with modified component", {
  sys <- dfr_dist_series(list(
    dfr_exponential(1),
    dfr_exponential(1)
  ))
  # Swap component 2 with a Weibull dist from algebraic.dist.
  new_comp <- algebraic.dist::weibull_dist(shape = 2, scale = 1)
  sys2 <- dist.structure::substitute_component(sys, 2L, new_comp)
  expect_true(dist.structure::is_dist_structure(sys2))
  c2 <- dist.structure::component(sys2, 2L)
  expect_s3_class(c2, "weibull_dist")
})


test_that("compose_systems nests series-of-series into a bigger series", {
  inner_a <- dfr_dist_series(list(dfr_exponential(1), dfr_exponential(1)))
  inner_b <- dfr_dist_series(list(dfr_exponential(1), dfr_exponential(1)))
  outer_comps <- lapply(1:2, function(i) algebraic.dist::exponential(1))
  outer <- dist.structure::series_dist(outer_comps)
  composed <- dist.structure::compose_systems(outer, list(inner_a, inner_b))
  expect_equal(dist.structure::ncomponents(composed), 4L)
  paths <- dist.structure::min_paths(composed)
  expect_equal(length(paths), 1L)
  expect_equal(sort(paths[[1]]), 1:4)
})


test_that("serieshaz fit+coef result can be replayed through dist.structure", {
  # Stage 3 pattern: fit with serieshaz, build DGP via algebraic.dist
  # + dist.structure::series_dist for post-fit topology queries.
  set.seed(1)
  sys <- dfr_dist_series(list(
    dfr_exponential(1),
    dfr_exponential(1)
  ))
  n <- 100
  samp <- algebraic.dist::sampler(sys)(n)
  df <- data.frame(t = samp)
  # Exponential series is not identifiable from system-level data; only
  # sum(rates) is. The Hessian may be near-singular, producing a
  # known warning from likelihood.model::fisher_mle. Suppress it here.
  fit_result <- suppressWarnings(fit(sys)(df, par = c(1, 1)))
  # Verify fit produced finite coefficients.
  rates_hat <- stats::coef(fit_result)
  expect_length(rates_hat, 2L)
  expect_true(all(is.finite(rates_hat)))
  # Replay as pure dist.structure: build a series of exponentials from the
  # fitted rates. Note: exp_series gives closed-form surv; series_dist gives
  # the general path.
  sys_replayed <- dist.structure::exp_series(rates_hat)
  expect_true(dist.structure::is_dist_structure(sys_replayed))
})
