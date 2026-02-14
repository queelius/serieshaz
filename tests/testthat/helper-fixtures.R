# Test fixtures for serieshaz

library(flexhaz)

# ---- Series system fixtures ----

#' Three exponential components -> equivalent to single exponential
make_exp_series <- function(rates = c(0.1, 0.2, 0.3)) {
    comps <- lapply(rates, dfr_exponential)
    dfr_dist_series(comps)
}

#' Two Weibull components
make_weibull_series <- function(shapes = c(2, 1.5), scales = c(100, 200)) {
    comps <- list(
        dfr_weibull(shape = shapes[1], scale = scales[1]),
        dfr_weibull(shape = shapes[2], scale = scales[2])
    )
    dfr_dist_series(comps)
}

#' Mixed series: Weibull + Gompertz
make_mixed_series <- function() {
    comps <- list(
        dfr_weibull(shape = 2, scale = 100),
        dfr_gompertz(a = 0.01, b = 0.05)
    )
    dfr_dist_series(comps)
}

# ---- Data frame fixtures ----

make_exact_data <- function(times) {
    data.frame(t = times, delta = rep(1L, length(times)))
}

make_censored_data <- function(times) {
    data.frame(t = times, delta = rep(0L, length(times)))
}

make_mixed_data <- function(exact_times, censored_times) {
    data.frame(
        t = c(exact_times, censored_times),
        delta = c(rep(1L, length(exact_times)), rep(0L, length(censored_times)))
    )
}
