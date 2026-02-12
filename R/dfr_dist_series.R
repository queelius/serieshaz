#' Series System Distribution from DFR Components
#'
#' Composes m \code{dfr_dist} component distributions into a series system
#' distribution. A series system fails when any component fails, so the
#' system hazard is the sum of component hazards: \eqn{h_{sys}(t) = \sum_j h_j(t)}.
#'
#' The resulting object inherits from \code{dfr_dist}, so all existing methods
#' (hazard, survival, CDF, density, sampling, log-likelihood, MLE fitting)
#' work automatically.
#'
#' @param components A list of \code{dfr_dist} objects representing system
#'   components.
#' @param par Optional concatenated parameter vector
#'   \eqn{\theta = (\theta_1, \ldots, \theta_m)}. If \code{NULL}, parameters
#'   are concatenated from component objects.
#' @param n_par Optional integer vector giving the number of parameters per
#'   component. Inferred from component \code{par} if not supplied; required
#'   when any component has \code{NULL} parameters.
#'
#' @return A \code{dfr_dist_series} object (inherits \code{dfr_dist}).
#'   Extra fields: \code{$components}, \code{$layout}, \code{$m}, \code{$n_par}.
#'
#' @details
#' \strong{Parameter layout}: Parameters are stored as a single concatenated
#' vector. The \code{$layout} field maps global indices to component indices.
#' For example, if component 1 has 2 parameters and component 2 has 1, then
#' \code{layout = list(1:2, 3)}.
#'
#' \strong{Analytical cumulative hazard}: If \emph{all} components provide
#' \code{cum_haz_rate}, the series system gets an analytical
#' \eqn{H_{sys}(t) = \sum_j H_j(t)}. Otherwise, falls back to numerical
#' integration.
#'
#' \strong{Score and Hessian}: Fall back to \code{numDeriv::grad} and
#' \code{numDeriv::hessian} on the (correct) composed log-likelihood.
#'
#' @examples
#' \donttest{
#' library(dfr.dist)
#'
#' # Three exponential components -> equivalent to single exponential
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_exponential(0.2),
#'     dfr_exponential(0.3)
#' ))
#' # System hazard = 0.6 (constant)
#' h <- hazard(sys)
#' h(10)  # 0.6
#'
#' # Mixed Weibull + Gompertz series system
#' sys2 <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_gompertz(a = 0.01, b = 0.1)
#' ))
#'
#' # Fit to data
#' solver <- fit(sys)
#' # result <- solver(df, par = c(0.1, 0.2, 0.3))
#' }
#'
#' @importFrom dfr.dist dfr_dist is_dfr_dist
#' @export
dfr_dist_series <- function(components, par = NULL, n_par = NULL) {
    stopifnot(is.list(components), length(components) >= 1L)
    for (j in seq_along(components)) {
        if (!is_dfr_dist(components[[j]]))
            stop("All components must be dfr_dist objects")
    }

    m <- length(components)

    # Determine parameter counts per component
    if (is.null(n_par)) {
        n_par <- vapply(components, function(comp) {
            if (is.null(comp$par))
                stop("Cannot infer parameter count: component has NULL par. ",
                     "Provide n_par or set par in each component.")
            length(comp$par)
        }, integer(1))
    }
    if (length(n_par) != m)
        stop(sprintf("n_par length (%d) must equal number of components (%d)",
                     length(n_par), m))

    # Compute parameter layout: which global indices map to each component
    layout <- vector("list", m)
    idx <- 1L
    for (j in seq_len(m)) {
        layout[[j]] <- seq.int(idx, length.out = n_par[j])
        idx <- idx + n_par[j]
    }
    total_np <- idx - 1L

    # Resolve parameters
    if (is.null(par)) {
        comp_pars <- lapply(components, function(comp) comp$par)
        if (!any(vapply(comp_pars, is.null, logical(1)))) {
            par <- unlist(comp_pars, use.names = FALSE)
        }
    }
    if (!is.null(par) && length(par) != total_np)
        stop(sprintf("Expected %d parameters but got %d",
                     total_np, length(par)))

    # Build composite hazard rate: h_sys(t) = sum_j h_j(t, theta_j)
    sys_rate <- function(t, par, ...) {
        h <- numeric(length(t))
        for (j in seq_len(m)) {
            par_j <- par[layout[[j]]]
            h <- h + components[[j]]$rate(t, par_j, ...)
        }
        h
    }

    # Build analytical cumulative hazard if ALL components provide one
    all_analytical <- all(vapply(
        components, function(comp) !is.null(comp$cum_haz_rate), logical(1)))
    sys_cum_haz <- if (all_analytical) {
        function(t, par, ...) {
            H <- 0
            for (j in seq_len(m)) {
                par_j <- par[layout[[j]]]
                H <- H + components[[j]]$cum_haz_rate(t, par_j, ...)
            }
            H
        }
    } else {
        NULL
    }

    # Construct as dfr_dist — inherits all methods automatically
    obj <- dfr_dist(
        rate = sys_rate,
        par = par,
        cum_haz_rate = sys_cum_haz,
        score_fn = NULL,
        hess_fn = NULL
    )

    # Attach series-specific metadata
    obj$components <- components
    obj$layout <- layout
    obj$m <- m
    obj$n_par <- n_par

    class(obj) <- c("dfr_dist_series", class(obj))
    obj
}

#' Test whether an object is a dfr_dist_series
#'
#' @param x Object to test.
#' @return Logical.
#' @export
is_dfr_dist_series <- function(x) {
    inherits(x, "dfr_dist_series")
}

#' Print method for series system distributions
#'
#' @param x A \code{dfr_dist_series} object.
#' @param ... Additional arguments (unused).
#' @export
print.dfr_dist_series <- function(x, ...) {
    cat(sprintf("Series system distribution with %d components\n", x$m))
    for (j in seq_len(x$m)) {
        np <- x$n_par[j]
        par_j <- if (!is.null(x$par)) x$par[x$layout[[j]]] else NULL
        par_str <- if (!is.null(par_j)) {
            paste(format(par_j, digits = 4), collapse = ", ")
        } else {
            "unknown"
        }
        cat(sprintf("  Component %d: %d param(s) [%s]\n", j, np, par_str))
    }
    cat("System hazard: h_sys(t) = sum_j h_j(t, theta_j)\n")
    cat("Survival: S_sys(t) = exp(-H_sys(t)) = prod_j S_j(t, theta_j)\n")
    invisible(x)
}

#' Assumptions for series system distributions
#'
#' @param model A \code{dfr_dist_series} object.
#' @param ... Additional arguments (unused).
#' @return Character vector of model assumptions.
#' @importFrom likelihood.model assumptions
#' @method assumptions dfr_dist_series
#' @export
assumptions.dfr_dist_series <- function(model, ...) {
    c(
        "Series system: system fails when any component fails",
        "Component independence: component lifetimes are independent",
        "Non-negative hazard: h_j(t) >= 0 for all j, t > 0",
        "Cumulative hazard diverges: lim(t->Inf) H_sys(t) = Inf",
        "Support is positive reals: t in (0, Inf)",
        "Observations are independent",
        "Censoring indicator: 1=exact, 0=right-censored, -1=left-censored",
        "Non-informative censoring"
    )
}
