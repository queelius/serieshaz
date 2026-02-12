#' @describeIn ncomponents Number of components in a series system.
#' @export
ncomponents.dfr_dist_series <- function(x, ...) {
    x$m
}

#' @describeIn component Extract component j as a standalone \code{dfr_dist}
#'   with its current parameters from the series system's parameter vector.
#' @export
component.dfr_dist_series <- function(x, j, ...) {
    if (j < 1L || j > x$m)
        stop(sprintf("Component index j=%d out of range [1, %d]", j, x$m))

    comp <- x$components[[j]]

    # Update component's par from the series system's current par vector
    if (!is.null(x$par)) {
        comp$par <- x$par[x$layout[[j]]]
    }
    comp
}

#' @describeIn param_layout Parameter index mapping for a series system.
#' @export
param_layout.dfr_dist_series <- function(x, ...) {
    x$layout
}

#' @describeIn component_hazard Hazard closure for component j of a series
#'   system. Returns \code{function(t, par_j = NULL, ...)} where \code{par_j}
#'   are the component-local parameters.
#' @importFrom algebraic.dist hazard
#' @export
component_hazard.dfr_dist_series <- function(x, j, ...) {
    if (j < 1L || j > x$m)
        stop(sprintf("Component index j=%d out of range [1, %d]", j, x$m))

    comp <- x$components[[j]]
    default_par <- if (!is.null(x$par)) x$par[x$layout[[j]]] else comp$par

    function(t, par = NULL, ...) {
        if (is.null(par)) par <- default_par
        comp$rate(t, par, ...)
    }
}

#' @describeIn sample_components Sample component lifetimes from a series
#'   system. Returns an n x m matrix where column j holds samples from
#'   component j. The system lifetime is \code{apply(mat, 1, min)}.
#' @importFrom algebraic.dist sampler
#' @export
sample_components.dfr_dist_series <- function(x, n, par = NULL, ...) {
    if (is.null(par)) par <- x$par
    if (is.null(par))
        stop("Parameters required: provide via 'par' or set in components")

    mat <- matrix(NA_real_, nrow = n, ncol = x$m)
    for (j in seq_len(x$m)) {
        par_j <- par[x$layout[[j]]]
        samp_j <- sampler(x$components[[j]])
        mat[, j] <- samp_j(n, par = par_j, ...)
    }
    colnames(mat) <- paste0("comp", seq_len(x$m))
    mat
}
