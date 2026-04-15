#' @method ncomponents dfr_dist_series
#' @export
ncomponents.dfr_dist_series <- function(x, ...) {
    x$m
}

#' @method component dfr_dist_series
#' @export
component.dfr_dist_series <- function(x, j, ...) {
    if (!is.numeric(j) || length(j) != 1L || is.na(j) || j != as.integer(j))
        stop("Component index j must be a single integer")
    j <- as.integer(j)
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
    if (!is.numeric(j) || length(j) != 1L || is.na(j) || j != as.integer(j))
        stop("Component index j must be a single integer")
    j <- as.integer(j)
    if (j < 1L || j > x$m)
        stop(sprintf("Component index j=%d out of range [1, %d]", j, x$m))

    comp <- x$components[[j]]

    function(t, par = NULL, ...) {
        if (is.null(par))
            par <- if (!is.null(x$par)) x$par[x$layout[[j]]] else comp$par
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


# ==========================================================================
# dist.structure topology methods (series specialization)
# ==========================================================================
#
# A series system has trivial topology: phi = AND, min_paths = {1..m},
# min_cuts = {{1}, ..., {m}}, signature = (1, 0, ..., 0). Registering
# specialized methods avoids the dist.structure default-method enumeration.
# All other dist.structure generics (critical_states, structural_importance,
# reliability, dual, is_coherent, system_lifetime, system_censoring) work
# automatically via dispatch through the dist_structure class.
# ==========================================================================


#' @method phi dfr_dist_series
#' @export
phi.dfr_dist_series <- function(x, state) {
    stopifnot(length(state) == x$m)
    as.integer(all(state == 1L))
}


#' @method min_paths dfr_dist_series
#' @export
min_paths.dfr_dist_series <- function(x) {
    list(seq_len(x$m))
}


#' @method min_cuts dfr_dist_series
#' @export
min_cuts.dfr_dist_series <- function(x) {
    lapply(seq_len(x$m), function(j) j)
}


#' @method system_signature dfr_dist_series
#' @export
system_signature.dfr_dist_series <- function(x) {
    c(1, rep(0, x$m - 1L))
}
