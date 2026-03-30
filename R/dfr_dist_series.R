#' @keywords internal
#' @details
#' The \pkg{serieshaz} package composes multiple dynamic failure rate
#' (\code{\link[flexhaz]{dfr_dist}}) distributions into a series system
#' distribution. A series system fails when \emph{any} component fails, so the
#' system hazard is the sum of component hazards:
#' \deqn{h_{sys}(t) = \sum_{j=1}^{m} h_j(t, \theta_j)}
#' and the system survival is the product of component survivals:
#' \deqn{S_{sys}(t) = \prod_{j=1}^{m} S_j(t, \theta_j)}
#'
#' The series system object inherits from \code{dfr_dist}, which in turn
#' inherits from \code{likelihood_model}, \code{univariate_dist}, and
#' \code{dist}. This means all existing methods --- hazard, survival, CDF,
#' density, quantile function, sampling, log-likelihood, score, Hessian, and
#' MLE fitting --- work automatically on series systems without
#' reimplementation.
#'
#' Parameters across all components are stored as a single flat vector, with a
#' \emph{layout} that maps global indices to per-component indices. This
#' design enables standard optimizers (e.g., \code{\link[stats]{optim}}) to
#' work directly on the concatenated parameter vector.
#'
#' @section Package functions:
#' \describe{
#'   \item{\code{\link{dfr_dist_series}}}{Constructor: compose components into
#'     a series system}
#'   \item{\code{\link{is_dfr_dist_series}}}{Type predicate}
#'   \item{\code{\link{ncomponents}}}{Number of components}
#'   \item{\code{\link{component}}}{Extract a single component}
#'   \item{\code{\link{param_layout}}}{Parameter index mapping}
#'   \item{\code{\link{component_hazard}}}{Component-level hazard closure}
#'   \item{\code{\link{sample_components}}}{Sample component lifetimes}
#' }
#'
#' @seealso
#' \code{\link{dfr_dist_series}} for the constructor,
#' \code{\link[flexhaz]{dfr_dist}} for the parent class,
#' \code{\link[algebraic.dist]{hazard}} for distribution generics,
#' \code{\link[likelihood.model]{loglik}} for statistical inference generics
#'
#' \code{vignette("series-overview")} for a quick-start guide,
#' \code{vignette("series-math")} for mathematical foundations,
#' \code{vignette("series-fitting")} for MLE fitting workflows
"_PACKAGE"

#' Series System Distribution from DFR Components
#'
#' Composes m \code{dfr_dist} component distributions into a series system
#' distribution. A series system fails when any component fails, so the
#' system hazard is the sum of component hazards:
#' \eqn{h_{sys}(t) = \sum_j h_j(t)}.
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
#' \strong{Score and Hessian}: Both use a decomposed per-component approach
#' rather than \code{numDeriv} on the full system log-likelihood. When
#' components provide analytical \code{score_fn}/\code{hess_fn}, the
#' cumulative hazard derivatives are computed analytically via the
#' all-censored trick; the hazard derivatives use \code{numDeriv::jacobian}
#' (score) or \code{numDeriv::hessian} (Hessian) per-component. The Hessian
#' exploits block structure: cross-component blocks use only the rate
#' Jacobians already computed for the score.
#'
#' \strong{Identifiability}: Exponential series systems are \emph{not}
#' identifiable from system-level data alone --- only the sum of rates is
#' identifiable. When fitting to data, check \code{sum(coef(result))} rather
#' than individual rate parameters. Mixed-type series systems (e.g., Weibull +
#' Gompertz) are generally identifiable because the components have different
#' hazard shapes.
#'
#' \strong{Nested series}: A \code{dfr_dist_series} is itself a
#' \code{dfr_dist}, so it can be used as a component in another series system.
#' The resulting nested system's hazard is the sum of all leaf-component
#' hazards.
#'
#' \strong{Class hierarchy}: \code{dfr_dist_series} inherits from
#' \code{dfr_dist} -> \code{likelihood_model} -> \code{univariate_dist} ->
#' \code{dist}. All methods from these parent classes work automatically.
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' # --- Basic exponential series ---
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
#' # System survival at t = 5
#' S <- surv(sys)
#' S(5)   # exp(-0.6 * 5)
#'
#' # --- Mixed Weibull + Gompertz series ---
#' sys2 <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_gompertz(a = 0.01, b = 0.1)
#' ))
#' h2 <- hazard(sys2)
#' h2(50)  # sum of Weibull and Gompertz hazards at t=50
#'
#' # --- Nested series ---
#' subsystem <- dfr_dist_series(list(
#'     dfr_exponential(0.05),
#'     dfr_exponential(0.10)
#' ))
#' full_system <- dfr_dist_series(list(
#'     subsystem,
#'     dfr_weibull(shape = 2, scale = 200)
#' ))
#'
#' # --- Fitting workflow ---
#' solver <- fit(sys)
#' # result <- solver(df, par = c(0.1, 0.2, 0.3))
#' # coef(result)   # fitted parameters
#' # vcov(result)   # variance-covariance matrix
#' # logLik(result) # maximized log-likelihood
#' }
#'
#' @seealso
#' \code{\link{is_dfr_dist_series}} for the type predicate,
#' \code{\link{ncomponents}} and \code{\link{component}} for introspection,
#' \code{\link{param_layout}} for parameter index mapping,
#' \code{\link{component_hazard}} for per-component hazard closures,
#' \code{\link{sample_components}} for sampling component lifetimes,
#' \code{\link[flexhaz]{dfr_dist}} for the parent class constructor,
#' \code{\link[algebraic.dist]{hazard}} for distribution generics
#'
#' @family series system
#' @importFrom flexhaz dfr_dist is_dfr_dist
#' @importFrom numDeriv jacobian grad hessian
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
    if (any(n_par < 1L))
        stop("Each component must have at least one parameter")

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
    # Supports vectorized t input via numeric(length(t)) initialization
    sys_rate <- function(t, par, ...) {
        h <- numeric(length(t))
        for (j in seq_len(m)) {
            h <- h + components[[j]]$rate(t, par[layout[[j]]], ...)
        }
        h
    }

    # Build analytical cumulative hazard if ALL components provide one
    has_cum_haz <- vapply(
        components, function(comp) !is.null(comp$cum_haz_rate), logical(1))
    sys_cum_haz <- if (all(has_cum_haz)) {
        function(t, par, ...) {
            H <- numeric(length(t))
            for (j in seq_len(m)) {
                H <- H + components[[j]]$cum_haz_rate(t, par[layout[[j]]], ...)
            }
            H
        }
    } else {
        NULL
    }

    # Build decomposed system score function.
    # Instead of numDeriv on the full system log-likelihood (O(P * n * m)),
    # decompose into per-component problems (O(P * n)):
    #   score_k = sum_{exact}(dh_k/dtheta_k / h_sys(t_i)) - sum(dH_k/dtheta_k)
    # The cumulative hazard derivative uses the component's analytical score_fn
    # (via all-censored trick) when available; otherwise numDeriv per-component.
    # Hazard derivative uses numDeriv::jacobian on component rate.
    sys_score_fn <- function(df, par, ob_col = "t", delta_col = "delta", ...) {
        t_obs <- df[[ob_col]]
        delta <- if (delta_col %in% names(df))
            df[[delta_col]] else rep(1L, nrow(df))
        n_obs <- length(t_obs)
        exact_idx <- which(delta == 1)
        n_exact <- length(exact_idx)

        # Compute system hazard at exact observation times
        if (n_exact > 0) {
            t_exact <- t_obs[exact_idx]
            h_sys_exact <- numeric(n_exact)
            for (j in seq_len(m)) {
                h_sys_exact <- h_sys_exact +
                    components[[j]]$rate(t_exact, par[layout[[j]]], ...)
            }
        }

        score_vec <- numeric(total_np)

        for (j in seq_len(m)) {
            par_j <- par[layout[[j]]]
            comp_j <- components[[j]]

            # Part 1: Cumulative hazard derivative = sum_i dH_j(t_i)/dtheta_j
            if (!is.null(comp_j$score_fn)) {
                # Analytical: score_fn(all_censored) = -sum dH_j/dtheta_j
                df_cens <- df
                df_cens[[delta_col]] <- rep(0L, n_obs)
                cum_haz_deriv <- -comp_j$score_fn(
                    df_cens, par_j,
                    ob_col = ob_col, delta_col = delta_col, ...)
            } else {
                # Numerical fallback per-component (lower-dim than full system)
                cum_haz_sum_fn <- function(p) {
                    if (!is.null(comp_j$cum_haz_rate)) {
                        sum(comp_j$cum_haz_rate(t_obs, p, ...))
                    } else {
                        sum(vapply(t_obs, function(ti) {
                            stats::integrate(
                                function(s) comp_j$rate(s, p, ...), 0, ti
                            )$value
                        }, numeric(1)))
                    }
                }
                cum_haz_deriv <- numDeriv::grad(cum_haz_sum_fn, par_j)
            }

            # Part 2: Hazard derivative for exact observations
            # = sum_{i:exact} dh_j(t_i)/dtheta_j / h_sys(t_i)
            if (n_exact > 0) {
                rate_at_exact <- function(p) comp_j$rate(t_exact, p, ...)
                dh_j <- numDeriv::jacobian(rate_at_exact, par_j)
                haz_deriv_term <- colSums(dh_j / h_sys_exact)
            } else {
                haz_deriv_term <- numeric(n_par[j])
            }

            score_vec[layout[[j]]] <- haz_deriv_term - cum_haz_deriv
        }

        score_vec
    }

    # Build decomposed system Hessian function.
    # Block structure:
    #   Within-component (k = k):
    #     sum_{exact} [d2h_k/dtheta_k^2 / h_sys - (dh_k)(dh_k)^T / h_sys^2]
    #     - sum_i d2H_k/dtheta_k^2
    #   Cross-component (k != l):
    #     -sum_{exact} (dh_k/dtheta_k)(dh_l/dtheta_l)^T / h_sys^2
    sys_hess_fn <- function(df, par, ob_col = "t", delta_col = "delta", ...) {
        t_obs <- df[[ob_col]]
        delta <- if (delta_col %in% names(df))
            df[[delta_col]] else rep(1L, nrow(df))
        n_obs <- length(t_obs)
        exact_idx <- which(delta == 1)
        n_exact <- length(exact_idx)

        hess_mat <- matrix(0, nrow = total_np, ncol = total_np)

        # Compute system hazard and per-component rate Jacobians at exact obs
        dh_list <- vector("list", m)  # rate Jacobians for cross-component terms
        if (n_exact > 0) {
            t_exact <- t_obs[exact_idx]
            h_sys_exact <- numeric(n_exact)
            for (j in seq_len(m)) {
                par_j <- par[layout[[j]]]
                h_sys_exact <- h_sys_exact +
                    components[[j]]$rate(t_exact, par_j, ...)
                rate_fn_j <- function(p) components[[j]]$rate(t_exact, p, ...)
                dh_list[[j]] <- numDeriv::jacobian(rate_fn_j, par_j)
            }
            inv_h <- 1 / h_sys_exact
            inv_h2 <- inv_h * inv_h
        }

        # Within-component blocks
        for (j in seq_len(m)) {
            par_j <- par[layout[[j]]]
            np_j <- n_par[j]
            comp_j <- components[[j]]
            idx_j <- layout[[j]]

            # Cumulative hazard Hessian: d2(sum_i H_j)/dtheta_j^2
            if (!is.null(comp_j$hess_fn)) {
                df_cens <- df
                df_cens[[delta_col]] <- rep(0L, n_obs)
                cum_haz_hess <- -comp_j$hess_fn(
                    df_cens, par_j,
                    ob_col = ob_col, delta_col = delta_col, ...)
            } else {
                cum_haz_sum_fn <- function(p) {
                    if (!is.null(comp_j$cum_haz_rate)) {
                        sum(comp_j$cum_haz_rate(t_obs, p, ...))
                    } else {
                        sum(vapply(t_obs, function(ti) {
                            stats::integrate(
                                function(s) comp_j$rate(s, p, ...), 0, ti
                            )$value
                        }, numeric(1)))
                    }
                }
                cum_haz_hess <- numDeriv::hessian(cum_haz_sum_fn, par_j)
            }

            if (n_exact > 0) {
                dh_j <- dh_list[[j]]

                # Rate Hessian: sum_i d2h_j/dtheta_j^2 / h_sys(t_i)
                # Compute per-observation Hessians of the rate function
                rate_hess_term <- matrix(0, np_j, np_j)
                for (i in seq_len(n_exact)) {
                    ti <- t_exact[i]
                    d2h_i <- numDeriv::hessian(
                        function(p) components[[j]]$rate(ti, p, ...), par_j)
                    rate_hess_term <- rate_hess_term + d2h_i * inv_h[i]
                }

                # Outer product: -sum_i (dh_j)(dh_j)^T / h_sys^2
                outer_term <- crossprod(dh_j * inv_h2, dh_j)

                hess_mat[idx_j, idx_j] <-
                    rate_hess_term - outer_term - cum_haz_hess
            } else {
                hess_mat[idx_j, idx_j] <- -cum_haz_hess
            }
        }

        # Cross-component blocks (k != l)
        if (n_exact > 0 && m > 1) {
            for (k in seq_len(m - 1)) {
                for (l in (k + 1):m) {
                    idx_k <- layout[[k]]
                    idx_l <- layout[[l]]
                    # -sum_i (dh_k)(dh_l)^T / h_sys^2
                    cross <- -crossprod(dh_list[[k]] * inv_h2, dh_list[[l]])
                    hess_mat[idx_k, idx_l] <- cross
                    hess_mat[idx_l, idx_k] <- t(cross)
                }
            }
        }

        hess_mat
    }

    # Construct as dfr_dist — inherits all methods automatically
    obj <- dfr_dist(
        rate = sys_rate,
        par = par,
        cum_haz_rate = sys_cum_haz,
        score_fn = sys_score_fn,
        hess_fn = sys_hess_fn
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
#' Returns \code{TRUE} if \code{x} inherits from \code{"dfr_dist_series"},
#' \code{FALSE} otherwise.
#'
#' @param x Object to test.
#' @return Logical scalar.
#'
#' @details
#' Since \code{dfr_dist_series} inherits from \code{dfr_dist}, an object
#' that passes \code{is_dfr_dist_series()} will also pass
#' \code{\link[flexhaz]{is_dfr_dist}()}. Use this function when you need to
#' distinguish series systems from ordinary \code{dfr_dist} objects.
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_exponential(0.2)
#' ))
#' is_dfr_dist_series(sys)  # TRUE
#' is_dfr_dist(sys)         # also TRUE (inherits dfr_dist)
#'
#' single <- dfr_exponential(0.5)
#' is_dfr_dist_series(single)  # FALSE
#' is_dfr_dist(single)         # TRUE
#'
#' is_dfr_dist_series(42)  # FALSE
#' }
#'
#' @seealso \code{\link{dfr_dist_series}} for the constructor,
#'   \code{\link[flexhaz]{is_dfr_dist}} for the parent class predicate
#' @family series system
#' @export
is_dfr_dist_series <- function(x) {
    inherits(x, "dfr_dist_series")
}

#' Print method for series system distributions
#'
#' Displays a human-readable summary of a series system distribution,
#' including the number of components, per-component parameter counts and
#' values, and the system hazard/survival formulas.
#'
#' @param x A \code{dfr_dist_series} object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#'
#' @details
#' The output includes:
#' \itemize{
#'   \item Header with the number of components
#'   \item One line per component showing its parameter count and current
#'     parameter values (or "unknown" if parameters are \code{NULL})
#'   \item The system hazard formula: \eqn{h_{sys}(t) = \sum_j h_j(t, \theta_j)}
#'   \item The system survival formula:
#'     \eqn{S_{sys}(t) = \prod_j S_j(t, \theta_j)}
#' }
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_weibull(shape = 2, scale = 100)
#' ))
#' print(sys)
#' # Series system distribution with 2 components
#' #   Component 1: 1 param(s) [0.1]
#' #   Component 2: 2 param(s) [2, 100]
#' # System hazard: h_sys(t) = sum_j h_j(t, theta_j)
#' # Survival: S_sys(t) = exp(-H_sys(t)) = prod_j S_j(t, theta_j)
#' }
#'
#' @seealso \code{\link{dfr_dist_series}} for the constructor
#' @family series system
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
#' Returns the statistical and structural assumptions underlying a series
#' system model, which are important for the validity of MLE-based inference.
#'
#' @param model A \code{dfr_dist_series} object.
#' @param ... Additional arguments (unused).
#' @return Character vector of model assumptions.
#'
#' @details
#' The assumptions returned are:
#' \itemize{
#'   \item \strong{Series structure}: The system fails when any component
#'     fails (weakest-link model)
#'   \item \strong{Component independence}: Component lifetimes are
#'     statistically independent
#'   \item \strong{Non-negative hazard}: Each component hazard satisfies
#'     \eqn{h_j(t) \geq 0} for all \eqn{t > 0}
#'   \item \strong{Proper distribution}: The cumulative hazard diverges,
#'     ensuring \eqn{S_{sys}(t) \to 0} as \eqn{t \to \infty}
#'   \item \strong{Positive support}: The time domain is \eqn{(0, \infty)}
#'   \item \strong{Independent observations}: The observed lifetimes are
#'     independent
#'   \item \strong{Censoring convention}: \code{delta = 1} for exact,
#'     \code{0} for right-censored, \code{-1} for left-censored
#'   \item \strong{Non-informative censoring}: The censoring mechanism
#'     carries no information about the failure process
#' }
#'
#' These assumptions are required for the MLE fitting procedure
#' (\code{\link[generics]{fit}}) to produce valid estimates. Violation of
#' component independence, in particular, invalidates the hazard-sum property
#' that defines series systems.
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_weibull(shape = 2, scale = 100)
#' ))
#' assumptions(sys)
#' }
#'
#' @seealso
#' \code{\link[likelihood.model]{assumptions}} for the generic,
#' \code{\link{dfr_dist_series}} for the constructor,
#' \code{vignette("series-fitting")} for how assumptions affect inference
#' @family series system
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
