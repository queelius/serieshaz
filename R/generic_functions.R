#' Get the parameter layout for a system
#'
#' Returns the mapping from global (flat) parameter indices to per-component
#' parameter indices, enabling the series system to distribute a single
#' parameter vector across its components.
#'
#' @param x A system object (e.g., \code{\link{dfr_dist_series}}).
#' @param ... Additional arguments passed to methods.
#' @return A list of integer vectors, one per component, containing global
#'   parameter indices.
#'
#' @details
#' Parameters across all components are stored as a single concatenated
#' vector \eqn{\theta = (\theta_1, \ldots, \theta_m)}. The layout maps
#' global indices back to each component. For example, with:
#' \itemize{
#'   \item Component 1: Weibull (shape, scale) --- 2 parameters
#'   \item Component 2: Exponential (rate) --- 1 parameter
#'   \item Component 3: Gompertz (a, b) --- 2 parameters
#' }
#' the layout is \code{list(1:2, 3, 4:5)}, so the global parameter vector
#' \code{c(shape1, scale1, rate2, a3, b3)} gets sliced as \code{par[1:2]}
#' for component 1, \code{par[3]} for component 2, and \code{par[4:5]} for
#' component 3.
#'
#' This design enables standard optimizers to work on a flat vector while
#' the series system internally distributes parameters to the correct
#' components.
#'
#' @examples
#'
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_exponential(0.05),
#'     dfr_gompertz(a = 0.01, b = 0.1)
#' ))
#' param_layout(sys)
#' # list(1:2, 3, 4:5)
#'
#' @seealso
#' \code{\link{component}} to extract a component with its parameters,
#' \code{\link[algebraic.dist]{params}} to get the full parameter vector,
#' \code{\link{dfr_dist_series}} for the constructor
#' @family system introspection
#' @export
param_layout <- function(x, ...) {
    UseMethod("param_layout")
}

#' Get the hazard function for a specific component
#'
#' Returns a closure that computes the hazard rate for component \code{j}
#' of a series system. Useful for plotting hazard decompositions and
#' understanding each component's contribution to system risk.
#'
#' @param x A system object (e.g., \code{\link{dfr_dist_series}}).
#' @param j Component index (integer, \code{1 <= j <= ncomponents(x)}).
#' @param ... Additional arguments passed to methods.
#' @return A closure \code{function(t, par = NULL, ...)} that evaluates
#'   component \code{j}'s hazard rate. If \code{par} is \code{NULL}, the
#'   component's default parameters (from the system) are used.
#'
#' @details
#' The returned closure evaluates \eqn{h_j(t, \theta_j)} for component
#' \code{j}. The \code{par} argument accepts \emph{component-local}
#' parameters (not the full system parameter vector). This is useful for:
#' \itemize{
#'   \item Plotting individual hazard contributions
#'   \item Verifying that \eqn{\sum_j h_j(t) = h_{sys}(t)}
#'   \item Sensitivity analysis on a single component
#' }
#'
#' @examples
#'
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_exponential(0.2)
#' ))
#'
#' h1 <- component_hazard(sys, 1)
#' h2 <- component_hazard(sys, 2)
#' h_sys <- hazard(sys)
#'
#' # Verify hazard sum property
#' t <- 10
#' h1(t) + h2(t)  # 0.3
#' h_sys(t)        # 0.3 (same!)
#'
#' @seealso
#' \code{\link{component}} to extract the full component object,
#' \code{\link[algebraic.dist]{hazard}} for the system-level hazard,
#' \code{\link{dfr_dist_series}} for the constructor
#' @family system introspection
#' @export
component_hazard <- function(x, j, ...) {
    UseMethod("component_hazard")
}

#' Sample component lifetimes from a system
#'
#' Generates an \eqn{n \times m} matrix where column \eqn{j} contains
#' independent samples from component \eqn{j}'s lifetime distribution.
#' The system lifetime is the row-wise minimum.
#'
#' @param x A system object (e.g., \code{\link{dfr_dist_series}}).
#' @param n Number of samples (rows).
#' @param par Optional parameter vector override.
#' @param ... Additional arguments passed to methods.
#' @return An \eqn{n \times m} numeric matrix of component lifetimes, with
#'   columns named \code{comp1}, \code{comp2}, etc.
#'
#' @details
#' Each column is sampled independently using the component's own sampler.
#' Since the series system fails when \emph{any} component fails, the system
#' lifetime for each observation is:
#'
#' \preformatted{
#' t_sys <- apply(mat, 1, min)
#' }
#'
#' The failing component for each observation can be identified via:
#'
#' \preformatted{
#' failing <- apply(mat, 1, which.min)
#' }
#'
#' This enables failure attribution analysis: what proportion of system
#' failures are caused by each component?
#'
#' @examples
#'
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_exponential(0.2),
#'     dfr_exponential(0.3)
#' ))
#'
#' set.seed(42)
#' mat <- sample_components(sys, n = 1000)
#' dim(mat)  # 1000 x 3
#'
#' # System lifetimes
#' t_sys <- apply(mat, 1, min)
#'
#' # Which component caused each failure?
#' failing <- apply(mat, 1, which.min)
#' table(failing) / 1000
#' # Proportions ~= c(1/6, 2/6, 3/6) for rates (0.1, 0.2, 0.3)
#'
#' @seealso
#' \code{\link[algebraic.dist]{sampler}} for system-level sampling,
#' \code{\link{component}} to extract individual component objects,
#' \code{\link{dfr_dist_series}} for the constructor
#' @family system introspection
#' @export
sample_components <- function(x, n, ...) {
    UseMethod("sample_components")
}
