#' Get the number of components in a system
#'
#' Returns the number of components \eqn{m} in a series system, corresponding
#' to the number of terms in \eqn{h_{sys}(t) = \sum_{j=1}^{m} h_j(t)}.
#'
#' @param x A system object (e.g., \code{\link{dfr_dist_series}}).
#' @param ... Additional arguments passed to methods.
#' @return Integer, the number of components.
#'
#' @details
#' For a \code{\link{dfr_dist_series}} object created from a list of
#' \code{m} components, this simply returns \code{m}. This is useful for
#' programmatically iterating over components, e.g., for plotting hazard
#' decompositions or computing failure attribution.
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_exponential(0.1),
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_gompertz(a = 0.01, b = 0.05)
#' ))
#' ncomponents(sys)  # 3
#' }
#'
#' @seealso
#' \code{\link{component}} to extract individual components,
#' \code{\link{dfr_dist_series}} for the constructor
#' @family system introspection
#' @export
ncomponents <- function(x, ...) {
    UseMethod("ncomponents")
}

#' Extract a component from a system
#'
#' Extracts component \code{j} from a series system as a standalone
#' \code{\link[flexhaz]{dfr_dist}} object, with its parameters set to the
#' current values from the system's parameter vector (via the layout).
#'
#' @param x A system object (e.g., \code{\link{dfr_dist_series}}).
#' @param j Component index (integer, \code{1 <= j <= ncomponents(x)}).
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[flexhaz]{dfr_dist}} object for component \code{j}.
#'
#' @details
#' The returned component object is a copy of the original component with its
#' \code{par} field updated to reflect the current system-level parameter
#' vector. This means you can evaluate the extracted component's hazard,
#' survival, etc. directly:
#'
#' \preformatted{
#' comp1 <- component(sys, 1)
#' h1 <- hazard(comp1)
#' h1(10)  # evaluates using parameters from the system
#' }
#'
#' Changes to the extracted component do \emph{not} propagate back to the
#' original series system.
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_exponential(0.05)
#' ))
#'
#' # Extract the Weibull component
#' wb <- component(sys, 1)
#' params(wb)  # c(2, 100)
#'
#' # Evaluate its hazard independently
#' h_wb <- hazard(wb)
#' h_wb(50)
#' }
#'
#' @seealso
#' \code{\link{ncomponents}} for the component count,
#' \code{\link{component_hazard}} for getting just the hazard closure,
#' \code{\link{param_layout}} for parameter index mapping,
#' \code{\link{dfr_dist_series}} for the constructor
#' @family system introspection
#' @export
component <- function(x, j, ...) {
    UseMethod("component")
}

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
#' \donttest{
#' library(flexhaz)
#'
#' sys <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_exponential(0.05),
#'     dfr_gompertz(a = 0.01, b = 0.1)
#' ))
#' param_layout(sys)
#' # list(1:2, 3, 4:5)
#' }
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
#' \donttest{
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
#' }
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
#' \donttest{
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
#' }
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
