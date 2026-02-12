#' Get the number of components in a system
#'
#' @param x A system object.
#' @param ... Additional arguments.
#' @return Integer, the number of components.
#' @export
ncomponents <- function(x, ...) {
    UseMethod("ncomponents")
}

#' Extract a component from a system
#'
#' @param x A system object.
#' @param j Component index.
#' @param ... Additional arguments.
#' @return A distribution object for component j.
#' @export
component <- function(x, j, ...) {
    UseMethod("component")
}

#' Get the parameter layout for a system
#'
#' @param x A system object.
#' @param ... Additional arguments.
#' @return A list mapping component indices to global parameter indices.
#' @export
param_layout <- function(x, ...) {
    UseMethod("param_layout")
}

#' Get the hazard function for a specific component
#'
#' @param x A system object.
#' @param j Component index.
#' @param ... Additional arguments.
#' @return A closure computing component j's hazard.
#' @export
component_hazard <- function(x, j, ...) {
    UseMethod("component_hazard")
}

#' Sample component lifetimes from a system
#'
#' Generates an n x m matrix where column j contains independent samples
#' from component j's lifetime distribution.
#'
#' @param x A system object.
#' @param n Number of samples.
#' @param par Optional parameter vector override.
#' @param ... Additional arguments.
#' @return An n x m numeric matrix of component lifetimes.
#' @export
sample_components <- function(x, n, ...) {
    UseMethod("sample_components")
}
