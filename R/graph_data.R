#' Network data formats
#'
#' List of accepted graph formats
#'
#' @name netdiffuseR-graphs
#' @details The \pkg{netdiffuseR} package can handle different types of graph
#' objects. Two general classes are defined accross the package's functions:
#' static graphs, and dynamic graphs.
#' \itemize{
#'  \item{In the case of \strong{static graphs}, these are represented as adjacency
#'  matrices of size \eqn{n\times n}{n * n} and can be either \code{\link{matrix}}
#'  (dense matrices) or \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'  (sparse matrix from the \pkg{\link[Matrix:Matrix]{Matrix}} package). While
#'  most of the package functions are defined for both classes, the default output
#'  graph is sparse, i.e. \code{dgCMatrix}.}
#'  \item{With respect to \strong{dynamic graphs}, these are represented by either
#'  a \code{\link{diffnet}} object, an \code{\link{array}} of size \eqn{n\times n \times T}{n * n * T}, or a list of size \eqn{T}
#'  with sparse matrices (class \code{dgCMatrix}) of size \eqn{n\times n}{n * n}.
#'  Just like the static graph case, while most of the functions accept both
#'  graph types, the default output is \code{dgCMatrix}.}
#' }
#' @section diffnet objects:
#'  In the case of \code{diffnet}-class objects, the following arguments can be omitted
#'  when calling fuictions suitable for graph objects:
#'  \itemize{
#'    \item{\code{toa}: Time of Adoption vector}
#'    \item{\code{adopt}: Adoption Matrix}
#'    \item{\code{cumadopt}: Cumulative Adoption Matrix}
#'    \item{\code{undirected}: Whether the graph is directed or not}
#'  }
#'
#' @section Objects' names:
#' When possible, \pkg{netdiffuseR} will try to reuse graphs dimensional names,
#' this is, \code{\link{rownames}}, \code{\link{colnames}}, \code{\link{dimnames}}
#' and \code{\link{names}} (in the case of dynamic graphs as lists). Otherwise,
#' when no names are provided, these will be created from scratch.
#' @include imports.R
#' @author George G. Vega Yon
NULL

stopifnot_graph <- function(x)
  stop("No method for graph of class -",class(x),"-. Please refer to the manual 'netdiffuseR-graphs'.")

#' Analyze an R object to identify the class of graph (if any)
#' @param graph Any class of accepted graph format (see \code{\link{netdiffuseR-graphs}}).
#' @details This function analyzes an R object and tries to classify it among the
#' accepted classes in \pkg{netdiffuseR}. If the object fails to fall in one of
#' the types of graphs the function returns with an error indicating what (and
#' when possible, where) the problem lies.
#'
#' The function was designed to be used with \code{\link{as_diffnet}}.
#' @seealso \code{\link{as_diffnet}}, \code{\link{netdiffuseR-graphs}}
#' @return Whe the object fits any of the accepted graph formats, a list of attributes including
#' \item{type}{Character scalar. Whether is a static or a dynamic graph}
#' \item{class}{Character scalar. The class of the original object}
#' \item{ids}{Character vector. Labels of the vertices}
#' \item{pers}{Integer vector. Labels of the time periods}
#' \item{nper}{Integer scalar. Number of time periods}
#' \item{n}{Integer scalar. Number of vertices in the graph}
#' Otherwise returns with error.
#' @author George G. Vega Yon
#' @export
classify_graph <- function(graph) {

  # Diffnet object
  if (inherits(graph, "diffnet")) {
    return(classify_graph(graph$graph))
  } else if (inherits(graph, "matrix") || inherits(graph, "dgCMatrix")) { # Static graphs
    # Step 0: Should have length
    d <- dim(graph)
    if (!d[1])
      stop("Nothing to do. Empty matrix.")

    # Step 1: Should be square
    if (d[1] != d[2])
      stop("-graph- must be a square matrix\n\tdim(graph) = c(",
           paste0(d, collapse=","),").")

    # Step 3: Should be numeric
    m <- mode(graph)
    if (!inherits(graph, "dgCMatrix") && !(m %in% c("numeric", "integer")))
      stop("-graph- should be either numeric or integer.\n\tmode(graph) =  \"",
           m, "\".")

    # Step 4: Dimension names
    ids <- rownames(graph)
    if (!length(ids)) ids <- 1:d[1]

    return(invisible(list(
      type="static",
      class="matrix",
      ids=ids,
      pers=1,
      nper=1,
      n=d[1]
    )))
  }
  # Dynamic graphs (list) ------------------------------------------------------
  else if (inherits(graph, "list")) {
    # Step 0: Should have length!
    t <- length(graph)
    if (t < 2)
      stop("-graph- must be at least of length 2.")

    # Step 1: All should be of class -dgCMatrix-
    c <- sapply(graph, inherits, "dgCMatrix")
    if (!all(c))
      stop("The following elements are not of class -dgCMatrix-:\n\t",
           paste0(which(!c), collapse=", "),".")

    # Step 2.1: All must be square matrices
    d <- lapply(graph, dim)
    s <- sapply(d, function(x) x[1] == x[2])

    # Step 2.2: It must have some people!
    if (!d[[1]][1])
      stop("Nothing to do. Empty graph.")

    if (!all(s))
      stop("The following adjmat are not square:\n\t",
           paste0(which(!s), collapse=", "),".")

    # Step 3: All must have the same dimension
    e <- unlist(d, TRUE) == d[[1]][1]
    if (!all(e))
      stop("The dimensions of all slices must be equal. ",
           "The following elements don't coincide with the first slice:\n\t",
           paste0(which(!e), collapse=", "),".")

    # Step 4.1: Individual's ids
    ids <- rownames(graph[[1]])
    if (!length(ids)) ids <- 1:d[[1]][1]

    # Step 4.2 Time ids
    suppressWarnings(pers <- as.integer(names(graph)))
    if (!length(pers)) pers <- 1:t
    else {
      # Step 4.2.1: Must be coercible into integer
      if (any(is.na(pers))) stop("names(graph) should be either numeric or integer.")

      # Step 4.2.1: Must keep uniqueness
      if (length(unique(pers)) != t) stop("When coercing names(graph) into integer,",
                                       "some slices acquired the same name.")
    }

    return(invisible(list(
      type="dynamic",
      class="list",
      ids=ids,
      pers=pers,
      nper=t,
      n=d[[1]][1])
    ))

  }
  # Dynamic graphs (array) -----------------------------------------------------
  else if (inherits(graph, "array")) {
    # Step 0: it should have length!
    d <- dim(graph)
    if (d[3] < 2)
      stop("-graph- must be at least of length 2.")

    # Step 1: there must be some people
    if (!d[1])
      stop("Nothing to do. Empty matrix.")

    # Step 2: It must be square
    if (d[1] != d[2])
      stop("Each adjmat in -graph- must be a square matrix\n\tdim(graph) = c(",
           paste0(d, collapse=","),").")

    # Step 3: Should be numeric
    m <- mode(graph)
    if (!(m %in% c("numeric", "integer")))
      stop("-graph- should be either numeric or integer.\n\tmode(graph) =  \"",
           m, "\".")

    # Step 4: Dimension names
    ids <- rownames(graph)
    if (!length(ids)) ids <- 1:d[1]

    pers <- as.numeric(dimnames(graph)[[3]])
    if (!length(pers)) pers <- 1:d[3]
    else {
      # Step 4.2.1: Must be coercible into integer
      suppressWarnings(alters <- as.integer(floor(pers)))
      if (any(is.na(alters))) stop("names(graph) should be either numeric or integer.")

      # Step 4.2.1: Must keep uniqueness
      if (length(unique(alters)) != length(pers))
        stop("When coercing names(graph) into integer,",
             "some slices acquired the same name.")
      pers <- alters
    }

    return(invisible(list(
      type="dynamic",
      class="array",
      ids=ids,
      pers=pers,
      nper=d[3],
      n=d[1])
    ))
  }

  # Other case (ERROR) ---------------------------------------------------------
  stop("Not an object allowed in netdiffuseR. It must be either:\n\t",
       "matrix, dgCMatrix, list or array.\n", "Please refer to ?\"netdiffuseR-graphs\" ")
}

# Auxiliar function to check if there's any attribute of undirectedness
checkingUndirected <- function(graph, warn=TRUE, default=getOption("diffnet.undirected")) {

  # Ifendifying the class of graph
  if (inherits(graph, "diffnet")) undirected <- graph$meta$undirected
  else undirected <- attr(graph, "undirected")

  if (warn)
    if (length(undirected) && undirected != FALSE)
      warning("The entered -graph- will now be directed.")

  if (!length(undirected)) undirected <- default

  invisible(undirected)

}
