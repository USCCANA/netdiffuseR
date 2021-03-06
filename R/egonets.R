#' Retrieve alter's attributes (network effects)
#'
#' For a given set of vertices V, retrieves each vertex's alter's
#' attributes. This function enables users to calculate exposure on variables
#' other than the attribute that is diffusing.  Further, it enables the
#' specification of alternative functions to use to characterize ego's
#' personal network including calculating the mean, maximum, minimum, median, or
#' sum of the alters' attributes. These measures may be static or dynamic over
#' the interval of diffusion and they may be binary or valued.
#'
#' @templateVar self TRUE
#' @templateVar valued TRUE
#' @template graph_template
#' @param attrs If \code{graph} is static, Numeric matrix with \eqn{n} rows, otherwise a list of numeric matrices with \eqn{n} rows.
#' @param V Integer vector. Set of vertices from which the attributes will be retrieved.
#' @param direction Character scalar. Either \code{"outgoing"}, \code{"incoming"}.
#' @param fun Function. Applied to each
#' @param as.df Logical scalar. When TRUE returns a data.frame instead of a list (see details).
#' @param ... Further arguments to be passed to \code{fun}.
#' @details
#'
#' By indexing inner/outer edges, this function retrieves ego network attributes
#' for all \eqn{v \in V}{v in V}, which by default is the complete set
#' of vertices in the graph.
#'
#' When \code{as.df=TRUE} the function returns a data.frame of size
#' \eqn{(|V|\times T)\times k}{(|V| * T)*k} where \eqn{T} is the number of time
#' periods and \eqn{k} is the number of columns generated by the function.
#'
#' The function can be used to create network effects as those in the \pkg{RSiena}
#' package. The difference here is that the definition of the statistic directly
#' relies on the user. For example, in the \pkg{RSiena} package, the dyadic covariate
#' effect \emph{37. covariate (centered) main effect (X)}
#'
#' \deqn{%
#' s_{i37}(x) = \sum_j x_{ij}(w_{ij}-\bar w)
#' }{%
#' s_i37(x) = sum(x[ij] * (w[ij] - mean(w)))
#' }
#'
#' Which, having a diffnet object with attributes named \code{x} and \code{w},
#' can be calculated as
#'
#' \preformatted{
#'     egonet_attrs(diffnet, as.df=TRUE, fun=function(dat) {
#'      sum(dat[, "x"]*(dat[, "w"] - mean(dat[, "w"])))
#'     })
#'     }
#'
#' Furthermore, we could use the \emph{median} centered instead, for example
#'
#' \preformatted{
#'     egonet_attrs(diffnet, as.df=TRUE, fun=function(dat) {
#'      sum(dat[, "x"]*(dat[, "w"] - median(dat[, "w"])))
#'     })
#'     }
#'
#' Where for each \eqn{i}, \code{dat} will be a matrix with as many rows
#' as individuals in his egonetwork. Such matrix holds the column names of the
#' attributes in the network.
#'
#' When \code{self = TRUE}, it will include ego's attributes, regardless
#' the network has loops or not.
#'
#' @return A list with ego alters's attributes. By default, if the graph is static, the
#' output is a list of length \code{length(V)} with matrices having the following
#' columns:
#'
#' \item{value}{Either the corresponding value of the tie.}
#' \item{id}{Alter's id}
#' \item{...}{Further attributes contained in \code{attrs}}
#'
#' On the other hand, if \code{graph} is dynamic, the output is list of length
#' \eqn{T} of lists of length \code{length(V)} with data frames having the following
#' columns:
#'
#' \item{value}{The corresponding value of the adjacency matrix.}
#' \item{id}{Alter's id}
#' \item{per}{Time id}
#' \item{...}{Further attributes contained in \code{attrs}}
#'
#' @examples
#' # Simple example with diffnet -----------------------------------------------
#' set.seed(1001)
#' diffnet <- rdiffnet(150, 5, seed.graph="small-world")
#'
#' # Adding attributes
#' indeg <- dgr(diffnet, cmode="indegree")
#' head(indeg)
#' diffnet[["indegree"]] <- indeg
#'
#' # Retrieving egonet's attributes (vertices 1 and 20)
#' egonet_attrs(diffnet, V=c(1,20))
#'
#' # Example with a static network ---------------------------------------------
#'
#' set.seed(1231)
#' n <- 20
#' net <- rgraph_ws(n = n, k = 4, p = .5)
#' someattr <- matrix(rnorm(n * 2), ncol= 2, dimnames = list(NULL, c("a", "b")))
#'
#' # Maximum of -a- in ego network
#' ans <- egonet_attrs(net, someattr, fun = function(x) max(x[,"a"]))
#' ans
#'
#' # checking it worked, taking a look at node 1, 2, and 3
#' max(someattr[which(net[1,] == 1),"a"]) == ans[1] # TRUE
#' max(someattr[which(net[2,] == 1),"a"]) == ans[2] # TRUE
#' max(someattr[which(net[3,] == 1),"a"]) == ans[3] # TRUE
#'
#'
#' @export
#' @include graph_data.r
#' @author George G. Vega Yon
#' @family data management functions
egonet_attrs <- function(
  graph, attrs, V=NULL,
  direction = "outgoing",
  fun = function(x) x,
  as.df = FALSE,
  self = getOption("diffnet.self"),
  valued = getOption("diffnet.valued"),
  ...
) {

  if (direction == "incoming") outer <- FALSE
  else if (direction == "outgoing") outer <- TRUE
  else stop("-direction- must be either 'incoming' or 'outgoing'")

  # Checking if no dim has been specified
  if (!inherits(graph, "diffnet"))
    if (missing(attrs))
      stop("If -graph- is not of class 'diffnet', -attrs- must be specified.")

  # Verifying dimensions
  if (inherits(graph, "diffnet")) {
    if (missing(attrs)) attrs <- diffnet.attrs(graph)
    else {
      if (!inherits(attrs, "list")) stop("-attrs- must be a list.")

      if (length(attrs) != graph$meta$nper)
        stop("-attrs-, ",length(attrs),
             " elements, must have as many elements as time periods ",
             graph$meta$nper,".")
    }
  } # else stopifnot_graph(graph)

  # Checking the set of vertices
  if (!length(V))
    V <- 1L:nnodes(graph)
  else {

    # Must be unique
    V <- unique(V)

    # V must be integer, and no NAs
    V <- as.integer(V)
    V <- V[complete.cases(V)]

    if (!length(V))
      stop("After removing incomplete cases, -V- is empty.")

    # Checking range
    test <- range(V)
    if (test[1] < 1L | test[2] > nnodes(graph))
      stop("Some values of -V- are out of range. Should be between ",
           1, " and ", nnodes(graph), ".")
  }

  cls <- class(graph)

  if ("diffnet" %in% cls) {
    egonet_attrs.list(
      graph  = graph$graph,
      attrs  = attrs,
      V      = V,
      outer  = outer,
      fun    = fun,
      as.df  = as.df,
      self   = self,
      valued = valued,
      ...
    )
  } else if ("list" %in% cls) {
    egonet_attrs.list(
      graph  = graph,
      attrs  = attrs,
      V      = V,
      outer  = outer,
      fun    = fun,
      as.df  = as.df,
      self   = self,
      valued = valued,
      ...
    )
  } else if ("matrix" %in% cls) {
    egonet_attrs.matrix(
      graph  = as_spmat(graph),
      attrs  = attrs,
      V      = V,
      outer  = outer,
      fun    = fun,
      as.df  = as.df,
      self   = self,
      valued = valued,
      ...
    )
  } else if ("dgCMatrix" %in% cls) {
    egonet_attrs.matrix(
      graph  = graph,
      attrs  = attrs,
      V      = V,
      outer  = outer,
      fun    = fun,
      as.df  = as.df,
      self   = self,
      valued = valued,
      ...
    )
  } else if ("array" %in% cls) {
    egonet_attrs.array(
      graph  = graph,
      attrs  = attrs,
      V      = V,
      outer  = outer,
      fun    = fun,
      as.df  = as.df,
      self   = self,
      valued = valued,
      ...
    )
  } else
    stopifnot_graph(graph)

}

egonet_attrs.matrix <- function(graph, attrs, V, outer, fun, as.df, self, valued, ...) {

  ids <- egonet_attrs_cpp(graph, V - 1L, outer, self, valued)
  sapply(lapply(ids, function(w) cbind(
    w, attrs[w[,"id"],,drop=FALSE]
  )), fun, ...)
}

egonet_attrs.array <- function(graph, attrs, V, outer, fun, as.df, self, valued, ...) {
  # Coercing into list
  dn <- dimnames(graph)[[3]]
  if (!length(dn)) dimnames(graph)[[3]] <- 1:dim(graph)[3]

  graph <- apply(graph, 3, methods::as, Class='dgCMatrix')
  egonet_attrs.list(graph, attrs, V, outer, fun, as.df, self, valued)

}

# For lists
egonet_attrs.list <- function(graph, attrs, V, outer, fun, as.df, self, valued, ...) {
  nper <- length(graph)
  if (nper != length(attrs))
    stop("-graph- and -attrs- must have the same length")

  # Period times
  tn <- names(graph)
  if (!length(tn)) tn <- 1:nper

  out <- lapply(1:nper, function(x) {
      ans <- egonet_attrs_cpp(graph[[x]], as.integer(V)-1L, outer, self, valued)
      lapply(ans, function(w) cbind(
        w, attrs[[x]][w[,"id"],,drop=FALSE]
        ))
  })

  # Adding names to V
  out <- lapply(out, `names<-`, V)

  # Naming
  names(out) <- tn

  # Applying the function
  out <- lapply(out, lapply, fun, ...)

  # In case of data.frame
  if (as.df) {
    out <- do.call(rbind, lapply(out, do.call, what=rbind))
    return(data.frame(id=rownames(out), out, row.names = NULL))
  }

  out
}
