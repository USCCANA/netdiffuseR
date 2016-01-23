#' Retrieve alter's attributes
#'
#' For a given set of vertices \eqn{E}, retrieves each vertex's ego network alter's
#' attributes.
#'
#' @param graph Any kind of graph
#' @param attrs Numeric matrix with \eqn{n} rows.
#' @param E Integer vector. Set of vertices from which the attributes will be retrieved.
#' @param outer Logical scalar. When TRUE builds ego nets using outer edges.
#' @param self Logical scalar. When FALSE ignores ego's own attributes.
#' @param valued Logical scalar. When TRUE stores the value of the edge, otherwise includes a one.
#' @details FUNCTION ON DEVELOPMENT, CURRENTLY WORKING ONLY FOR DIFFNET OBJECTS
#' @examples
#' # Creating a random graph
#' set.seed(1001)
#' diffnet <- rdiffnet(150, 20, seed.graph="small-world")
#'
#' # Adding attributes
#' indeg <- dgr(diffnet, cmode="indegree")
#' diffnet.attrs(diffnet, "vertex", "dyn") <-
#'  lapply(1:20, function(x) cbind(indeg=indeg[,x]))
#'
#' # Retrieving egonet's attributes (vertices 1 and 20)
#' egonet_attrs(diffnet, E=c(1,20))
#'
#' @export
egonet_attrs <- function(
  graph, attrs, E=NULL,
  outer = TRUE,
  self = getOption("diffnet.self"),
  valued = getOption("diffnet.valued")
) {

  # Checking if no dim has been specified
  if (!inherits(graph, "diffnet"))
    if (missing(attrs))
      stop("If -graph- is not of class 'diffnet', -attrs- must be specified.")

  # Verifying dimensions
  if (inherits(graph, "diffnet")) {
    if (missing(attrs)) attrs <- diffnet.attrs(graph)
    else {
      if (!inherit(attrs, "list")) stop("-attrs- must be a list.")

      if (length(attrs) != graph$meta$nper)
        stop("-attrs-, ",length(attrs),
             " elements, must have as many elements as time periods ",
             graph$meta$nper,".")
    }
  } else stop_ifnotgraph(graph)

  switch(
    class(graph),
    diffnet = egonet_attrs.list(
      graph$graph, attrs, if (!length(E)) graph$meta$ids else E, outer, self, valued),
    list = egonet_attrs.list(graph, attrs, E, outer, self, valued)
  )

}

# For lists
egonet_attrs.list <- function(graph, attrs, E, outer, self, valued) {
  nper <- length(graph)
  if (nper != length(attrs))
    stop("-graph- and -attrs- must have the same length")

  out <- lapply(1:nper, function(x) {
    egonet_attrs_cpp(graph[[x]], as.integer(E), attrs[[x]], outer, self, valued)
  })

  # Adding names to E
  out <- lapply(out, `names<-`, E)

  # Naming
  tn <- names(graph)
  if (!length(tn)) tn <- 1:length(graph)
  names(out) <- tn

  out
}
