#' Retrieve alter's attributes
#' @param graph Any kind of graph
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
  } else netdiffuseR:::stop_ifnotgraph(graph)

  switch(
    class(graph),
    diffnet = egonet_attrs.list(
      graph$graph, attrs, if (!length(E)) graph$meta$ids else E, outer, self, valued),
    list = egonet_attrs.list(graph, attrs, E, outer, self, valued)
  )

}

#' For lists
egonet_attrs.list <- function(graph, attrs, E, outer, self, valued) {
  nper <- length(graph)
  if (nper != length(attrs))
    stop("-graph- and -attrs- must have the same length")

  lapply(1:nper, function(x) {
    netdiffuseR:::egonet_attrs_cpp(graph[[x]], as.integer(E), attrs[[x]], outer, self, valued)
  })
}
