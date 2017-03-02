#' Convertion between graph classes
#' @param graph Either a \code{\link{diffnet}} or \code{\link[igraph:igraph]{igraph}} graph object.
#' @param slices An integer vector indicating the slices to subset.
#' @return Either a list of \code{length(slices)} \code{igraph}
#' (\code{diffnet_to_igraph}), or a \code{diffnet} object (\code{igraph_to_diffnet})
#' objects.
#' @export
#' @family graph formats
#' @examples
#' # Reading the meddical innovation data into igraph --------------------------
#' x <- diffnet_to_igraph(medInnovationsDiffNet)
#'
#' # Fetching the times of adoption
#' igraph::vertex_attr(x[[1]], "toa")
#'
diffnet_to_igraph <- function(graph, slices=1:nslices(graph)) {

  if (!inherits(graph, "diffnet")) stopifnot_graph(graph)

  pers <- graph$meta$pers
  nper <- graph$meta$nper

  # Listing attributes
  static.attrs <- colnames(graph$vertex.static.attrs)
  dynamic.attrs <- colnames(graph$vertex.dyn.attrs)

  # Creating container
  out <- vector("list", length(slices))
  names(out) <- pers[slices]


  for (p in 1:length(out)) {
    # Index
    s <- slices[p]
    dimnames(graph$graph[[s]]) <- with(graph$meta, list(ids,ids))

    # Graph
    tempgraph <- igraph::graph_from_adjacency_matrix(
      graph[,,s,drop=TRUE][[1]],
      mode      = ifelse(graph$meta$undirected, "undirected", "directed"),
      weighted  = TRUE,
      diag      = graph$meta$self
      )

    # Computing positions
    tempgraph <- igraph::permute(tempgraph, match(igraph::V(tempgraph)$name,
                                                   graph$meta$ids))

    # Vertex Static Attributes
    for (k in static.attrs)
      tempgraph <-
      igraph::set_vertex_attr(graph=tempgraph, name=k, value=graph[[k]])

    # Vertex Dyn Attributes
    for (k in dynamic.attrs)
      tempgraph <-
      igraph::set_vertex_attr(graph=tempgraph, name=k, graph[[k]][[s]])

    # Time of adoption
    tempgraph <- set_vertex_attr(
      graph=tempgraph, name="toa", value=graph$toa
      )

    # Attributes
    tempgraph<-igraph::set_graph_attr(tempgraph, "name", graph$meta$name)
    tempgraph<-igraph::set_graph_attr(tempgraph, "behavior", graph$meta$behavior)

    # Storing the value
    out[[p]] <- tempgraph
  }
  out
}

#' @export
#' @rdname diffnet_to_igraph
#' @param toavar Character scalar. Name of the attribute that holds the times of adoption.
#' @param t0 Integer scalar. Passed to \code{\link{as_diffnet}}.
#' @param t1 Integer scalar. Passed to \code{\link{as_diffnet}}.
#' @param ... Further arguments passed to \code{\link{as_diffnet}}.
igraph_to_diffnet <- function(
  graph, toavar,
  t0 = NULL,
  t1 = NULL, ...) {

  if (missing(toavar))
    stop("Please provide the name of the -toa- var.")

  if (igraph::is_igraph(graph)) {
    # Getting attributes
    toa <- igraph::vertex_attr(graph, toavar)

    # Checking toa
    if (!length(t0)) t0 <- min(toa, na.rm = TRUE)
    if (!length(t1)) t1 <- max(toa, na.rm = TRUE)

    mat <- igraph::as_adj(graph, attr="weight")
    mat <- lapply(seq_len(t1-t0+1), function(x) mat)

    vertex.static.attrs <- as.data.frame(
      igraph::vertex.attributes(graph))

    if (length(vertex.static.attrs)) {
      test <- which(!(colnames(vertex.static.attrs) %in% c(toavar, "name")))
      vertex.static.attrs <- vertex.static.attrs[,test,drop=FALSE]
    }

    return(as_diffnet(mat, toa, t0, t1,
               vertex.static.attrs = vertex.static.attrs,
               name = igraph::graph_attr(graph, "name"),
               behavior = igraph::graph_attr(graph, "behavior"),
               ...))
  }

  stopifnot_graph(graph)
}
