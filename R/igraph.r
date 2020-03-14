#' Coercion between graph classes
#'
#' @param graph Either a \code{\link{diffnet}} or \code{\link[igraph:igraph]{igraph}} graph object.
#' @param slices An integer vector indicating the slices to subset.
#' @return Either a list of \code{length(slices)} \code{igraph}
#' (\code{diffnet_to_igraph}), or a \code{diffnet} object (\code{igraph_to_diffnet})
#' objects.
#'
#' @family Foreign
#' @examples
#' # Reading the medical innovation data into igraph --------------------------
#' x <- diffnet_to_igraph(medInnovationsDiffNet[,,1:4])
#'
#' # Fetching the times of adoption
#' igraph::vertex_attr(x[[1]], "toa")
#' @name igraph
NULL

#' @rdname igraph
#' @export
diffnet_to_igraph <- function(graph, slices=1:nslices(graph)) {

  if (!inherits(graph, "diffnet"))
    stopifnot_graph(graph)

  # Filtering
  graph <- graph[,,slices]

  # Listing attributes
  static.attrs <- colnames(graph$vertex.static.attrs)
  dynamic.attrs <- colnames(graph$vertex.dyn.attrs[[1]])

  # Creating container
  out <- vector("list", nslices(graph))
  names(out) <- dimnames(graph)[[3]]

  for (p in 1:length(out)) {
    # Index
    tempgraph <- graph$graph[[p]]
    dimnames(tempgraph) <- with(graph$meta, list(ids,ids))

    # Graph
    tempgraph <- igraph::graph_from_adjacency_matrix(
      adjmatrix = tempgraph,
      mode      = ifelse(graph$meta$undirected, "undirected", "directed"),
      weighted  = TRUE,
      diag      = graph$meta$self
      )

    # Vertex Static Attributes
    for (k in static.attrs)
      tempgraph <-
      igraph::set_vertex_attr(graph=tempgraph, name=k, value=graph[[k]])

    # Vertex Dyn Attributes
    for (k in dynamic.attrs)
      tempgraph <-
      igraph::set_vertex_attr(graph=tempgraph, name=k, value=graph[[k]][[p]])

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
#' @rdname igraph
#' @param toavar Character scalar. Name of the attribute that holds the times of adoption.
#' @param t0 Integer scalar. Passed to \code{\link[=diffnet-class]{new_diffnet}}.
#' @param t1 Integer scalar. Passed to \code{\link[=diffnet-class]{new_diffnet}}.
#' @param graph.list A list of \code{igraph} objects.
#' @param ... Further arguments passed to \code{\link{as_diffnet}}.
igraph_to_diffnet <- function(
  graph = NULL,
  graph.list = NULL,
  toavar,
  t0 = NULL,
  t1 = NULL, ...) {

  # At least one
  if (!length(graph.list) & !length(graph))
    stop("Either -graph.list- or -graph- should be provided.")
  else if (length(graph.list) & length(graph))
    stop("Only one of -graph.list- or -graph- should be provided.")

  # Toavar cannot be missing
  if (missing(toavar))
    stop("Please provide the name of the -toa- var.")

  # Checking list
  islist <- FALSE
  if (length(graph)) {
    if (!inherits(graph, "igraph"))
      stop("-graph- should be of class -igraph-.")
  } else {
    if (!inherits(graph.list, "list"))
      stop("-graph.list- should be a list")
    else {

      # All must be of class igraph
      test <- which(!sapply(graph.list, igraph::is_igraph))
      if (length(test))
        stop("The following elements of -graph.list- are not of class -igraph-",
             paste(test, collapse=", "), ".")

      # All must have the same attributes
      test <- igraph::list.vertex.attributes(graph.list[[1]])
      test <- which(sapply(graph.list, function(x) {
        !all(igraph::list.vertex.attributes(x) %in% test)
      }))
      if (length(test))
        stop("All -igraph- objects in -graph.list- must have the same",
             " vertex attributes. The following differ from the first element:",
             paste(test, collapse=", "), ".")
    }

    islist <- TRUE
  }


  # Getting attributes
  toa <- if (!islist) igraph::vertex_attr(graph, toavar)
  else igraph::vertex_attr(graph.list[[1]], toavar)

  # Checking toa
  if (!length(t0)) t0 <- min(toa, na.rm = TRUE)
  if (!length(t1)) t1 <- max(toa, na.rm = TRUE)

  mat <- if (!islist) {
    wattr <- igraph::list.edge.attributes(graph)
    wattr <- if ("weights" %in% wattr) "weights" else NULL
    igraph::as_adj(graph, attr = wattr, sparse = TRUE)
  } else lapply(graph.list, function(g) {
    wattr <- igraph::list.edge.attributes(g)
    wattr <- if ("weights" %in% wattr) "weights" else NULL
    igraph::as_adj(g, attr=wattr, sparse=TRUE)
    })

  # Adjusting sizes
  if (!islist)
    mat <- lapply(seq_len(t1-t0+1), function(x) mat)

  # Getting the attributes
  vertex.static.attrs <- NULL
  vertex.dyn.attrs    <- NULL
  if (!islist) {
    vertex.static.attrs <- as.data.frame(
      igraph::vertex.attributes(graph))
  } else {
    vertex.dyn.attrs <- lapply(graph.list, function(x) {
      as.data.frame(igraph::vertex.attributes(x))
    })
  }

  # Removing toa
  if (length(vertex.static.attrs)) {
    test <- which(!(colnames(vertex.static.attrs) %in% c(toavar, "name")))
    vertex.static.attrs <- vertex.static.attrs[,test,drop=FALSE]
  }

  if (length(vertex.dyn.attrs)) {
    vertex.dyn.attrs <- lapply(vertex.dyn.attrs, function(x) {
      test <- which(!(colnames(x) %in% c(toavar, "name")))
      x[,test,drop=FALSE]
    })

  }

  if (!islist)
    graph.list <- list(graph)

  return(new_diffnet(mat, toa, t0, t1,
             vertex.static.attrs = vertex.static.attrs,
             vertex.dyn.attrs    = vertex.dyn.attrs,
             name                = igraph::graph_attr(graph.list[[1]], "name"),
             behavior            = igraph::graph_attr(graph.list[[1]], "behavior"),
             self                = any(igraph::is.loop(graph.list[[1]])),
             undirected          = !igraph::is.directed(graph.list[[1]]),
             ...))

}


# This is for creating themes for igraph plotting ------------------------------

# The default themes
igraph_plotting_defaults <- list(
  vertex.frame.color  = "#31688EFF",
  vertex.color        = "#21908CFF",
  edge.color          = grDevices::gray(4/7, .5),
  edge.arrow.size     = .25,
  vertex.label        = NA,
  vertex.label.color  = "black",
  vertex.label.family = "sans",
  vertex.label.dist   = 1,
  edge.curved         = TRUE,
  rescale             = TRUE,
  add                 = TRUE
)

plotting_defaults <- c(
  color  = "steelblue",
  family = "sans",
  lwd    = 2,
  border = igraph_plotting_defaults$edge.color,
  igraph_plotting_defaults,
  arrow.color = grDevices::gray(3/7, .75),
  background = grDevices::gray(.9)
)

# This function changes defaults accordignly
# The function takes -obj_name- which is the name of the object that holds
# the igraph parameters, and modifies it in the parent frame directly (so
# no copying over the functions, edition of the environment itself).
set_igraph_plotting_defaults <- function(obj_name) {

  env <- parent.frame()

  # All igraph plots are added
  if (length(env[[obj_name]]$add) && !env[[obj_name]]$add)
    stop("The argument -add- cannot be changed to FALSE.")

  # Checking all the reminder arguments
  for (default in names(igraph_plotting_defaults)) {
    if (!length(env[[obj_name]][[default]]))
      env[[obj_name]][[default]] <- igraph_plotting_defaults[[default]]
  }
}

# This function sets the default values for plotting parameters looking at
# parameters available in igraph_plotting_defaults
set_plotting_defaults <- function(params) {
  env <- parent.frame()
  for (param in params) {
    if (!exists(param, envir = env))
      stop("No such parameter!")

    if (!length(env[[param]]))
      env[[param]] <- plotting_defaults[[param]]


  }

}
