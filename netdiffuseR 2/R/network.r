#' Coercion between \code{diffnet}, \code{network} and \code{networkDynamic}
#'
#' @param graph An object of class \code{\link{diffnet}}
#' @param slices An integer vector indicating the slices to subset
#' @param toavar Character scalar. Name of the vertex attribute that holds the times of adoption.
#' @param ... Further arguments passed to \code{\link[networkDynamic:networkDynamic]{networkDynamic}}
#'
#' @section Caveats:
#' Since \code{diffnet} does not support edges attributes, these will be lost when
#' converting from \code{network}-type objects. The same applies to \code{network}
#' attributes.
#'
#' @name network
#' @aliases networkDynamic
NULL

#' @rdname network
#' @return \code{diffnet_to_network} returns a list of length \code{length(slices)} in which
#' each element is a \code{\link[network:network]{network}} object corresponding a slice of the
#' \code{graph} (\code{diffnet} object). The attributes list will include \code{toa} (time of
#' adoption).
#' @family Foreign
#' @export
#'
#' @examples
#' # Cohersing a diffnet to a list of networks ---------------------------------
#' set.seed(1)
#' ans <- diffnet_to_network(rdiffnet(20, 2))
#' ans
#'
#' # and back
#' network_to_diffnet(graph.list = ans, toavar="toa")
#'
#' # If it was static, we can use -graph- instead
#' network_to_diffnet(ans[[1]], toavar="toa")
#'
diffnet_to_network <- function(graph, slices = 1:nslices(graph), ...) {
  if (!inherits(graph, "diffnet"))
    stop("-graph- should be a diffnet object.")

  graph <- graph[,,slices]
  n     <- nnodes(graph)

  structure(
    Map(function(g, a, time) {
      dimnames(g) <- list(rownames(graph), rownames(graph))

      # Creating the network object
      ans <- network::network(
        x = as.matrix(g),
        vertex.attr = c(
          list(toa = graph$toa),
          unclass(a),
          unclass(graph$vertex.static.attrs)
        ),
        loops = graph$meta$self
      )

      # Setting attributes
      network::set.network.attribute(ans, "name", graph$meta$name)
      network::set.network.attribute(ans, "behavior", graph$meta$behavior)

      ans
    }, g = graph$graph, a = graph$vertex.dyn.attrs)
    ,
    names = dimnames(graph)[[3]]
  )
}

#' @rdname network
#' @details \code{diffnet_to_networkDynamic} calls \code{diffnet_to_network} and
#' uses the output to call \code{networkDynamic}, passing the resulting list of
#' \code{network} objects as \code{network.list} (see \code{\link[networkDynamic:networkDynamic]{networkDynamic}}).
#'
#' By default, \code{diffnet_to_networkDynamic} passes \code{net.obs.period} as
#' \preformatted{
#'   net.obs.period = list(
#'     observations = list(range(graph$meta$pers)),
#'     mode="discrete",
#'     time.increment = 1,
#'     time.unit = "step"
#'   )
#' }
#'
#' @export
#' @param diffnet2net.args List of arguments passed to \code{diffnet_to_network}.
#' @param netdyn.args List of arguments passed to \code{\link[networkDynamic:networkDynamic]{networkDynamic}}
#' @return An object of class \code{networkDynamic}.
#' @examples
#' # A random diffusion network ------------------------------------------------
#' set.seed(87)
#' dn  <- rdiffnet(50, 4)
#' ans <- diffnet_to_networkDynamic(dn)
#'
#' # and back
#' networkDynamic_to_diffnet(ans, toavar = "toa")
#'
diffnet_to_networkDynamic <- function(
  graph,
  slices = 1:nslices(graph),
  diffnet2net.args = list(),
  netdyn.args = list()
) {

  # Checking class and cutting
  if (!inherits(graph, "diffnet"))
    stop("-graph- should be a diffnet object.")

  graph <- graph[,,slices]

  # Checking arguments in net.obs.period
  if ("net.obs.period" %in% names(netdyn.args)) {
    net.obs.period <- netdyn.args[["net.obs.period"]]
    netdyn.args[["net.obs.period"]] <- NULL
  }
  else
    net.obs.period <- list(
      observations = list(range(graph$meta$pers)),
      mode="discrete",
      time.increment = 1,
      time.unit = "step"
    )

  # Generating list of networks
  ans <- do.call(diffnet_to_network, c(list(graph=graph), diffnet2net.args))

  # Calling networkDynamic
  do.call(
    networkDynamic::networkDynamic,
    c(
      list(
        network.list   = ans,
        onsets         = graph$meta$pers,
        termini        = graph$meta$pers,
        net.obs.period = net.obs.period
      ),
      netdyn.args
    )
  )


}

#' @rdname network
#' @export
#' @details
#' By default, \code{networkDynamic_to_diffnet} uses the first slice as reference for
#' vertex attributes and times of adoption.
#'
networkDynamic_to_diffnet <- function(graph, toavar) {

  # Getting timeframe
  pers <- range(unlist(network::get.network.attribute(graph, "net.obs.period")[["observations"]]))
  pers <- pers[1]:pers[2]

  # Extracting according to time
  ans <- lapply(pers, function(i) {
    networkDynamic::network.extract(
      x = graph,
      at = i
      )
  })

  # Collapsing
  ans <- lapply(ans, networkDynamic::network.collapse)

  # Extracting adjacency matrix
  vnames <- network::get.vertex.attribute(ans[[1]], "vertex.names")
  net    <- lapply(ans, function(x) {
    ans <- methods::as(as.matrix(x), "dgCMatrix")
    dimnames(ans) <- list(vnames,vnames)
    ans
    })
  names(net) <- pers

  # Extracting toa
  toa <- network::get.vertex.attribute(ans[[1]], attrname = toavar)

  # Extracting attributes
  vattrs <- setdiff(network::list.vertex.attributes(ans[[1]]), c(toavar, "vertex.names"))
  vattrs <- lapply(ans, function(x) {
      ans <- lapply(vattrs, function(y) network::get.vertex.attribute(x, y))
      names(ans) <- vattrs
      as.data.frame(ans)
    }
  )

  # Coercing into a diffnet object
  new_diffnet(
    graph            = net,
    toa              = toa,
    vertex.dyn.attrs = vattrs,
    t0               = min(min(toa, na.rm=TRUE), pers[1]),
    t1               = max(max(toa, na.rm=TRUE), pers[length(pers)]),
    name             = network::get.network.attribute(ans[[1]], "name"),
    behavior         = network::get.network.attribute(ans[[1]], "behavior"),
    self             = network::has.loops(ans[[1]])
    )
}

#' @rdname network
#' @param graph.list A list of \code{network} objects.
#' @param t0 Integer scalar. Passed to \code{\link[=diffnet-class]{new_diffnet}}.
#' @param t1 Integer scalar. Passed to \code{\link[=diffnet-class]{new_diffnet}}.
#' @export
#' @details
#' By default, \code{network_to_diffnet} uses the first element of \code{graph}
#' (a list) as reference for vertex attributes and times of adoption.
network_to_diffnet <- function(graph = NULL, graph.list = NULL, toavar,
                               t0 = NULL, t1 = NULL) {

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
    if (!inherits(graph, "network"))
      stop("-graph- should be of class -network-.")
  } else {
    if (!inherits(graph.list, "list"))
      stop("-graph.list- should be a list")
    else {

      # All must be of class igraph
      test <- which(!sapply(graph.list, inherits, what="network"))
      if (length(test))
        stop("The following elements of -graph.list- are not of class -network-",
             paste(test, collapse=", "), ".")

      # All must have the same attributes
      test <- network::list.vertex.attributes(graph.list[[1]])
      test <- which(sapply(graph.list, function(x) {
        !all(network::list.vertex.attributes(x) %in% test)
      }))
      if (length(test))
        stop("All -network- objects in -graph.list- must have the same",
             " vertex attributes. The following differ from the first element:",
             paste(test, collapse=", "), ".")
    }

    islist <- TRUE
  }

  # # Getting timeframe
  # pers <- if (!islist) range(unlist(network::get.vertex.attribute(graph, toavar)), na.rm=TRUE)
  # else range(unlist(network::get.vertex.attribute(graph.list[[1]], toavar)), na.rm=TRUE)
  #
  # pers <- pers[1]:pers[2]

  # Extracting adjacency matrix
  vnames <- if (!islist) network::get.vertex.attribute(graph, "vertex.names")
  else network::get.vertex.attribute(graph.list[[1]], "vertex.names")

  net <- if (!islist) methods::as(as.matrix(graph), "dgCMatrix")
  else lapply(graph.list, function(x) {
    ans <- methods::as(as.matrix(x), "dgCMatrix")
    dimnames(ans) <- list(vnames,vnames)
    ans
  })

  # Extracting toa
  toa <- if (!islist) network::get.vertex.attribute(graph, attrname = toavar)
  else network::get.vertex.attribute(graph.list[[1]], attrname = toavar)

  # Extracting attributes
  vattrs <- if (!islist) setdiff(network::list.vertex.attributes(graph), c(toavar, "vertex.names"))
  else setdiff(network::list.vertex.attributes(graph.list[[1]]), c(toavar, "vertex.names"))

  vertex.static.attrs <- NULL
  vertex.dyn.attrs    <- NULL

  if (!islist) {
    # Static attributes
    vertex.static.attrs <- lapply(
      vattrs, function(x) {
        network::get.vertex.attribute(graph, x)
        }
      )

    names(vertex.static.attrs) <- vattrs
    vertex.static.attrs <- as.data.frame(vertex.static.attrs)
  } else {

    # Dynamic attributes
    vertex.dyn.attrs <- lapply(
      graph.list, function(g) {

        # Getting attributes
        ans <- lapply(vattrs, function(x) {
          network::get.vertex.attribute(g, x)
        })

        # Cohersing as data.frame
        names(ans) <- vattrs
        as.data.frame(ans)
      }
    )
  }

  # Coercing into a diffnet object
  if (!islist) graph.list <- list(graph)

  new_diffnet(
    graph            = net,
    toa              = toa,
    vertex.dyn.attrs = vertex.dyn.attrs,
    vertex.static.attrs = vertex.static.attrs,
    t0               = ifelse(length(t0), t0, min(toa, na.rm=TRUE)),
    t1               = ifelse(length(t1), t1, max(toa, na.rm=TRUE)),
    name             = network::get.network.attribute(graph.list[[1]], "name"),
    behavior         = network::get.network.attribute(graph.list[[1]], "behavior"),
    self             = network::has.loops(graph.list[[1]])
  )
}
