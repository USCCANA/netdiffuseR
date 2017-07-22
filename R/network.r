#' Convertion between graph classes
#'
#' @param graph An object of class \code{\link{diffnet}}
#' @param slices An integer vector indicating the slices to subset
#' @param ... Further arguments passed to \code{\link[networkDynamic:networkDynamic]{networkDynamic}}
#'
#' @name network
NULL

#' @rdname network
#' @return \code{diffnet_to_network} returns a list of length \code{length(slices)} in which
#' each element is a \code{\link[network:network]{network}} object corresponding a slice of the
#' \code{graph} (\code{diffnet} object). The attributes list will include \code{toa} (time of
#' adoption).
#' @family Foreign
#' @export
diffnet_to_network <- function(graph, slices = 1:nslices(graph), ...) {
  if (!inherits(graph, "diffnet"))
    stop("-graph- should be a diffnet object.")

  graph <- graph[,,slices]
  n     <- nnodes(graph)

  structure(
    Map(function(g, a, time) {
      dimnames(g) <- list(rownames(graph), rownames(graph))

      network::network(
        x = as.matrix(g),
        vertex.attr = c(
          list(toa = graph$toa),
          unclass(a),
          unclass(graph$vertex.static.attrs)
        )
      )
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
#' data(medInnovationsDiffNet)
#' diffnet_to_networkDynamic(medInnovationsDiffNet)
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
#' @param toavar Character scalar. Name of the vertex attribute that holds the times of adoption.
#' @export
#' @details
#' By default, \code{networkDynamic_to_diffnet} uses the first slice as reference for
#' vertex attributes and times of adoption.
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
  as_diffnet(graph = net, toa = toa, vertex.dyn.attrs = vattrs,
             t0 = min(min(toa, na.rm=TRUE), pers[1]),
             t1 = max(max(toa, na.rm=TRUE), pers[length(pers)])
             )
}
