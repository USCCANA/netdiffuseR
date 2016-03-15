#' Combine diffnet objects
#'
#' Combining \code{\link{diffnet}} objects that share time periods and attributes names, but
#' vertices ids (only valid for diffnet objects that have an empty intersection
#' between vertices ids)
#'
#' @param ... diffnet objects to be concatenated.
#' @param recursive Ignored.
#' @details The diffnet objects in \code{...} must fulfill the following conditions:
#' \enumerate{
#'  \item Have the same time range,
#'  \item have the same vertex attributes, and
#'  \item have an empty intersection of vertices ids,
#' }
#'
#' The meta data regarding \code{undirected}, \code{value}, and \code{multiple}
#' are set to \code{TRUE} if any of the concatenating diffnet objects has that
#' meta equal to \code{TRUE}.
#'
#' The resulting diffnet object's columns in the vertex attributes ordering (both
#' dynamic and static) will coincide with the first diffnet's ordering.
#'
#' @return A new \code{diffnet} object with as many vertices as the sum of each
#' concatenated diffnet objects' number of vertices.
#'
#'
#' @examples
#' # Calculate structural equivalence exposure by city -------------------------
#' data(medInnovationsDiffNet)
#'
#' # Subsetting diffnets
#' city1 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 1]
#' city2 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 2]
#' city3 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 3]
#' city4 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 4]
#'
#' # Computing exposure in each one
#' city1[["expo_se"]] <- exposure(city1, alt.graph="se", valued=TRUE)
#' city2[["expo_se"]] <- exposure(city2, alt.graph="se", valued=TRUE)
#' city3[["expo_se"]] <- exposure(city3, alt.graph="se", valued=TRUE)
#' city4[["expo_se"]] <- exposure(city4, alt.graph="se", valued=TRUE)
#'
#' # Concatenating all
#' diffnet <- c(city1, city2, city3, city4)
#' diffnet
#'
#'
#' @export
#' @family diffnet methods
c.diffnet <- function(..., recursive=FALSE) {
  diffnets <- list(...)

  # Checking class
  test <- which(!sapply(diffnets, inherits, what="diffnet"))
  if (length(test))
    stop("Some objects are not of class -diffnet-:\n\t",
         paste0(test, collapse="', "), ".")

  # Comparing arguments
  N   <- 0
  IDS <- NULL
  for (dni in 1:length(diffnets)) {
    for (dnj in dni:length(diffnets)) {

      # Can't compare me to myself!
      if (dni == dnj) next

      # Checking ids -----------------------------------------------------------
      test <- intersect(diffnets[[dni]]$meta$ids, diffnets[[dnj]]$meta$ids)
      if (length(test))
        stop("No pair of diffnets can share vertices. The following diffnets have ",
             "common vertex ids:\n\t", dni, " and ", dnj, " with the ids ",
             paste0(test, collapse=", "), ".")

      # Checking pers ----------------------------------------------------------
      trani <- range(diffnets[[dni]]$meta$pers)
      tranj <- range(diffnets[[dnj]]$meta$pers)
      test  <-  which(trani != tranj)
      if (length(test))
        stop("All diffnets have to have the same time range. The following pair of",
             " diffnets have different time range:\n\t",
             dni, " (", trani[1] ," to ", trani[2],") and ",
             dnj, " (", tranj[1] ," to ", tranj[2],").")

      # Checking static attributes ---------------------------------------------
      static.a.i <- colnames(diffnets[[dni]]$vertex.static.attrs)
      static.a.j <- colnames(diffnets[[dnj]]$vertex.static.attrs)

      # Number of attributes
      test <- length(static.a.i) != length(static.a.j)
      if (test)
        stop("The number of -vertex.static.attrs- differs between diffnets ",dni,
             " and ",dnj,".")

      # Attributes names
      if (length(static.a.i)) {
        static.a.i <- sort(static.a.i)
        static.a.j <- sort(static.a.j)
        test <- which(static.a.i != static.a.j)
        if (length(test))
          stop("All attributes names must match. The following diffnet objects",
               " have different static attributes colnames:\n\t", dni, " and ", dnj,".")

        # Ordening
        static.a.i <- colnames(diffnets[[dni]]$vertex.static.attrs)
        diffnets[[dni]]$vertex.static.attrs <-
          diffnets[[dni]]$vertex.static.attrs[,static.a.i, drop=FALSE]

        # Ordening
        diffnets[[dnj]]$vertex.static.attrs <-
          diffnets[[dnj]]$vertex.static.attrs[,static.a.i, drop=FALSE]
      }

      # Checking dynamic attributes --------------------------------------------
      dyn.a.i <- colnames(diffnets[[dni]]$vertex.dyn.attrs[[1]])
      dyn.a.j <- colnames(diffnets[[dnj]]$vertex.dyn.attrs[[1]])

      # Number of attributes
      test <- length(dyn.a.i) != length(dyn.a.j)
      if (test)
        stop("The number of -vertex.dyn.attrs- differs between diffnets ",dni,
             " and ",dnj,".")

      # Attributes names
      if (length(dyn.a.i)) {
        dyn.a.i <- sort(dyn.a.i)
        dyn.a.j <- sort(dyn.a.j)
        test <- which(dyn.a.i != dyn.a.j)
        if (length(test))
          stop("All attributes names must match. The following diffnet objects",
               " have different dynamic attributes colnames:\n\t", dni, " and ", dnj,".")

        # Ordening
        dyn.a.i <- colnames(diffnets[[dni]]$vertex.dyn.attrs[[1]])

        for (per in 1:diffnets[[dni]]$meta$nper)
          diffnets[[dni]]$vertex.dyn.attrs[[per]] <-
          diffnets[[dni]]$vertex.dyn.attrs[[per]][,dyn.a.i, drop=FALSE]

        # Ordening
        for (per in 1:diffnets[[dni]]$meta$nper)
          diffnets[[dnj]]$vertex.dyn.attrs[[per]] <-
          diffnets[[dnj]]$vertex.dyn.attrs[[per]][,dyn.a.i, drop=FALSE]
      }
    }

    # Incrementing the number of individuals in the net
    N   <- N + diffnets[[dni]]$meta$n
    IDS <- c(IDS, diffnets[[dni]]$meta$ids)
  }

  # Building big matrix --------------------------------------------------------
  N    <- as.integer(N)
  nper <- as.integer(diffnets[[1]]$meta$nper)

  A        <- methods::new("dgCMatrix", Dim=c(N,N), p=rep(0L,N + 1L),
                       Dimnames = list(IDS, IDS))
  A        <- lapply(1:nper, function(x) A)
  names(A) <- diffnets[[1]]$meta$pers

  i    <- 1
  for (dni in 1:length(diffnets)) {

    # Range
    n   <- diffnets[[dni]]$meta$n
    ran <- i:(i + n - 1)

    # Adding adjacency matrices
    for (p in 1:nper)
      A[[p]][ran, ran] <- diffnets[[dni]]$graph[[p]]

    # Incrementing
    i <- i + n
  }

  # Building attributes --------------------------------------------------------
  vertex.static.attrs <- do.call(
    rbind, lapply(diffnets, function(x) x$vertex.static.attrs))

  vertex.dyn.attrs <- vector("list", nper)
  for (i in 1:nper)
    vertex.dyn.attrs[[i]] <- do.call(
      rbind, lapply(diffnets, function(x) x$vertex.dyn.attrs[[i]]))
  names(vertex.dyn.attrs) <- diffnets[[1]]$meta$pers

  graph.attrs         <- diffnets[[1]]$graph.attrs

  # Building time of adoption and adopt mats -----------------------------------
  toa      <- do.call(c, lapply(diffnets, function(x) x$toa))
  adopt    <- do.call(rbind, lapply(diffnets, function(x) x$adopt))
  cumadopt <- do.call(rbind, lapply(diffnets, function(x) x$cumadopt))

  # Building meta
  meta <- list(
    type  = "dynamic",
    class = diffnets[[1]]$meta$class,
    ids   = IDS,
    pers  = diffnets[[1]]$meta$pers,
    nper  = nper,
    n     = N,
    # Others
    self       = any(sapply(diffnets, function(x) x$meta$self)),
    undirected = any(sapply(diffnets, function(x) x$meta$undirected)),
    multiple   = any(sapply(diffnets, function(x) x$meta$multiple))
  )

  graph <- structure(list(graph=A, toa=toa, adopt=adopt, cumadopt=cumadopt,
                vertex.static.attrs =vertex.static.attrs,
                vertex.dyn.attrs = vertex.dyn.attrs,
                graph.attrs = graph.attrs, meta=meta), class="diffnet")
  return(graph)
}

# data("medInnovationsDiffNet")
#
# dn1 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 1]
# dn2 <- medInnovationsDiffNet[medInnovationsDiffNet[["city"]] == 2]
#
# x <- c(dn1, dn2)
