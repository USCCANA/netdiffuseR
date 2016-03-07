#' Creates a \code{diffnet} class object
#'
#' \code{diffnet} objects contain difussion of innovation networks. With adjacency
#' matrices and time of adoption (toa) vector as its main components, most of the
#' package's functions have methods for this class of objects.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param gmode Character scalar. Passed to \code{\link[sna:gplot]{gplot}.}
#' @param toa Numeric vector of size \eqn{n}. Times of adoption.
#' @param t0 Integer scalar. Passed to \code{\link{toa_mat}}.
#' @param t1 Integer scalar. Passed to \code{\link{toa_mat}}.
#' @param undirected Logical scalar.
#' @param self Logical scalar.
#' @param multiple Logical scalar.
#' @param ... In the case of \code{plot}, further arguments passed to \code{gplot}, otherwise
#' is ignored.
#' @param x A \code{diffnet} object.
#' @param object A \code{diffnet} object.
#' @param y Ignored.
#' @param t Integer scalar indicating the time slice to plot.
#' @param displaylabels Logical scalar. When TRUE, \code{plot} shows vertex labels.
#' @param vertex.col Character scalar/vector. Color of the vertices.
#' @param vertex.cex Numeric scalar/vector. Size of the vertices.
#' @param edge.col Character scalar/vector. Color of the edges.
#' @param mode Character scalar. Name of the layout algorithm to implement (see details).
#' @param layout.par Layout parameters (see details).
#' @param main Character. A title template to be passed to sprintf.
#' @param i Indices specifying elements to replace. See \code{\link[base:Extract]{Extract}}.
#' @param value In the case of \code{diffnet.toa}, replacement, otherwise see below.
#' @param vertex.dyn.attrs List of length \eqn{T}. Contains either matrices or data
#' frames with vertex attributes.
#' @param vertex.static.attrs Either a data frame or a matrix with \eqn{n} rows.
#' @param graph.attrs Numeric matrix with a single row.
#' @param slices Either an integer or character vector. While integer vectors are used as
#' indexes, character vectors are used jointly with the time period labels.
#' @param element Character vector/scalar. Indicates what to retrieve/alter.
#' @param attr.class Character vector/scalar. Indicates the class of the attribute, either dynamic (\code{"dyn"}),
#' or static (\code{"static"}).
#' @param as.df Logical scalar. When TRUE returns a data.frame.
#' @param no.print Logical scalar. When TRUE suppress screen messages.
#' @param skip.moran Logical scalar. When TRUE Moran's I is not reported (see details).
#' @param valued Logical scalar. When FALSE non-zero values in the adjmat are set to one.
#' @export
#' @seealso Default options are listed at \code{\link{netdiffuseR-options}}
#' @details Plotting is done via the function \code{\link[sna:gplot]{gplot}},
#' and its layout via \code{\link[sna:gplot.layout]{gplot.layout}}, both from
#' the (\pkg{sna}) package.
#'
#' \code{vertex.cex} can either be a numeric scalar, a numeric vector or a character
#' scalar taking any of the following values \code{"degree"}, \code{"indegree"}, or
#' \code{"outdegree"}. The later will be passed to \code{\link{dgr}} to calculate
#' degree of the selected slice and will be normalized as
#'
#' \deqn{vertex.cex = d/[max(d) - min(d)]\times 2 + .5}{vertex.cex = d/[max(d) - min(d)]* 2 + .5}
#'
#' where \code{d=sqrt(dgr(graph))}.
#'
#' In the case of the \code{summary} method, Moran's I is calculated over the
#' cumulative adoption matrix using as weighting matrix the inverse of the geodesic
#' distance matrix. All this via \code{\link{moran}}. For each time period \code{t},
#' this is calculated as:
#'
#' \preformatted{
#'  m = moran(C[,t], G^(-1))
#' }
#'
#' Where \code{C[,t]} is the t-th column of the cumulative adoption matrix,
#' \code{G^(-1)} is the element-wise inverse of the geodesic matrix at time \code{t},
#' and \code{moran} is \pkg{netdiffuseR}'s moran's I routine. When \code{skip.moran=TRUE}
#' Moran's I is not reported. This can be useful when the graph is particuarly
#' large (tens of thousands of vertices) as when doing so geodesic distances are
#' not calculated, which avoids allocating a square matrix of size \eqn{n} on
#' the memory. As a difference from the adjacency matrices, the matrix with the
#' geodesic distances can't be stored as a sparse matrix (saving space).
#'
#' @section Auxiliary functions:
#'
#' \code{diffnet.attrs} Allows retriving network attributes. In particular, by default
#' returns a list of length \eqn{T} with data frames with the following columns:
#'
#' \enumerate{
#'  \item \code{per} Indicating the time period to which the observation corresponds.
#'  \item \code{toa} Indicating the time of adoption of the vertex.
#'  \item Further columns depending on the vertex and graph attributes.
#' }
#'
#' Each vertex static attributes' are repeated \eqn{T} times in total so that these
#' can be binded (\code{rbind}) to dynamic attributes.
#'
#' When \code{as.df=TRUE}, this convenience function is useful as it can be used
#' to create event history (panel data) datasets used for model fitting.
#'
#' Conversely, the replacement method allows including new vertex or graph
#' attributes either dynamic or static (see examples below).
#'
#' \code{diffnet.toa(graph)} works as an alias of \code{graph$toa}.
#' The replacement method, \code{diffnet.toa<-} used as \code{diffnet.toa(graph)<-...},
#' is the right way of modifying times of adoption as when doing so it
#'  performs several checks on the time ranges, and
#' recalculates adoption and cumulative adoption matrices using \code{toa_mat}.
#'
#'
#' @family diffnet methods
#' @aliases diffnet diffnet-class
#' @examples
#'
#' # Creating a random graph
#' set.seed(123)
#' graph <- rgraph_ba(t=9)
#' graph <- lapply(1:5, function(x) graph)
#'
#' # Pretty TOA
#' names(graph) <- 2001L:2005L
#' toa <- sample(c(2001L:2005L,NA), 10, TRUE)
#'
#' # Creating diffnet object
#' diffnet <- as_diffnet(graph, toa)
#' diffnet
#' summary(diffnet)
#'
#' # Plotting slice 4
#' plot(diffnet, t=4)
#'
#' # ATTRIBUTES ----------------------------------------------------------------
#'
#' # Retrieving attributes
#' diffnet.attrs(diffnet, "vertex", "static")
#'
#' # Now as a data.frame (only static)
#' diffnet.attrs(diffnet, "vertex", "static", as.df = TRUE)
#'
#' # Now as a data.frame (all of them)
#' diffnet.attrs(diffnet, as.df = TRUE)
#'
#' @return
#' A list of class \code{diffnet} with the following elements:
#' \item{graph}{A list of length \eqn{T}. Containing sparse square matrices of size \eqn{n}
#' and class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.}
#' \item{toa}{An integer vector of size \eqn{T} with times of adoption.}
#' \item{adopt, cumadopt}{Numeric matrices of size \eqn{n\times T}{n*T} as those returned
#' by \code{\link{toa_mat}}.}
#' \item{vertex.static.attrs}{If not NULL, a data frame with \eqn{n} rows with vertex static
#' attributes.}
#' \item{vertex.dyn.attrs}{A list of length \eqn{T} with data frames containing vertex attributes
#' throught time (dynamic).}
#' \item{graph.attrs}{If not NULL, a numeric matrix with 1 row containing graph attributes.}
#' \item{meta}{A list of length 9 with the following elements:
#' \itemize{
#'  \item \code{type}: Character scalar equal to \code{"dynamic"}.
#'  \item \code{class}: Character scalar equal to \code{"list"}.
#'  \item \code{ids}: Character vector of size \eqn{n} with vertices' labels.
#'  \item \code{pers}: Integer vector of size \eqn{T}.
#'  \item \code{nper}: Integer scalar equal to \eqn{T}.
#'  \item \code{n}: Integer scalar equal to \eqn{n}.
#'  \item \code{self}: Logical scalar.
#'  \item \code{undirected}: Logical scalar.
#'  \item \code{multiple}: Logical scalar.
#' }
#' }
#' @author Vega Yon
as_diffnet <- function(graph, toa, t0=min(toa, na.rm = TRUE), t1=max(toa, na.rm = TRUE),
                       vertex.dyn.attrs = NULL, vertex.static.attrs= NULL,
                       graph.attrs = NULL,
                       undirected=getOption("diffnet.undirected"),
                       self=getOption("diffnet.self"), multiple=getOption("diffnet.multiple")) {

  # Step 0.0: Check if its diffnet!
  if (inherits(graph, "diffnet")) {
    message("Nothing to do, the graph is already of class \"diffnet\".")
    return(graph)
  }

  # Step 1.1: Check graph ------------------------------------------------------
  meta <- classify_graph(graph)
  if (meta$type=="static") stop("-graph- should be dynamic.")

  # Step 1.2: Checking that lengths fit
  if (length(toa)!=meta$n) stop("-graph- and -toa- have different lengths (",
                                meta$n, " and ", length(toa), " respectively). ",
                                "-toa- should be of length n (number of vertices).")

  # Vertex attrs
  if (length(vertex.dyn.attrs)) {

    attlen <- lapply(vertex.dyn.attrs, nrow)
    if (any(attlen != meta$n)) stop("-graph- and -vertex.dyn.attrs- have different lengths (",
                                    meta$n, " and ", paste(attlen, collapse=", "), "respectively). ",
                                    "-vertex.dyn.attrs- should have n rows.")

    # Coercing into data.frames
    cnames <- colnames(vertex.dyn.attrs[[1]])
    if (!length(cnames))
      cnames <- sprintf("v.dyn.att%03d", 1:ncol(vertex.dyn.attrs[[1]]))

    # Checking if it is a data.frame or not
    if (!inherits(vertex.dyn.attrs[[1]], "data.frame"))
      warning("-vertex.dyn.attrs- will be coerced to a data.frame.")

    vertex.dyn.attrs <- lapply(vertex.dyn.attrs, function(x) {
      if (!inherits(x, "data.frame")) x <- as.data.frame(x)
      dimnames(x) <- list(meta$ids, cnames)
      x
    })
  }

  # Checking static attributes
  if (length(vertex.static.attrs)) {

    # Must be either a data.frame or a matrix
    if (!any(c("data.frame", "matrix") %in% class(vertex.static.attrs)))
      stop("-vertex.static.attrs- must be either a data.frame or a matrix.")

    attlen <- nrow(vertex.static.attrs)
    if (attlen != meta$n) stop("-graph- and -vertex.static.attrs- have different lengths (",
                               meta$n, " and ", attlen, "respectively). ",
                               "-vertex.static.attrs- should have n rows.")

    # Coercing into matrix
    cnames <- colnames(vertex.static.attrs)
    if (!length(cnames))
      cnames <- sprintf("v.static.att%03d", 1:ncol(vertex.static.attrs))

    # Checking if it is a data.frame or not
    if (!inherits(vertex.static.attrs, "data.frame")) {
      warning("-vertex.static.attrs- will be coerced to a data.frame.")
      vertex.static.attrs <- as.data.frame(vertex.static.attrs)
    }

    dimnames(vertex.static.attrs) <- list(meta$ids, cnames)
  }

  # Step 2.1: Checking class of TOA and coercing if necesary
  if (!inherits(toa, "integer")) {
    warning("Coercing -toa- into integer.")
    toa <- as.integer(toa)
  }

  # Step 2.2: Checking names of toa
  if (!length(names(toa)))
    names(toa) <- meta$ids

  # Step 3.1: Creating Time of adoption matrix -----------------------------------
  mat <- toa_mat(toa, labels = meta$ids, t0=t0, t1=t1)

  # Step 3.2: Verifying dimensions and fixing meta$pers
  tdiff <- meta$nper - ncol(mat[[1]])
  if (tdiff < 0)
    stop("Range of -toa- is bigger than the number of slices in -graph- (",
         ncol(mat[[1]]), " and ", length(graph) ," respectively). ",
         "There must be at least as many slices as range of toa.")
  else if (tdiff > 0)
    stop("Range of -toa- is smaller than the number of slices in -graph- (",
         ncol(mat[[1]]), " and ", length(graph) ," respectively). ",
         "Please provide lower and upper boundaries for the values in -toa- ",
         "using -t0- and -t- (see ?toa_mat).")

  meta$pers <- as.integer(colnames(mat$adopt))

  # Step 4.1: Change the class (or set the names) of the graph
  if (meta$class=="array") {
    graph <- lapply(1:meta$nper, function(x) {
      x <- methods::as(graph[,,x], "dgCMatrix")
      dimnames(x) <- with(meta, list(ids, ids))
      x
    })
    names(graph) <- meta$pers
  } else { # Setting names (if not before)
    if (!length(names(graph)))
      names(graph) <- meta$pers
    for(i in 1:meta$nper)
      dimnames(graph[[i]]) <- with(meta, list(ids, ids))
  }

  # Step 5: Compleating attributes and building the object and returning
  meta$self       <- self
  meta$undirected <- undirected
  meta$multiple   <- multiple

  # Checking attributes (last step)
  if (!length(vertex.dyn.attrs)) {
    vertex.dyn.attrs <-
      lapply(meta$pers, function(x) {
        as.data.frame(
          matrix(ncol=0, nrow=meta$n, dimnames = list(meta$ids, NULL))
          )
      })

    # Labeling
    names(vertex.dyn.attrs) <- meta$pers
  }

  if (!length(vertex.static.attrs)) {
    vertex.static.attrs <- as.data.frame(
      matrix(ncol=0, nrow=meta$n, dimnames = list(meta$ids, NULL))
    )
  }

  return(structure(list(
    graph = graph,
    toa   = toa,
    adopt = mat$adopt,
    cumadopt = mat$cumadopt,
    # Attributes
    vertex.static.attrs = vertex.static.attrs,
    vertex.dyn.attrs    = vertex.dyn.attrs,
    graph.attrs         = graph.attrs,
    meta = meta
  ), class="diffnet"))
}

#' @export
#' @rdname as_diffnet
diffnet.attrs <- function(graph, element=c("vertex","graph"), attr.class=c("dyn","static"),as.df=FALSE) {
  nper <- graph$meta$nper
  pers <- graph$meta$pers
  n    <- graph$meta$n

  # Only for diffnet objects
  if (!inherits(graph, "diffnet")) stopifnot_graph(graph)

  # Checking elements
  if (any(!(element %in% c("vertex", "graph"))))
    stop("-element- should only have 'vertex', and/or 'graph'.")

  # Checking classes
  if (any(!(attr.class %in% c("dyn", "static"))))
    stop("-attr.class- should only have 'dyn', and/or 'static'.")

  # Expanding graph static attr
  g.static <- NULL
  v.dyn    <- NULL
  v.static <- NULL
  if ("graph" %in% element) g.static <- graph$graph.attrs
  if ("vertex" %in% element) {
    if ("dyn"    %in% attr.class) v.dyn    <- graph$vertex.dyn.attrs
    if ("static" %in% attr.class) v.static <- graph$vertex.static.attrs
  }

  # Parsing attributes
  if (!length(g.static)) g.static <- as.data.frame(matrix(ncol=0, nrow=n))
  if (!length(v.static)) v.static <- as.data.frame(matrix(ncol=0, nrow=n))
  if (!length(v.dyn[[1]])) v.dyn <- lapply(1:nper, function(y) as.data.frame(matrix(ncol=0, nrow=n)))

  attrs <- cbind(toa=graph$toa, v.static, g.static)
  out <- lapply(1:nper, function(y) {
    cbind(per=rep(pers[y], n),attrs, v.dyn[[y]])
  })

  if (as.df) {
    out <- do.call(rbind, out)
    return(data.frame(out, id=rep(graph$meta$ids, nper), row.names = NULL))
  }

  names(out) <- graph$meta$pers
  out
}

#' @rdname as_diffnet
#' @export
`diffnet.attrs<-` <- function(graph, element="vertex", attr.class="static", value) {

  .Deprecated("[[<-.diffnet")

  # Checking class
  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")

  # Checking what to add
  if (any(!(element %in% c("vertex", "graph"))))
    stop("-element- should be either 'vertex' or 'graph'.")

  if (any(!(attr.class %in% c("static", "dyn"))))
    stop("-attr.class- should be either 'dyn' or 'static'.")

  if (("vertex" == element) && ("static" == attr.class)) {
    # Checking object class
    if (!(class(value) %in%  c("data.frame","matrix")))
      stop("-value- should be either a matrix or a data.frame.")

    # Checking dimensions
    attlen <- nrow(value)
    if (attlen != graph$meta$n) stop("-graph- and -value- have different lengths (",
                                 graph$meta$n, " and ", attlen, " respectively). ",
                                 "-value- should have n rows.")

    # Checking dimnames and coercing into matrix before including them
    # into the data
    cnames <- colnames(value)

    if (length(graph$vertex.static.attrs)) k <- ncol(graph$vertex.static.attrs)
    else k <- 0

    if (!length(cnames))
      cnames <- sprintf("v.static.att%03d", 1:ncol(value) + k)

    # Checking if it is a data.frame or not
    if (!inherits(value, "data.frame")) {
      warning("-value- will be coerced to a data.frame.")
      value <- as.data.frame(value)
    }

    dimnames(value) <- list(graph$meta$ids, cnames)

    # Adding the values
    if (length(graph$vertex.static.attrs)) graph$vertex.static.attrs <- cbind(graph$vertex.static.attrs, value)
    else graph$vertex.static.attrs <- value

  } else if (("vertex" == element) && ("dyn" == attr.class)) {

    # Act depending on the class of object
    gdim <- unlist(graph$meta[c("n", "nper")])
    if (inherits(value, "matrix") | inherits(value, "data.frame")) {

      # Checking dimensions
      test <- which(dim(value) != gdim)
      if (length(test))
        stop("Incorrect dimensions. The -value- must be a data.frame/matrix of size ",
             gdim[1], "x", gdim[2])

      # Coercing into a list
      value <- lapply(seq_len(gdim[2]), function(x) value[,x, drop=FALSE])

    } else if (inherits(value, "list")) {
      # Checking all elements are data.frames/matrices
      test <- which(sapply(value, function(x)
        (!inherits(x, "matrix") & !inherits(x, "data.frame"))))

      if (length(test))
        stop("Some elements of -value- have incorrect class:\n\t",
             paste(test, collapse = ", "), ".")

      # Checking if all elements have the right dimension
      test <- which(sapply(value, function(x) nrow(x) != gdim[1]))
      if (length(test))
        stop("Some of the elements of -value- have incorrect number of rows:\n\t",
             paste(test, collapse=", "), ".")

    } else {
      stop("-value- should be either a matrix/data.frame of size n*T, or a list ",
           "of size T with vectors of length n.")
    }

    # # Checking the length of the attributes
    # if (length(value) != graph$meta$nper)
    #   stop("The length -value-, ",length(value),
    #        ", must coincide with the number of periods, ",graph$meta$nper,".")
    #
    # attlen <- lapply(value, nrow)
    # if (any(attlen != graph$meta$n)) stop("-graph- and -value- have different lengths (",
    #                                       graph$meta$n, " and ", paste(attlen, collapse=", "), "respectively). ",
    #                                 "-value- should have n rows.")

    # Coercing into matrices
    cnames <- colnames(value[[1]])

    if (length(graph$vertex.dyn.attrs[[1]])) k <- ncol(graph$vertex.dyn.attrs[[1]])
    else k <- 0

    if (!length(cnames))
      cnames <- sprintf("v.static.att%03d", k + 1)

    # Checking if it is a data.frame or not
    if (!inherits(value[[1]], "data.frame"))
      warning("-value- will be coerced to a data.frame.")

    value <- lapply(value, function(y) {
      if (!inherits(y, "data.frame")) y <- as.data.frame(y)
      dimnames(y) <- list(graph$meta$ids, cnames)
      y
    })

    # Adding the values
    if (k) {
      for (i in 1:graph$meta$nper)
        graph$vertex.dyn.attrs[[i]] <- cbind(graph$vertex.dyn.attrs[[i]], value[[i]])
    } else {
      for (i in 1:graph$meta$nper)
        graph$vertex.dyn.attrs[[i]] <- value[[i]]
    }
  }

  graph

}

#' @rdname as_diffnet
#' @export
diffnet.toa <- function(graph) {
  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")
  graph$toa
}

#' @rdname as_diffnet
#' @export
`diffnet.toa<-` <- function(graph, i, value) {
  if (!inherits(graph, "diffnet")) stop("-graph- must be a 'diffnet' object")
  if (missing(i)) i <- 1:graph$meta$n

  # Checking values of the data: normalizing
  test <- !(value %in% c(graph$meta$pers, NA))
  if (any(test)) stop("Some elements of -value- (",
                      paste0(head(value[test], 20), collapse=", "),
                      ifelse(length(value[test]) > 20,", ...", "")
                      ,") are not within the range of the original graph.")

  # Changing the value of toa
  graph$toa[i] <- value

  # Recalculating adopt and cumadopt
  mat <- toa_mat(graph$toa, t0=graph$meta$pers[1], t1=graph$meta$pers[graph$meta$nper])

  # checking stack
  nper <- ncol(mat[[1]])
  graph$adopt    <- mat$adopt
  graph$cumadopt <- mat$cumadopt

  graph

}

