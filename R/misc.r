#' Recodes an edgelist such that ids go from 1 to n
#' @param data Edgelist as either a matrix or dataframe with ego and alter
#' @param ... Further arguments for the method (ignored)
#' @return A recoded edgelist as a two-column matrix/data.frame depending
#' on the class of \code{data}. The output includes an attribute called "recode"
#' which contains a two column data.frame providing a mapping between the
#' previous code and the new code (see the examples)
#' @export
#' @details Required for using most of the package's functions, as ids are used
#' as a reference for accessing elements in adjacency matrices.
#' @seealso \code{\link{edgelist_to_adjmat}}
#' @examples
#' # Simple example ------------------------------------------------------------
#' edgelist <- cbind(c(1,1,3,6),c(4,3,200,1))
#' edgelist
#' recoded_edgelist <- recode(edgelist)
#' recoded_edgelist
#'
#' # Retrieving the "recode" attribute
#' attr(recoded_edgelist, "recode")
#' @keywords misc
#' @author George G. Vega Yon
recode <- function(data, ...) UseMethod("recode")

#' @rdname recode
#' @export
recode.data.frame <- function(data, ...) {
  dn <- dimnames(data)
  data <- recode.matrix(as.matrix(data), ...)
  output <- as.data.frame(data)
  dimnames(output) <- dn
  attr(output, "recode") <- attr(data, "recode")
  output
}

#' @rdname recode
#' @export
recode.matrix <- function(data, ...) {

  # Checking the size of the matrix
  dn <- dimnames(data)
  data <- as.factor(as.character(as.vector(data)))
  n <- length(data)
  output <- cbind(data[1:(n/2)], data[(n/2+1):n])
  data <- unique(data)

  # Previous order w/ codes
  rc <- data.frame(
    code=as.integer(data), label=as.character(data),
    stringsAsFactors = FALSE)

  rc <- rc[order(rc[,1]),]

  # Removing NA
  rc <- rc[!is.na(rc$code),]

  dimnames(output) <- dn

  attr(output, "recode") <- rc
  output
}

#' Pretty numbers within a range.
#'
#' A wrapper for \code{\link[base:pretty]{pretty}}.
#'
#' @param x Numeric vector passed to \code{\link[base:pretty]{pretty}}.
#' @param min.n Integer scalar passed to \code{\link[base:pretty]{pretty}}.
#' @param xrange Numeric vector of length 2. Indicates the range in which the
#'  output vector should lie on.
#' @param ... Further arguments passed to the method.
#'
#' The only difference with \code{pretty} is that this function subsets the
#' resulting vector as
#'
#' \code{tick[(tick >= xrange[1]) & (tick <= xrange[2])]}
#'
#' @examples
#' # Simple example ------------------------------------------------------------
#' set.seed(3331)
#' x <- runif(10)
#' pretty(x)
#' pretty_within(x)
#' range(x)
#' @return A vector sequence of `n + 1` round values in the specified range.
#' @export
#' @keywords misc
pretty_within <- function(x, min.n=5, xrange=range(x, na.rm = TRUE), ...) {
    tick <- pretty(x, min.n = min.n, ...)
    tick[(tick >= xrange[1]) & (tick <= xrange[2])]
}

#' Draw a color key in the current device
#'
#' @param x A numeric vector with the data (it is used to extract the range).
#' @param main Character scalar. Title of the key.
#' @param tick.marks A numeric vector indicating the levels to be included in the axis.
#' @param key.pos A numeric vector of length 4 with relative coordinates of the
#'  key (as \% of the plotting area, see  \code{\link[graphics:par]{par("usr")}})
#' @param pos Integer scalar. Position of the axis as in \code{\link[graphics:text]{text}}.
#' @param nlevels Integer scalar. Number of levels (colors) to include in the color key.
#' @param color.palette Color palette of \code{length(nlevels)}.
#' @param tick.width Numeric vector of length 2 indicating the length of the inner
#'  and outer tick marks as percentage of the axis.
#' @param labels Character vector. When provided, specifies using different
#' labels for the tick marks than those provided by \code{tick.marjks}.
#' @param add.box Logical scalar. When \code{TRUE} adds a box around the key.
#' @param na.col Character scalar. If specified, adds an aditional box indicating the NA color.
#' @param na.height Numeric scalar. Relative height of the NA box. Only use if
#' \code{na.col} is not \code{NULL}.
#' @param na.lab Character scalar. Label of the \code{NA} block. Only use if
#' \code{na.col} is not \code{NULL}.
#' @param ... Further arguments to be passed to \code{\link[graphics:rect]{rect}}
#' @export
#' @return Invisible \code{NULL}.
#' @examples
#' set.seed(166)
#' x <- rnorm(100)
#' col <- colorRamp(c("lightblue", "yellow", "red"))((x - min(x))/(max(x) - min(x)))
#' col <- rgb(col, maxColorValue = 255)
#' plot(x, col=col, pch=19)
#' drawColorKey(x, nlevels = 100, border="transparent",
#'  main="Key\nLike A\nBoss")
#' @family visualizations
#' @author George G. Vega Yon
#' @keywords misc
drawColorKey <- function(
  x,
  tick.marks    = pretty_within(x),
  labels        = tick.marks,
  main          = NULL,
  key.pos       = c(.925,.975,.05,.95),
  pos           = 2,
  nlevels       = length(tick.marks),
  color.palette = viridisLite::viridis(nlevels),
  tick.width    = c(.01,.0075),
  add.box       = TRUE,
  na.col        = NULL,
  na.height     = .1,
  na.lab        = "n/a",
  ...) {

  # Checking the pos argument
  test <- which((key.pos > 1) | (key.pos < 0))
  if (length(test))
    stop("Invalid -key.pos-. All values should be in the interval [0,1].")

  # Getting the coordinates
  coords <- par("usr")
  ranges <- list(x = coords[2] - coords[1], y = coords[4] - coords[3])
  coords <- coords + with(ranges, c(
    x*key.pos[1], -x*(1-key.pos[2]),
    y*key.pos[3], -y*(1-key.pos[4])
  ))

  # Giving an space for the NA, putting the starting point
  # a 10% higher.
  if (length(na.col)) {
    na.coords <- coords[3]
    coords[3] <- coords[3] + (coords[4] - coords[3])*na.height
  } else {
    na.coords <- 0
  }

  # Adjusting for text
  if (length(main)) {
    nlines <- length(strsplit(main, "\n")[[1]])
    coords[4] <- coords[4] - par("cxy")[2]*nlines
    # nchars <- min(max(nchar(tick.marks)),5)
    # coords[1] <- coords[1] + par("cxy")
  }

  s <- seq(coords[3], coords[4], length.out = nlevels + 1)
  rcoords <- data.frame(
    x0 = coords[1],
    x1 = coords[2],
    y0 = s[-length(s)],
    y1 = s[-1]
  )

  # Drawing rectangles
  with(rcoords, graphics::rect(x0, y0, x1, y1, col=color.palette, ...))

  # Drawing tickmarks
  tw <- tick.width*ranges$x

  # Adjusting the scale
  xran <- range(x, na.rm = TRUE)
  atick.marks <- (tick.marks - xran[1])/(xran[2] - xran[1])*(coords[4] - coords[3]) + coords[3]

  i <- switch (pos,NA,1,NA,2)
  sgn <- ifelse(pos==2, 1, -1)

  # Drawing the box
  if (!add.box)
    graphics::segments(
      x0 = coords[i],
      y0 = ifelse(length(na.col), na.coords, coords[3]),
      x1 = coords[i], y1 = coords[4])

  graphics::segments(coords[i] - tw[1]*sgn, atick.marks, coords[i] + tw[2]*sgn, atick.marks)
  graphics::text(coords[i] - sgn*tw[1], atick.marks, labels = labels, pos=pos)

  # Adding box
  if (add.box)
    graphics::rect(
      xleft = coords[1],
      ybottom = ifelse(length(na.col), na.coords, coords[3]),
      xright = coords[2], ytop = coords[4])

  # Adding the NA if any
  if (length(na.col)) {

    # Drawing the rectangle
    graphics::rect(
      xleft   = coords[1],
      ybottom = na.coords,
      xright  = coords[2],
      ytop    = coords[3],
      col     = na.col,
      border  = "transparent"
    )

    # Adding text
    # graphics::segments(
    #   x0 = coords[i] - tw[1]*sgn,
    #   y0 = (coords[3] + na.coords)/2,
    #   x1 = coords[i] + tw[2]*sgn,
    #   y1 = (coords[3] + na.coords)/2
    # )

    # graphics::text(coords[i] - sgn*tw[1], atick.marks, labels = labels, pos=pos)
    graphics::text(
      x = coords[i] - sgn*tw[1],
      y = (coords[3] + na.coords)/2,
      labels = na.lab,
      pos = pos
    )
  }

  if (length(main))
    graphics::text(
      x=(coords[2] + coords[1])/2,
      y=coords[4], pos=3, labels=main
      )

  invisible(NULL)
}


#' Rescale vertex size to be used in \code{\link[igraph:plot.igraph]{plot.igraph}}.
#'
#' This function rescales a vertex size before passing it to
#' \code{\link[igraph:plot.igraph]{plot.igraph}} so that the resulting vertices
#' have the desired size relative to the x-axis.
#' @param vertex.size Numeric vector of unscaled vertices' sizes. This is unit-free.
#' @param par.usr Integer vector of length 4 with the coordinates of plotting region.
#'  by default uses \code{par("usr")}.
#' @param minmax.relative.size A numeric vector of length 2. Represents the
#'  desired min and max vertex sizes relative to the x-axis in terms of percentage
#'  (see details).
#' @param adjust Numeric scalar. Adjustment made to the resulting adjusted size
#'  (see details).
#'
#' @details
#' \code{minmax.relative.size} limits the minimum and maximum size that a vertex
#' can take in the plot relative to the x-axis scale. The values for the x-axis
#' scale are by default retrieved by accessing to \code{par("usr")}. By default
#' the vertex are rescaled to be at least 1\% of the size of the plotting region
#' and no more than 5\% of the plotting region, \code{minmax.relative.size=c(.01, .05)}.
#'
#' The default value for \code{adjust} is taken from \code{\link[igraph:igraph]{igraph}}
#' version 1.0.1. In particular, the function \code{igraph:::.igraph.shape.circle.plot},
#' in which before passing the \code{vertex.size} to the function
#' \code{\link[graphics:symbols]{symbols}}, the vertex size is reduced by 200.
#'
#' The rescaling is as follows:
#' \deqn{%
#'  v' = \frac{v - \underbar v}{\bar v - \underbar v}\times (\bar s - \underbar s) + \underbar s
#' }{%
#'  v' = (v - v_min)/(v_max - v_min) * (s_max - s_min) + s_min
#' }
#'
#' Where \eqn{v} is the vertex size, \eqn{\bar v}{v_max} and \eqn{\underbar v}{v_min} are
#' the max and min values of \eqn{v} respectively, and \eqn{\bar s}{s_max} and
#' \eqn{\underbar s}{s_min} are the max and min size that vertices take in terms
#' of \code{minmax.relative.size} and \code{par.usr}. The adjusted value \eqn{v'}
#' is then multiplied by \code{adjust}.
#'
#' \code{igraph_vertex_rescale} and \code{vertex_rescale_igraph} are aliases.
#'
#' @return An integer vector of the same length as \code{vertex.size} with
#' rescaled values.
#' @export
#' @examples
#'
#' library(igraph)
#'
#' # Random graph and coordinates
#' set.seed(2134)
#' g <- barabasi.game(10)
#' coords <- layout_nicely(g)
#'
#' # Random size and figures
#' size <- runif(10)
#' size <- cbind(size, size)
#' shap <- sample(c("circle", "square"),10,TRUE)
#'
#' # Plotting
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow=c(2,2), mai=rep(.5,4))
#' for (i in seq(1, 1000, length.out = 4)) {
#'   # New plot-window
#'   plot.new()
#'   plot.window(xlim=range(coords[,1]*i), ylim=range(coords[,2]*i))
#'
#'   # plotting graph
#'   plot(g, layout=coords*i, add=TRUE, rescale=FALSE,
#'        vertex.shape = shap,
#'        vertex.size  = rescale_vertex_igraph(size) # HERE WE RESCALE!
#'   )
#'
#'   # Adding some axis
#'   axis(1, lwd=0, lwd.ticks = 1)
#'   axis(2, lwd=0, lwd.ticks = 1)
#'   box()
#' }
#'
#' par(oldpar)
#'
#' @family visualizations
#' @author George G. Vega Yon
rescale_vertex_igraph <- function(
  vertex.size,
  par.usr=par("usr"),
  minmax.relative.size= getOption("diffnet.minmax.relative.size", c(0.01, 0.04)),
  adjust=200
) {
  if (!length(vertex.size)) return(
    rescale_vertex_igraph(1, par.usr, minmax.relative.size, adjust))

  # Adjusting x
  xrange <- range(vertex.size)
  xscale <- (par.usr[2] - par.usr[1])*minmax.relative.size
  vertex.size     <- (vertex.size - xrange[1] + 1e-15)/(xrange[2] - xrange[1] + 1e-15)*
    (xscale[2] - xscale[1]) + xscale[1]

  return(vertex.size*adjust/2)
}

#' @rdname rescale_vertex_igraph
#' @export
igraph_vertex_rescale <- rescale_vertex_igraph

#' @rdname rescale_vertex_igraph
#' @export
vertex_rescale_igraph <- rescale_vertex_igraph

# Function to create vertex size accordingly to the
compute_vertex_size <- function(x, vertex.size, slice=1L) {

  # If it is null
  if (!length(vertex.size))
    return(rep(1L, nnodes(x)))

  # If it is of length 1
  if (length(vertex.size) == 1 && is.character(vertex.size)) {

    # Matching degree
    cmodes <- c("indegree", "degree", "outdegree")
    cmode <- cmodes[pmatch(vertex.size, cmodes)]

    if (is.na(cmode))
      stop("Invalid -vertex.size-.\"",vertex.size,"\" is not supported, it should be either \"",
           paste(cmodes, collapse = "\", \""), ".")

    # Repeating the values
    return(dgr(x, cmode = cmode)[,slice])
  }

  # Applyging the function accordignly
  if (is.numeric(vertex.size)) {

    if (length(vertex.size) == 1) return(rep(vertex.size, nnodes(x)))
    else return(vertex.size)

  } else
    stop("Invalid -vertex.size-. It cannot be of class -", class(vertex.size),
         "-. Valid values are either a numeric vector/scalar, or any of ",
         "\"degree\", \"indegree\", or \"outdegree\".")

}

add_dimnames.mat <- function(x) {
  # Getting the parameters
  env <- parent.frame()
  x   <- as.character(match.call()$x)

  # Updating row and colnames if necesary
  if (!length(rownames(env[[x]])))
    rownames(env[[x]]) <- 1L:nrow(env[[x]])
  if (!length(colnames(env[[x]])))
    colnames(env[[x]]) <- 1L:ncol(env[[x]])
  if (length(dim(env[[x]])) == 3 && length(dimnames(env[[x]])) != 3)
    dimnames(env[[x]])[[3]] <- 1L:nslices(env[[x]])

}

add_dimnames.list <- function(x) {

  # Getting the call and the environments
  env <- parent.frame()
  x   <- as.character(match.call()$x)

  # Updating row and colnames if necesary
  for (i in 1:length(env[[x]])) {
    if (!length(rownames(env[[x]][[i]])))
      rownames(env[[x]][[i]]) <- 1L:nrow(env[[x]][[i]])
    if (!length(colnames(env[[x]][[i]])))
      colnames(env[[x]][[i]]) <- 1L:ncol(env[[x]][[i]])
  }

  if (!length(names(env[[x]])))
    names(env[[x]]) <- 1L:nslices(env[[x]])

}




#' Coerce a matrix-like objects to \code{dgCMatrix} (sparse matrix)
#'
#' This helper function allows easy coercion to sparse matrix objects
#' from the \pkg{Matrix} package, \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
#'
#' @param x An object to be coerced into a sparse matrix.
#' @param make.dimnames Logical scalar. When \code{TRUE}, it makes sure that the
#' returned object has dimnames.
#' @param ... Further arguments passed to the method.
#'
#' @details
#' In the case of the \code{igraph} and \code{network} methods, \code{...} is passed to
#' \code{\link[igraph:as_adj]{as_adj}} and \code{\link[network:as.matrix.network]{as.matrix.network}}
#' respectively.
#'
#' @return Either a list with \code{dgCMatrix} objects or a \code{dgCMatrix} object.
#'
#' @export
#'
#' @examples
#'
#' set.seed(1231)
#' x <- rgraph_er(10)
#'
#' # From matrix object
#' as_dgCMatrix(as.matrix(x))
#'
#' # From a network object
#' as_dgCMatrix(network::as.network(as.matrix(x)))
#'
#' # From igraph object
#' as_dgCMatrix(igraph::graph_from_adjacency_matrix(x))
#'
#' # From array
#' myarray <- array(dim=c(10,10,2))
#' myarray[,,1] <- as.matrix(x)
#' myarray[,,2] <- as.matrix(x)
#'
#' myarray
#' as_dgCMatrix(myarray)
#'
#' # From a diffnet object
#' ans <- as_dgCMatrix(medInnovationsDiffNet)
#' str(ans)
#'
#'
as_dgCMatrix <- function(x, make.dimnames = TRUE, ...) {
  UseMethod("as_dgCMatrix")
}

#' @export
#' @rdname as_dgCMatrix
as.dgCMatrix <- as_dgCMatrix

#' @export
#' @rdname as_dgCMatrix
as_spmat <- as_dgCMatrix

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.default <- function(x, make.dimnames = TRUE, ...) {

  ans <- methods::as(x, "dgCMatrix")

  # Updating row and colnames if necesary
  if (make.dimnames)
    add_dimnames.mat(ans)

  ans
}

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.diffnet <- function(x, make.dimnames = TRUE, ...) {
  ans <- x$graph
  ans <- lapply(ans, `dimnames<-`, list(x$meta$ids,x$meta$ids))
  ans
}

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.array <- function(x, make.dimnames = TRUE, ...) {

  # Matrices need a special treatment
  if (inherits(x, "matrix"))
    return(as_dgCMatrix.default(x))

  ans <- apply(x, 3, methods::as, Class="dgCMatrix")

  # Updating row and colnames if necesary
  if (make.dimnames)
    add_dimnames.list(ans)

  ans

}

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.igraph <- function(x, make.dimnames = TRUE, ...) {
  as_dgCMatrix(igraph::as_adj(x, ...), make.dimnames = make.dimnames)
}

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.network <- function(x, make.dimnames = TRUE, ...) {
  as_dgCMatrix(network::as.matrix.network(x), make.dimnames = make.dimnames)
}

#' @export
#' @rdname as_dgCMatrix
as_dgCMatrix.list <- function(x, make.dimnames=TRUE, ...) {

  lapply(x, as_dgCMatrix, make.dimnames = make.dimnames)

}


subtitle <- function(
  x,
  coords = graphics::par("usr")[c(1, 4)] - c(0, graphics::strheight(x)*1.5/2),
  pos    = 4,
  cex    = 1
  ) {

  text(x = coords[1], y = coords[2], labels = x, pos = pos, cex=cex, offset=0)

}
