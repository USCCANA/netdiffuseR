# vertex.attrs <- function(x)
# edge.attrs <- function(x)
#
# `[<-.diffnet` <- function(graph, v, value) {
#   graph$vertex.attrs[v] <- value
#   graph
# }

#' @export
#' @rdname as_diffnet
plot.diffnet <- function(
  x,y=NULL, t=1, displaylabels = FALSE, vertex.col = c("blue", "grey"),
  gmode=ifelse(x$meta$undirected, "graph", "digraph"),
  vertex.cex = "degree", edge.col = "gray", mode = "fruchtermanreingold",
  layout.par = NULL, main = "Diffusion network in time %d", ...) {

  # Checking that the time period is actually within
  if (!(t %in% 1:x$meta$nper))
    stop("-t- must be an integer within 1 and ",x$meta$nper,".")

  # Extracting the graph to be plotted
  graph <- methods::as(x$graph[[t]], "matrix.csc")

  # Setting the colors
  cols <- with(x, ifelse(cumadopt[,t], vertex.col[1], vertex.col[2]))

  # Calcularing layout
  fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")
  coords <- fun(graph, layout.par)

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character"))
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(x$graph[[t]], undirected = x$meta$undirected)
      vertex.cex <- sqrt(vertex.cex)
      r <- range(vertex.cex)

      # If all the vertices have the same degree
      vertex.cex <- (vertex.cex - r[1]+ .1)/(r[2] - r[1] + .1)*2
    } else {
      stop("Invalid -vertex.cex-")
    }

  sna::gplot(graph, displaylabels=displaylabels, vertex.col=cols,
             coord=coords, edge.col=edge.col, label=x$meta$ids,
             main=sprintf(main, x$meta$pers[t]), vertex.cex=vertex.cex, gmode=gmode, ...)

  invisible(coords)
}

#' @export
#' @rdname as_diffnet
print.diffnet <- function(x, ...) {
  with(x, cat(
    "Dynamic network of class -diffnet-",
    paste(" # of nodes        :", meta$n),
    paste(" # of time periods :", meta$nper),
    paste(" Final prevalence  :",
          formatC(sum(cumadopt[,meta$nper])/meta$n, digits = 2, format="f")
          ),
    paste(" Type              :", ifelse(meta$undirected, "undirected", "directed")),
    sep="\n"
    )
  )
  invisible(x)
}

#' @export
#' @rdname as_diffnet
summary.diffnet <- function(object, slices=NULL, no.print=FALSE,
                            skip.moran=FALSE, valued=getOption("diffnet.valued",FALSE), ...) {
  # Subsetting
  if (!length(slices)) slices <- 1:object$meta$nper

  # If no valued
  if (!valued)
    for (i in 1:object$meta$nper)
      object$graph[[i]]@x <- rep(1, length(object$graph[[i]]@x))

  # Checking that the time period is actually within
  test <- !(slices %in% 1:object$meta$nper)
  if (any(test))
    stop("-slices- must be an integer range within 1 and ",object$meta$nper,".")

  slices <- sort(slices)

  # To make notation nicer
  meta <- object$meta

  # Computing density
  d <- unlist(lapply(object$graph[slices], function(x) {
    nelements <- length(x@x)
    x <-nelements/(meta$n * (meta$n-1))
  }))

  # Computing moran's I
  if (!skip.moran) {
    m <- vector("numeric", length(slices))
    for (i in 1:length(slices)) {
      g <- 1/(1e-15 + sna::geodist(as.matrix(object$graph[[slices[i]]]), meta$n,
              count.paths=FALSE)$gdist)
      diag(g) <- 0
      m[i] <- moran(object$cumadopt[,slices[i]], g)
    }
  }
  # Computing adopters, cumadopt and hazard rate
  ad <- colSums(object$adopt[,slices,drop=FALSE])
  ca <- t(cumulative_adopt_count(object$cumadopt))[slices,-3, drop=FALSE]
  hr <- t(hazard_rate(object$cumadopt, no.plot = TRUE))[slices,,drop=FALSE]

  # Left censoring
  lc <- sum(object$toa == meta$pers[1], na.rm = TRUE)
  rc <- sum(is.na(object$toa), na.rm=TRUE)

  out <- data.frame(
    adopt = ad,
    cum_adopt = ca[,1],
    cum_adopt_pcent = ca[,2],
    hazard = hr,
    density=d
  )

  if (!skip.moran) out$moran <- m

  if (no.print) return(out)

  # Function to print data.frames differently
  header <- c(" Period ","Adopters","Cum Adopt.", "Cum Adopt. %",
              "Hazard Rate","Density", if (!skip.moran) "Moran's I" else NULL)
  slen   <- nchar(header)
  hline  <- paste(sapply(sapply(slen, rep.int, x="-"), paste0, collapse=""),
                  collapse=" ")
  rule   <- paste0(rep("-", sum(slen) + length(slen) - 1), collapse="")

  # Quick Formatting function
  qf <- function(x, digits=2) sprintf(paste0("%.",digits,"f"), x)

  cat("Diffusion network summary statistics\n",rule,"\n",sep="")
  cat(header,"\n")
  cat(hline, "\n")
  for (i in 1:nrow(out)) {
    cat(sprintf(
      paste0("%",slen,"s", collapse=" "),
      qf(meta$pers[slices[i]],0), qf(out[i,1],0), qf(out[i,2],0), qf(out[i,3]),
      ifelse(i==1, "-",qf(out[i,4])), qf(out[i,5]), if (!skip.moran) qf(out[i,6]) else ""
    ), "\n")
  }


  # print(out, digits=2)

  cat(
    rule,
    paste(" Left censoring  :", sprintf("%3.2f (%d)", lc/meta$n, lc)),
    paste(" Right centoring :", sprintf("%3.2f (%d)", rc/meta$n, rc)),
    paste(" # of nodes      :", sprintf("%d",meta$n)),
    sep="\n"
  )

  invisible(out)
}

#' Plot the diffusion process
#'
#' Creates a colored network plot showing the structure of the graph through time
#' (one network plot for each time period)  and the set of adopter and non-adopters
#' in the network.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param cumadopt \eqn{n\times T}{n*T} matrix.
#' @param displaylabels Logical scalar. When TRUE vertex labels are displayed (see \code{\link[sna:gplot]{gplot}})
#' @param slices Integer vector. Indicates what slices to plot. By default all are plotted.
#' @param undirected Logical scalar.
#' @param vertex.col A character vector of size 3 with colors names.
#' @param vertex.cex Numeric vector of size \eqn{n}. Size of the vertices.
#' @param label Character vector of size \eqn{n}. If no provided, rownames of
#' the graph are used.
#' @param edge.col Character scalar/vector. Color of the edge.
#' @param mode Character scalar. Name of the layout algorithm to implement (see details).
#' @param layout.par Layout parameters (see details).
#' @param mfrow.par Vector of size 2 with number of rows and columns to be passed to \code{\link{par}.}
#' @param main Character scalar. A title template to be passed to \code{\link{sprintf}.}
#' @param mai Numeric vector of size 4. To be passed to \code{\link{par}.}
#' @param mar Numeric vector of size 4. To be passed to \code{\link{par}.}
#' @param gmode Character scalar. See \code{\link[sna:gplot]{gplot}.}
#' @param ... Further arguments to be passed to \code{\link[sna:gplot]{gplot}.}
#' @param lgd List of arguments to be passed to \code{\link{legend}}.
#'
#' @details Plotting is done via the function \code{\link[sna:gplot]{gplot}},
#' and its layout via \code{\link[sna:gplot.layout]{gplot.layout}}, both from
#' the (\pkg{sna}) package.
#'
#' In order to center the attention on the diffusion process itself, the
#' positions of each vertex are computed only once by aggregating the networks
#' through time, this is, instead of computing the layout for each time \eqn{t},
#' the function creates a new graph accumulating links through time.
#'
#' The \code{mfrow.par} sets how to arrange the plots on the device. If \eqn{T=5}
#' and \code{mfrow.par=c(2,3)}, the first three networks will be in the top
#' of the device and the last two in the bottom.
#'
#' The argument \code{vertex.col} contains the colors of non-adopters, new-adopters,
#' and adopters respectively. The new adopters (default color \code{"red"}) have a different
#' color that the adopters when the graph is at their time of adoption, hence,
#' when the graph been plotted is in \eqn{t=2} and \eqn{toa=2} the vertex will
#' be plotted in red.
#'
#' \code{vertex.cex} can either be a numeric scalar, a numeric vector or a character
#' scalar taking any of the following values \code{"degree"}, \code{"indegree"}, or
#' \code{"outdegree"}. The later will be passed to \code{\link{dgr}} to calculate
#' degree of the cumulated graph and will be normalized as
#'
#' \deqn{vertex.cex = [d - \min(d) + .1]/[\max(d) - \min(d) + .1]\times 2}{vertex.cex = [d - min(d) + .1]/[max(d) - min(d) + .1]* 2}
#'
#' where \eqn{d=\sqrt{dgr(graph)}}{d=sqrt(dgr(graph))}.
#'
#' @examples
#' # Generating a random graph
#' set.seed(1234)
#' n <- 6
#' nper <- 5
#' graph <- rgraph_er(n,nper, p=.3, undirected = FALSE)
#' toa <- sample(2000:(2000+nper-1), n, TRUE)
#' adopt <- toa_mat(toa)
#'
#' plot_diffnet(graph, adopt$cumadopt)
#' @return Calculated coordinates for the grouped graph (invisible).
#' @family visualizations
#' @keywords hplot
#' @export
#' @author Vega Yon
plot_diffnet <- function(
  graph, cumadopt,
  slices=NULL,
  displaylabels=FALSE,
  undirected=TRUE,
  vertex.col=c("grey", "red", "blue"),
  vertex.cex="degree",
  label=rownames(graph[[1]]),
  edge.col="gray",
  mode="fruchtermanreingold", layout.par=NULL,
  mfrow.par=NULL, main="Network in period %d",
  mai=c(0,0,1,0),
  mar=rep(1,4) + 0.1, gmode=ifelse(undirected, "graph", "digraph"),
  lgd = list(x="center", legend=c("Non adopters", "New adopters","Adopters"), pch=21,
             bty="n", cex=1.2, horiz=TRUE), ...
) {
  switch (class(graph),
    array = plot_diffnet.array(
      graph, cumadopt, slices, displaylabels, undirected, vertex.col, vertex.cex, label,
      edge.col, mode, layout.par, mfrow.par, main, mai, mar, gmode, lgd, ...),
    list = plot_diffnet.list(
      graph, cumadopt, slices, displaylabels, undirected, vertex.col, vertex.cex, label,
      edge.col, mode, layout.par, mfrow.par, main, mai, mar, gmode, lgd, ...),
    diffnet = plot_diffnet.list(
      graph$graph, graph$cumadopt, slices, displaylabels, graph$meta$undirected,
      vertex.col, vertex.cex, label=graph$meta$ids,
      edge.col, mode, layout.par, mfrow.par, main, mai, mar, gmode, lgd,...),
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_diffnet
plot_diffnet.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) graph[,,x])
  names(graph) <- dn
  plot_diffnet.list(graph, ...)
}

# @export
# @rdname plot_diffnet
plot_diffnet.list <- function(graph, cumadopt, slices,
                         displaylabels=FALSE,
                         undirected=TRUE,
                         vertex.col=c("grey", "red", "blue"),
                         vertex.cex="degree",
                         label=rownames(graph[[1]]),
                         edge.col="gray",
                         mode="fruchtermanreingold", layout.par=NULL,
                         mfrow.par=NULL, main="Network in period %d",
                         mai=c(0,0,1,0),
                         mar=rep(1,4) + 0.1,
                         gmode=ifelse(undirected, "graph", "digraph"),
                         lgd = list(x="center", legend=c("Non adopters", "New adopters","Adopters"), pch=21,
                                    bty="n", cex=1.2, horiz=TRUE), ...) {

  # Checking slices
  if (!length(slices)) slices <- 1:ncol(cumadopt)
  graph    <- graph[slices]
  cumadopt <- cumadopt[,slices,drop=FALSE]

  t <- length(graph)
  n <- nrow(graph[[1]])
  cumgraph <- Matrix::sparseMatrix(i={}, j={}, dims=c(n, n))
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[[i]]
  }

  # Getting the coords
  fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")

  # In order to be using SNA functions, we need to coerse the graph
  # into an object from SparseM
  coords <- fun(methods::as(cumgraph, "matrix.csc"), layout.par)

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character"))
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(cumgraph, undirected=undirected)
      vertex.cex <- sqrt(vertex.cex)
      r <- range(vertex.cex)
      vertex.cex <- (vertex.cex - r[1]+ .1)/(r[2] - r[1] + .1)*2
    } else {
      stop("Invalid -vertex.cex-")
    }

#   if ( (length(vertex.cex)==1) && (n > 1) )
#     vertex.cex <- rep(vertex.cex,n)

  # Figuring out the dimension
  if (!length(mfrow.par)) {
    if (t<4) mfrow.par <- c(1,t)
    else if (t==4) mfrow.par <- c(2,2)
    else if (t==5) mfrow.par <- c(2,3)
    else if (t==6) mfrow.par <- c(2,3)
    else if (t==7) mfrow.par <- c(2,4)
    else if (t==8) mfrow.par <- c(2,4)
    else if (t==9) mfrow.par <- c(3,4)
    else if (t==10) mfrow.par <- c(3,4)
    else if (t==11) mfrow.par <- c(3,4)
    else if (t==12) mfrow.par <- c(3,4)
    else mfrow.par <- c(ceiling(t/4),4)
  }

  test <- prod(mfrow.par)-t
  if (test) {
    marlayout <- matrix(1:prod(mfrow.par), ncol=mfrow.par[2], byrow = TRUE)
    marlayout[nrow(marlayout),(ncol(marlayout) - test + 1):ncol(marlayout)] <- t + 1
    lgd$horiz <- FALSE
  } else {
    marlayout <- rbind(matrix(1:prod(mfrow.par), ncol=mfrow.par[2], byrow = TRUE),
                       prod(mfrow.par)+1)
  }

  # Plotting
  curseed <- .Random.seed
  oldpar <- par(no.readonly = TRUE)
  # par(mfrow=mfrow.par, mai=mai, mar=mar)
  par(mai=mai, mar=mar)
  layout(marlayout)

  times <- as.integer(names(graph))
  for(i in 1:t)  {
    # Colors, new adopters are painted differently
    cols <- ifelse(!cumadopt[,i], vertex.col[1],
                   ifelse(!cumadopt[,i-1*(i!=1)] | rep(i,n) == 1, vertex.col[2], vertex.col[3]))

    set.seed(curseed)
    # cgraph <- sna::as.sociomatrix.sna(adjmat_to_edgelist(graph[[i]], undirected))
    if (inherits(graph[[i]], "dgCMatrix")) g <- methods::as(graph[[i]], "matrix.csc")
    else g <- graph[[i]]
    sna::gplot(g,
               displaylabels = displaylabels, vertex.col = cols, coord=coords,
               edge.col = edge.col,vertex.cex = vertex.cex, label=label,
               main=sprintf(main, times[i]), gmode=gmode, ...)
  }

  # Legend
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))

  lgd$pt.bg  <- vertex.col

  do.call(legend, lgd)

#   with(lgd,
#     legend(pos, pt.bg = vertex.col,
#            legend = text, pch=pch, bty=bty, cex=cex, horiz = horiz)
#   )

  par(oldpar)
  invisible(coords)

}

#' Threshold levels through time
#'
#' Draws a graph where the coordinates are given by time of adoption, x-axis,
#' and threshold level, y-axis.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param expo \eqn{n\times T}{n * T} matrix. Esposure to the innovation obtained from \code{\link{exposure}}
#' @param toa Integer vector of size \eqn{n}. Times of Adoption
#' @param t0 Integer scalar. Passed to \code{\link{threshold}}.
#' @param undirected Logical scalar.
#' @param no.contemporary Logical scalar. When TRUE, edges for vertices with the same
#' \code{toa} won't be plotted.
#' @param main Character scalar. Title of the plot.
#' @param xlab Character scalar. x-axis label.
#' @param ylab Character scalar. y-axis label.
#' @param vertex.cex Numeric vector of size \eqn{n}. Relative size of the vertices.
#' @param vertex.col Either a vector of size \eqn{n} or a scalar indicating colors of the vertices.
#' @param vertex.label Character vector of size \eqn{n}. Labels of the vertices.
#' @param vertex.lab.pos Integer value to be passed to \code{\link{text}} via \code{pos}.
#' @param edge.width Numeric. Width of the edges.
#' @param edge.col Character. Color of the edges.
#' @param arrow.length Numeric value to be passed to \code{\link{arrows}}.
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn.
#' @param bty See \code{\link{par}}.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @family visualizations
#' @seealso Use \code{\link{threshold}} to retrieve the corresponding threshold
#' obtained returned by \code{\link{exposure}}.
#' @keywords hplot
#' @examples
#'
#' # Generating a random graph
#' set.seed(1234)
#' n <- 6
#' nper <- 5
#' graph <- rgraph_er(n,nper, p=.3, undirected = FALSE)
#' toa <- sample(2000:(2000+nper-1), n, TRUE)
#' adopt <- toa_mat(toa)
#'
#' # Computing exposure
#' expos <- exposure(graph, adopt$cumadopt)
#'
#' plot_threshold(graph, expos, toa)
#'
#' # Calculating degree (for sizing the vertices)
#' plot_threshold(graph, expos, toa, vertex.cex = "indegree")
#'
#' @export
#' @author Vega Yon
plot_threshold <- function(
  graph, expo, toa, t0=min(toa, na.rm = TRUE),
  undirected=getOption("diffnet.undirected"), no.contemporary=TRUE,
  main="Time of Adoption by Network Threshold", xlab="Time", ylab="Threshold",
  vertex.cex="degree", vertex.col="blue", vertex.label=NULL, vertex.lab.pos=3,
  edge.width = 2, edge.col = "gray", arrow.length=.20,
  include.grid = TRUE, bty="n", ...
) {

  if (missing(expo))
    if (!inherits(graph, "diffnet")) {
      stop("-expo- should be provided when -graph- is not of class 'diffnet'")
    } else {
      expo <- exposure(graph)
    }

  switch (class(graph),
    array = plot_threshold.array(
      graph, expo, toa, t0, undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label, vertex.lab.pos, edge.width, edge.col,
      arrow.length, include.grid, bty, ...),
    list = plot_threshold.list(
      graph, expo, toa, t0, undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label, vertex.lab.pos, edge.width, edge.col,
      arrow.length, include.grid, bty, ...),
    diffnet = {
      # If graph is diffnet, then we should do something different (because the
      # first toa may not be the firts one as toa may be stacked to the right.
      # see ?as_diffnet)
      # graph$toa <- graph$toa - min(graph$meta$pers) + 1L

      plot_threshold.list(
      graph$graph, expo,
      graph$toa, t0=graph$meta$pers[1], graph$meta$undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label, vertex.lab.pos, edge.width, edge.col,
      arrow.length, include.grid, bty, ...)
      },
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_threshold
plot_threshold.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) graph[,,x])
  names(graph) <- dn
  plot_threshold.list(graph, ...)
}

# @export
# @rdname plot_threshold
plot_threshold.list <- function(
  graph, expo=NULL, toa, t0=min(toa, na.rm=TRUE),
  undirected=getOption("diffnet.undirected"), no.contemporary=TRUE,
  main="Time of Adoption by Network Threshold", xlab="Time", ylab="Threshold",
  vertex.cex="degree", vertex.col="blue", vertex.label=NULL, vertex.lab.pos=3,
  edge.width = 2, edge.col = "gray", arrow.length=.20,
  include.grid = TRUE, bty="n", ...) {
  # Step 0: Getting basic info
  t <- length(graph)
  n <- nrow(graph[[1]])

  # Step 1: Creating the cumulative graph
  cumgraph <- Matrix::sparseMatrix(i={}, j={}, dims=c(n, n))
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[[i]]
  }

  # Creating the pos vector
  y <- threshold(expo, toa, t0)

  # Jitter to the xaxis and limits
  jit <- jitter(toa, amount = .25)
  xran <- range(toa, na.rm = TRUE)
  xlim <- xran + c(-1,1)
  yran <- c(0,1)
  ylim <- yran + (yran[2] - yran[1])*.1*c(-1,1)

  # Step 2: Checking colors and sizes

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character")) {
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(cumgraph, undirected = undirected)
      vertex.cex <- sqrt(vertex.cex)
      r <- range(vertex.cex)

      # If all the vertices have the same degree
      vertex.cex <- (vertex.cex - r[1]+ .1)/(r[2] - r[1] + .1)/4
    } else {
      stop("Invalid -vertex.cex-")
    }
  } else if (length(vertex.cex)==1) {
    vertex.cex <- rep(vertex.cex, n)
  } else if (inherits(vertex.cex, "character")) stop("Invalid value for -vertex.cex-.")

  # Plotting
  # oldpar <- par(no.readonly = TRUE)
  plot(NULL, xlim=xlim, ylim=ylim, bty=bty, xlab=xlab, ylab=ylab, main=main, ...)

  # Should there be a grid??
  if (include.grid) grid()

  # Now, for y (it should be different)
  xran <- range(xlim, na.rm = TRUE)
  yran <- range(ylim, na.rm = TRUE)
  vertex.cex.y <- vertex.cex *(yran[2]-yran[1])/(xran[2]-xran[1])

  # Drawing arrows, first we calculate the coordinates of the edges, for this we
  # use the function edges_coords. This considers aspect ratio of the plot.
  edges <- netdiffuseR::edges_coords(cumgraph, toa, jit, y, vertex.cex, undirected, no.contemporary)
  edges <- as.data.frame(edges)

  with(edges, arrows(x0, y0, x1, y1, lwd = edge.width, col = edge.col,
                                length=arrow.length))

  # Drawing the vertices and its labels
  symbols(jit, y, circles=vertex.cex, inches=FALSE, bg=vertex.col, add=TRUE)

  # Positioning labels can be harsh, so we try with this algorithm
  if (!length(vertex.label)) vertex.label <- 1:n
  text(x=jit, y=y+vertex.cex.y, labels = vertex.label, pos=vertex.lab.pos)

  # par(oldpar)

  invisible(data.frame(toa=toa,threshold=y, jit=jit))

}

#' Plot distribution of infect/suscep
#'
#' After calculating infectiousness and susceptibility of each individual on the
#' network, it creates an \code{nlevels} by \code{nlevels} matrix indicating the
#' number of individuals that lie within each cell, and draws a heatmap.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param toa Integer vector of size \eqn{T}. Passed to infection/susceptibility.
#' @param t0 Integer scalar. See \code{\link{toa_mat}}.
#' @param normalize Logical scalar.  Passed to infection/susceptibility.
#' @param K Integer scalar.  Passed to infection/susceptibility.
#' @param r Numeric scalar.  Passed to infection/susceptibility.
#' @param expdiscount Logical scalar.  Passed to infection/susceptibility.
#' @param bins Integer scalar. Size of the grid (\eqn{n}).
#' @param nlevels Integer scalar. Number of levels to plot (see \code{\link{filled.contour}}).
#' @param logscale Logical scalar. When TRUE the axis of the plot will be presented in log-scale.
#' @param main Character scalar. Title of the graph.
#' @param xlab Character scalar. Title of the x-axis.
#' @param ylab Character scalar. Title of the y-axis.
#' @param sub Character scalar. Subtitle of the graph.
#' @param color.palette a color palette function to be used to assign colors in the plot (see \code{\link{filled.contour}}).
#' @param include.grid Logical scalar. When TRUE, the grid of the graph is drawn.
#' @param ... Additional parameters to be passed to \code{\link{filled.contour}.}
#' @param exclude.zeros Logical scalar. When TRUE, observations with zero values
#' @param valued Logical scalar. When FALSE non-zero values in the adjmat are set to one.
#' in infect or suscept are excluded from the graph. This is done explicitly when \code{logscale=TRUE}.
#' @details
#'
#' This plotting function was inspired by Aral, S., & Walker, D. (2012).
#'
#' @return A list with three elements:
#' \item{infect}{A numeric vector of size \eqn{n} with infectiousness levels}
#' \item{suscep}{A numeric vector of size \eqn{n} with susceptibility levels}
#' \item{coords}{A list containing the class marks and counts used to draw the
#' plot via \code{\link{filled.contour}} (see \code{\link{grid_distribution}})}
#' \item{complete}{A logical vector with \code{TRUE} when the case was included in
#' the plot. (this is relevant whenever \code{logscale=TRUE})}
#' @family visualizations
#' @seealso Infectiousness and susceptibility are computed via \code{\link{infection}} and
#' \code{\link{susceptibility}}.
#' @keywords hplot
#' @references
#' Aral, S., & Walker, D. (2012). "Identifying Influential and Susceptible Members
#' of Social Networks". Science, 337(6092), 337–341.
#' \url{http://doi.org/10.1126/science.1215842}
#' @export
#' @examples
#' # Generating a random graph
#' set.seed(1234)
#' n <- 100
#' nper <- 20
#' graph <- rgraph_er(n,nper, p=.2, undirected = FALSE)
#' toa <- sample(1:(1+nper-1), n, TRUE)
#'
#' # Visualizing distribution of suscep/infect
#' out <- plot_infectsuscep(graph, toa, K=3, logscale = TRUE)
#' @author Vega Yon
plot_infectsuscep <- function(
  graph, toa, t0=NULL,normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE, bins=50,nlevels=round(bins/2),
  logscale=TRUE, main="Distribution of Infectiousness and\nSusceptibility",
  xlab="Infectiousness of ego", ylab="Susceptibility of ego",
  sub=ifelse(logscale, "(in log-scale)", NA), color.palette=function(n) grey(n:1/n),
  include.grid=TRUE, exclude.zeros=FALSE, valued=getOption("diffnet.valued",FALSE), ...
) {

  # Checking the times argument
  if (missing(toa))
    if (!inherits(graph, "diffnet")) {
      stop("-toa- should be provided when -graph- is not of class 'diffnet'")
    } else {
      toa <- graph$toa
      t0    <- min(graph$meta$pers)
    }

  if (!length(t0)) t0 <- min(toa, na.rm = TRUE)

  switch (class(graph),
    array = plot_infectsuscep.array(
      graph, toa, t0, normalize, K, r, expdiscount, bins, nlevels, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued, ...),
    list = plot_infectsuscep.list(
      graph, toa, t0, normalize, K, r, expdiscount, bins, nlevels, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued,...),
    diffnet = plot_infectsuscep.list(
      graph$graph, graph$toa, t0, normalize, K, r, expdiscount, bins, nlevels, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued,...),
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_infectsuscep
plot_infectsuscep.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) methods::as(graph[,,x], "dgCMatrix"))
  names(graph) <- dn
  plot_infectsuscep.list(graph, ...)
}

# @export
# @rdname plot_infectsuscep
plot_infectsuscep.list <- function(graph, toa, t0, normalize,
                              K, r, expdiscount,
                              bins,nlevels,
                              logscale,
                              main,
                              xlab,
                              ylab,
                              sub,
                              color.palette,
                              include.grid, exclude.zeros, valued,
                              ...) {
  # Computing infect and suscept
  infect <- infection(graph, toa, t0, normalize, K, r, expdiscount, valued)
  suscep <- susceptibility(graph, toa, t0, normalize, K, r, expdiscount, valued)

  # Performing classification (linear)
  if (logscale) {
    infectp<-log(infect)
    suscepp<-log(suscep)

    # Only keeping complete cases
    complete <- is.finite(infectp) & is.finite(suscepp)
    if (any(!complete)) warning("When applying logscale some observations are missing.")
    infectp <- infectp[complete,]
    suscepp <- suscepp[complete,]
  }
  else {
    infectp <- infect
    suscepp <- suscep
    complete <- vector(length=length(infectp))
  }

  if ((!length(infectp) | !length(suscepp)) & logscale)
    stop("Can't apply logscale (undefined values).")

  # If excluding zeros
  include <- rep(TRUE,length(toa))
  if (exclude.zeros) {
    include[!infectp | !suscepp] <- FALSE
  }

  # Computing infect & suscept
  coords <- netdiffuseR::grid_distribution(x=infectp[include], y=suscepp[include], bins)


  # Nice plot
  n <- sum(coords$z)
  with(coords, filled.contour(
    x,y,
    z/n, bty="n", main=main, xlab=xlab, ylab=ylab, sub=sub, color.palette =color.palette,
    xlim=range(pretty(x)), ylim=range(pretty(y)),
    plot.axes={

      # Preparing the tickmarks for the axis
      xticks  <- pretty(x)
      yticks  <- pretty(y)
      if (logscale) {
        xlticks <- exp(xticks)
        ylticks <- exp(yticks)
      } else {
        xlticks <- xticks
        ylticks <- yticks
      }

      # Drawing the axis
      axis(1, xticks, sprintf("%.2f",xlticks))
      axis(2, yticks, sprintf("%.2f",ylticks))

      # Putting the grid
      if (include.grid) grid()
    }, nlevels=nlevels, ...))
  # if (include.grid) grid()

  invisible(list(infect=infect, suscept=suscep, coords=coords,
                 complete=complete))
}

#' Visualize adopters and cumulative adopters
#' @param obj Either a diffnet object or a cumulative a doption matrix.
#' @param freq Logical scalar. When TRUE frequencies are plotted instead of proportions.
#' @param what Character vector of length 2. What to plot.
#' @param add Logical scalar. When TRUE lines and dots are added to the current graph.
#' @param include.legend Logical scalar. When TRUE a legend of the graph is plotted.
#' @param pch Integer vector of length 2. See \code{\link{matplot}}.
#' @param type Character vector of length 2. See \code{\link{matplot}}.
#' @param ylim Numeric vector of length 2. Sets the plotting limit for the y-axis.
#' @param lty Numeric vector of length 2. See \code{\link{matplot}}.
#' @param col Character vector of length 2. See \code{\link{matplot}}.
#' @param bg Character vector of length 2. See \code{\link{matplot}}.
#' @param xlab Character scalar. Name of the x-axis.
#' @param ylab Character scalar. Name of the y-axis.
#' @param main Character scalar. Title of the plot
#' @param ... Further arguments passed to \code{\link{matplot}}.
#' @param include.grid Logical scalar. When TRUE, the grid of the graph is drawn
#' @family visualizations
#' @examples
#' # Generating a random diffnet
#' set.seed(8321)
#' diffnet <- rdiffnet(20, 5, seed.graph="small-world", seed.nodes="central")
#'
#' plot_adopters(diffnet)
#'
#' # Alternatively, we can use a TOA Matrix
#' toa <- sample(c(NA, 2010,2015), 20, TRUE)
#' mat <- toa_mat(toa)
#' plot_adopters(mat$cumadopt)
#' @return List of matrices as described in \code{\link{cumulative_adopt_count}}
#' @export
#' @author Vega Yon
plot_adopters <- function(obj, freq=FALSE, what=c("adopt","cumadopt"),
                          add=FALSE, include.legend=TRUE, include.grid=TRUE,
                          pch=c(21,24), type=c("b", "b"),
                          ylim=if (!freq) c(0,1) else NULL, lty=c(1,1), col=c("black","black"),
                          bg = c("lightblue","gray"),
                          xlab="Time", ylab=ifelse(freq, "Frequency", "Proportion"),
                          main="Adopters and Cumulative Adopters", ...) {

  # Checking what
  if (any(!(what %in% c("adopt", "cumadopt"))))
    stop("Invalid curve to plot. -what- must be in c(\"adopt\",\"cumadopt\").")

  # Computing the TOA mat
  if (inherits(obj, "diffnet")) {
    cumadopt <- cumulative_adopt_count(obj$cumadopt)
    adopt    <- colSums(obj$adopt)
    n        <- obj$meta$n
  }
  else {
    cumadopt <- cumulative_adopt_count(obj)
    adopt    <- cumadopt["num",] - c(0,cumadopt["num",1:(ncol(cumadopt)-1)])
    n        <- nrow(obj)
  }

  out <- cumadopt

  # In the case that the user wants pcent (the default)
  if (!freq) {
    cumadopt <- cumadopt/n
    adopt    <- adopt/n
  }

  # Time names...
  times <- colnames(cumadopt)
  if ((length(ylim) == 1) && is.na(ylim))
    ylim <- NULL

  # Building matrix to plot
  k <- length(what)
  n <- length(times)
  mat <- matrix(ncol=k, nrow=n)
  if ("cumadopt" %in% what) mat[,1] <- cumadopt["num",]
  if ("adopt" %in% what) mat[,k] <- adopt

  # Fixing parameters
  test <- c("cumadopt" %in% what, "adopt" %in% what)
  if (length(type) > k) type <- type[test]
  if (length(lty)  > k) lty  <- lty[test]
  if (length(col)  > k) col  <- col[test]
  if (length(bg)   > k) bg   <- bg[test]
  if (length(pch)  > k) pch  <- pch[test]

  out <- matplot(times, y=mat, ylim=ylim, add=add, type=type,
          lty=lty, col=col, xlab=xlab, ylab=ylab, main=main, pch=pch,
          bg=bg,...)

  # If not been added
  if (!add) {
    if (include.legend)
      legend("topleft", bty="n",
             legend = c("Cumulative adopters", "Adopters")[test], fill = bg)

    if (include.grid)
      grid()
  }

  invisible(out)
}

# x <- cumulative_adopt_count(diffnet)
# z <- x["num",] - c(0,x["num",1:(ncol(x)-1)])
# cumsum(z)
# x["num",]

#' @rdname as_diffnet
#' @export
`[.diffnet` <- function(x, i, j, drop=FALSE) {

  # Checking range/list of ids or names
  if (inherits(i, "character")) {
    test <- which(!(i %in% x$meta$ids))
  } else if (inherits(i, "numeric") || inherits(i, "integer")) {
    i <- as.integer(i)
    test <- which(!(i %in% 1L:x$meta$n))
  }

  # If any of the elements is not in diffnet, then
  # return error.
  if (length(test))
    stop("Not all -i- in the set of vertices -ids-:\n\t",
         paste0(
           c(head(i[test], 20),
             if (length(test)>20) "..." else NULL),
           collapse=", "
         ))

  # Getting the indexes
  index <- if (inherits(i, "character")) {
    which(x$meta$id %in% i)
  } else {
    which(1:x$meta$n %in% i)
  }

  # Subsetting
  # 1.0: graph and attributes
  ndynattrs <- length(unlist(x$vertex.dyn.attrs)) > 0
  for (i in 1:x$meta$nper) {
    x$graph[[i]] <- x$graph[[i]][index, index, drop=FALSE]
    if (ndynattrs) x$vertex.dyn.attrs[[i]] <- x$vertex.dyn.attrs[[i]][index,,drop=FALSE]
  }

  # 2.0: Matrices
  x$adopt               <- x$adopt[index,,drop=FALSE]
  x$cumadopt            <- x$cumadopt[index,,drop=FALSE]
  x$vertex.static.attrs <- x$vertex.static.attrs[index,,drop=FALSE]
  x$toa                 <- x$toa[index]

  # 3.0: Attrubytes
  x$meta$ids <- x$meta$ids[index]
  x$meta$n   <- length(index)

  return(x)
}

#' \code{diffnet} Arithmetic Operators
#'
#' Addition, substraction and network power of diffnet objects
#'
#' @param x A \code{diffnet} class object.
#' @param y Integer scalar. Power of the network
#' @param valued Logical scalar. When FALSE all non-zero entries of the adjacency
#' matrices are set to one.
#'
#' @details Using binary operators, ease data management process with diffnet.
#'
#' By default the binary operator \code{^} assumes that the graph is valued,
#' hence the power is computed using a weighted edges. Otherwise, if more control
#' is needed, the user can use \code{graph_power} instead.
#'
#' @return A diffnet class object
#'
#' @examples
#' # Computing two-steps away threshold with the Brazilian farmers data --------
#' data(brfarmersDiffNet)
#'
#' expo1 <- threshold(brfarmersDiffNet)
#' expo2 <- threshold(brfarmersDiffNet^2)
#'
#' # Computing correlation
#' cor(expo1,expo2)
#'
#' # Drawing a qqplot
#' qqplot(expo1, expo2)
#'
#' # Working with inverse ------------------------------------------------------
#' brf2_step <- brfarmersDiffNet^2
#' brf2_step <- 1/brf2_step
#'
#' ba <- rdiffnet(10,5, seed.graph="scale-free", rgraph.args=list(m=4))
#'
#' @export
#' @name diffnet-arithmetic
`^.diffnet` <- function(x,y) {
  for (i in 1:x$meta$nper) {
    g <- x$graph[[i]]
    for (p in 1:y)
      x$graph[[i]] <- x$graph[[i]] %*% g
  }
  x
}

#' @rdname diffnet-arithmetic
#' @export
graph_power <- function(x, y, valued=getOption("diffnet.valued", FALSE)) {
  # If no valued
  if (!valued)
    for (i in 1:x$meta$nper)
      x$graph[[i]]@x <- rep(1, length(x$graph[[i]]@x))

  x^y
}

#' @rdname diffnet-arithmetic
#' @export
`/.diffnet` <- function(y, x) {
  for (i in 1:x$meta$nper)
    x$graph[[i]]@x <- y/(x$graph[[i]]@x)

  x
}

# #' @export
# #' @rdname as_diffnet
# subset.diffnet <- function(x, subset, ...) {
#
# }
