#' Creates a \code{diffnet} class object
#' @param data Either an adjacency matrix, an array or an edgelist
#' @param tao Numeric vector of size \eqn{n}. Times of adoption.
#' @param recode Logical scalar. Passed to \code{\link{toa_mat}}
as_diffnet <- function(data, toa,
                       weights=NULL, times=NULL, undirected=FALSE,
                       self=FALSE, multiple=FALSE,
                       use.incomplete=FALSE,
                       recode=TRUE, ...) {

  # Step 0: Figuring out if it is an edgelist
  d <- dim(data)
  n <- d[1]
  k <- d[2]
  t <- d[3]
  if ((n!=k) & (k>2)) stop("Invalid -data-. It should be either an edgelist or a square matrix.")
  else if ((k==2))

  # Step 1: Creating the graph, first we need to see the time length
  trange <- range(toa)
  t <- trange[2]-trange[1]

  if      (type == "adjmat")   graph <- array(rep(graph, t), dim=c(n, n, t))
  else if (type == "edgelist") graph <- edgelist_to_adjmat(
    data, weights, times, t, undirected=undirected, self=self, multiple=multiple)
  else {
    # If its already an array, we better check for the names!
    graph <- data
  }

  return(data)
}

#' Plot the diffusion process
#'
#' Creates a colored network plot showing the structure of the graph through time
#' (one network plot for each time period)  and the set of adopter and non-adopters
#' in the network.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param cumadopt \eqn{n\times T}{n*T} matrix
#' @param displaylabels Logical. When TRUE vertex labels are displayed (see \code{\link[sna:gplot]{gplot}})
#' @param vertex.col A character vector of size 2 with colors
#' @param vertex.cex Numeric vector of size \eqn{n}. Size of the vertices
#' @param label Character vector of size \eqn{n}. If no provided, rownames of
#' the graph are used.
#' @param edge.col Character. Color of the edge
#' @param mode Character. Name of the layout algorithm to implement (see details)
#' @param layout.par Layout parameters (see details)
#' @param mfrow.par Vector of size 2 with number of rows and columns to be passed to \code{\link{par}}
#' @param main Characetr. A title template to be passed to \code{\link{sprintf}}
#' @param mai Numeric vector of size 4. To be passed to \code{\link{par}}
#' @param mar Numeric vector of size 4. To be passed to \code{\link{par}}
#' @param ... Further arguments to be passed to \code{\link[sna:gplot]{gplot}}
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
#' @examples
#' #' # Generating a random graph
#' set.seed(1234)
#' n <- 6
#' nper <- 5
#' graph <- rand_graph(n,nper, p=.3, undirected = FALSE)
#' toa <- sample(2000:(2000+nper-1), n, TRUE)
#' adopt <- toa_mat(toa)
#'
#' plot_diffnet(graph, adopt$cumadopt)
#' @return Calculated coordinates (invisible).
#' @family visualizations
#' @keywords hplot
#' @export
plot_diffnet <- function(graph, ...) UseMethod("plot_diffnet")

#' @export
#' @rdname plot_diffnet
plot_diffnet.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) graph[,,x])
  names(graph) <- dn
  plot_diffnet.list(graph, ...)
}

#' @export
#' @rdname plot_diffnet
plot_diffnet.list <- function(graph, cumadopt,
                         displaylabels=FALSE,
                         undirected=TRUE,
                         vertex.col=c("blue","grey"),
                         vertex.cex=NA,
                         label=rownames(graph[[1]]),
                         edge.col="gray",
                         mode="fruchtermanreingold", layout.par=NULL,
                         mfrow.par=NULL, main="Network in time %d",
                         mai=c(0,0,1,0),
                         mar=rep(1,4) + 0.1, ...) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  cumgraph <- Matrix::Matrix(0, n, n, sparse=TRUE)
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[[i]]
  }

  # Getting the coords
  fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")

  # In order to be using SNA functions, we need to coerse the graph
  # into an object from SparseM
  coords <- fun(methods::as(cumgraph, "matrix.csc"), layout.par)

  # Adjusting vertex sizes
  if ((length(vertex.cex)==1) && (n > 1) && is.na(vertex.cex))
    vertex.cex <- (max(coords[,1]) - min(coords[,1]))/10*3

  if ( (length(vertex.cex)==1) && (n > 1) )
    vertex.cex <- rep(vertex.cex,n)

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

  # Plotting
  curseed <- .Random.seed
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=mfrow.par, mai=mai, mar=mar)

  times <- as.integer(names(graph))

  for(i in 1:t)  {
    cols <- ifelse(cumadopt[,i], vertex.col[1], vertex.col[2])
    set.seed(curseed)
    # cgraph <- sna::as.sociomatrix.sna(adjmat_to_edgelist(graph[[i]], undirected))
    sna::gplot(methods::as(graph[[i]], "matrix.csc"),
               displaylabels = displaylabels, vertex.col = cols, coord=coords,
               edge.col = edge.col,vertex.cex = vertex.cex, label=label,
               main=sprintf(main, times[i]), ...)
  }

  invisible(coords)

}

#' Threshold level through time
#'
#' Draws a graph where the coordinates are given by time of adoption, x-axis,
#' and threshold level, y-axis.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param exposure \eqn{n\times T}{n * T} matrix. Esposure to the innovation obtained from \code{\link{exposure}}
#' @param toa Integer vector of size \eqn{n}. Times of Adoption
#' @param times.recode Logical scalar. TRUE when time recoding must be done.
#' @param main Character. Title of the plot
#' @param xlab Character. x-axis label
#' @param ylab Character. y-axis label
#' @param vertex.cex Numeric vector of size \eqn{n}. Relative size of the vertices
#' @param vertex.col Either a vector of size \eqn{n} or a scalar indicating colors of the vertices
#' @param vertex.label Character vector of size \eqn{n}. Labels of the vertices
#' @param vertex.lab.pos Integer value to be passed to \code{\link{text}} via \code{pos}
#' @param edge.width Numeric. Width of the edges
#' @param edge.col Character. Color of the edges
#' @param arrow.length Numeric value to be passed to \code{\link{arrows}}
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn
#' @param bty See \code{\link{par}}
#' @param ... Additional arguments passed to \code{plot} via \code{\link[sna:gplot]{gplot}}
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
#' graph <- rand_graph(n,nper, p=.3, undirected = FALSE)
#' toa <- sample(2000:(2000+nper-1), n, TRUE)
#' adopt <- toa_mat(toa)
#'
#' # Computing exposure
#' expos <- exposure(graph, adopt$cumadopt, undirected = FALSE)
#'
#' plot_threshold(graph, expos, toa)
#'
#' # Calculating degree (for sizing the vertices)
#' indegree <- netdiffuseR::degree(graph, cmode="indegree")
#' indegree <- apply(indegree, 1, mean)
#' plot_threshold(graph, expos, toa, vertex.cex = indegree)
#'
#' @export
plot_threshold <- function(graph, ...) UseMethod("plot_threshold")

#' @export
#' @rdname plot_threshold
plot_threshold.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) graph[,,x])
  names(graph) <- dn
  plot_threshold.list(graph, ...)
}

#' @export
#' @rdname plot_threshold
plot_threshold.list <- function(graph, exposure, toa, times.recode=TRUE, undirected=TRUE,
                           main="Time of Adoption by Network Threshold", xlab="Time", ylab="Threshold",
                           vertex.cex=NA,
                           vertex.col="blue", vertex.label=NULL, vertex.lab.pos=3,
                           edge.width = 2, edge.col = "gray", arrow.length=.20,
                           include.grid = TRUE,
                           bty="n", ...) {
  # Step 0: Getting basic info
  t <- length(graph)
  n <- nrow(graph[[1]])

  # Step 1: Creating the cumulative graph
  cumgraph <- Matrix::Matrix(0, n, n, sparse=TRUE)
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[[i]]
  }

  # Creating the pos vector
  y <- threshold(exposure, toa, times.recode)

  # Jitter to the xaxis and limits
  jit <- jitter(toa, amount = .25)
  xran <- range(toa)
  xlim <- xran + c(-1,1)
  yran <- range(y)
  ylim <- yran + (yran[2] - yran[1])*.1*c(-1,1)

  # Step 2: Checking colors and sizes

  # Adjusting size of the vertices
  if ((length(vertex.cex)==1) && (n > 1) && is.na(vertex.cex))
    vertex.cex <- (xran[2] - xran[1])/10/3

  if ( (length(vertex.cex)==1) && (n > 1) )
    vertex.cex <- rep(vertex.cex,n)

  if ( (length(vertex.col)==1) && (n > 1) )
    vertex.col <- rep(vertex.col,n)

  # Plotting
  oldpar <- par(no.readonly = TRUE)
  plot(NULL, xlim=xlim, ylim=ylim, bty=bty, xlab=xlab, ylab=ylab, main=main, ...)

  # Should there be a grid??
  if (include.grid) grid()

  # Now, for y (it should be different)
  xran <- range(xlim)
  yran <- range(ylim)
  vertex.cex.y <- vertex.cex *(yran[2]-yran[1])/(xran[2]-xran[1])

  # Drawing arrows, first we calculate the coordinates of the edges, for this we
  # use the function edges_coords. This considers aspect ratio of the plot.
  edges <- netdiffuseR::edges_coords(cumgraph, toa, jit, y, vertex.cex, undirected)
  edges <- as.data.frame(edges)

  with(edges, arrows(x0, y0, x1, y1, lwd = edge.width, col = edge.col,
                                length=arrow.length))

  # Drawing the vertices and its labels
  symbols(jit, y, circles=vertex.cex, inches=FALSE, bg=vertex.col, add=TRUE)

  # Positioning labels can be harsh, so we try with this algorithm
  if (!length(vertex.label)) vertex.label <- 1:n
  text(x=jit, y=y+vertex.cex.y, labels = vertex.label, pos=vertex.lab.pos)

  par(oldpar)

  invisible(data.frame(toa=toa,threshold=y, jit=jit))

}

#' Plot distribution of infect/suscep
#'
#' After calculating infectiousness and susceptibility of each individual on the
#' network, it creates an \code{nlevels} by \code{nlevels} matrix indicating the
#' number of individuals that lie within each cell, and draws a heatmap.
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param toa Times of adoption
#' @param normalize Logical. TRUE
#' @param K Integer. Number of time periods to consider
#' @param r Double. Rate
#' @param expdiscount Logical.
#' @param bins Integer scalar. Size of the grid (\eqn{n}).
#' @param nlevels Integer. Number of levels to plot (see \code{\link{filled.contour}}).
#' @param logscale Logical. When TRUE the axis of the plot will be presented in log-scale.
#' @param main Character. Title of the graph.
#' @param xlab Character. Title of the x-axis.
#' @param ylab Character. Title of the y-axis.
#' @param sub Character. Subtitle of the graph.
#' @param color.palette a color palette function to be used to assign colors in the plot (see \code{\link{filled.contour}}).
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn.
#' @param ... Additional parameters to be passed to \code{\link{filled.contour}}
#' @details
#'
#' This plotting function was inspired by Aral, S., & Walker, D. (2012).
#'
#' @return A list with three elements:
#' \item{infect}{A numeric vector of size \eqn{n} with infectiousness levels}
#' \item{suscep}{A numeric vector of size \eqn{n} with susceptibility levels}
#' \item{coords}{A list containing the class marks and counts used to draw the
#' plot via \code{\link{filled.contour}} (see \code{\link{grid_distribution}})}
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
#' n <- 500
#' nper <- 30
#' graph <- rand_graph(n,nper, p=.2, undirected = FALSE)
#' toa <- sample(1:(1+nper-1), n, TRUE)
#'
#' # Visualizing distribution of suscep/infect
#' out <- plot_infectsuscep(graph, toa, K=3, logscale = TRUE)
plot_infectsuscep <- function(graph, ...) UseMethod("plot_infectsuscept")

#' @export
#' @rdname plot_infectsuscep
plot_infectsuscep.array <- function(graph, ...) {
  dn <- dimnames(graph)[[3]]
  graph <- lapply(1:dim(graph)[3], function(x) graph[,,x])
  names(graph) <- dn
  plot_infectsuscep.list(graph, ...)
}

#' @export
#' @rdname plot_infectsuscep
plot_infectsuscep.list <- function(graph, toa, normalize=TRUE,
                              K=1L, r=0.5, expdiscount=FALSE,
                              bins=50,nlevels=round(bins/2),
                              logscale=TRUE,
                              main="Distribution of Infectiousness and\nSusceptibility",
                              xlab="Infectiousness of ego",
                              ylab="Susceptibility of ego",
                              sub=ifelse(logscale, "(in log-scale)", NA),
                              color.palette=function(n) grey(n:1/n),
                              include.grid=TRUE,
                              ...) {
  # Computing infect and suscept
  infect <- infection(graph, toa, normalize, K, r, expdiscount)
  suscep <- susceptibility(graph, toa, normalize, K, r, expdiscount)

  # Performing classification (linear)
  if (logscale) {
    infect<-log(infect); infect[which(!is.finite(infect))] <- 0
    suscep<-log(suscep); suscep[which(!is.finite(suscep))] <- 0
  }
  coords <- netdiffuseR::grid_distribution(x=infect, y=suscep, bins)

  # Nice plot
  n <- sum(coords$z)
  with(coords, filled.contour(
    x,y,z/n, bty="n", main=main, xlab=xlab, ylab=ylab, sub=sub, color.palette =color.palette,
    plot.axes={
      axis(1);axis(2)
      if (include.grid) grid()
    }, nlevels=nlevels, ...))
  # if (include.grid) grid()

  invisible(list(infect=infect, suscept=suscep, coords=coords))
}
