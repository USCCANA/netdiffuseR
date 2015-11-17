#' Visualize diffusion process
#'
#' Creates a colored network plot showing the structure of the graph through time
#' (one network plot for each time period)  and the set of adopter and non-adopters
#' in the network.
#'
#' @param graph An array
#' @param cumadopt \eqn{n\times T}{n*T} matrix
#' @param vcols A vector of size 2 with colors
#' @param mode Character. Name of the layout algorithm to implement (see details)
#' @param layout.par Layout parameters (see details)
#' @param mfrow.par Vector of size 2 with number of rows and columns to be passed to \code{\link{par}}
#' @param main Characetr. A title template to be passed to \code{\link{sprintf}}
#' @param ... Further arguments to be passed to gplot
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
#' @return Calculated coordinates (invisible).
#' @export
plot_diffnet <- function(graph, cumadopt, vcols=c("blue","grey"), mode="fruchtermanreingold", layout.par=NULL,
                         mfrow.par=NULL, main="Network in time %d",...) {
  t <- dim(graph)[3]
  n <- dim(graph)[1]

  cols <- matrix(ncol=t, nrow=n)
  cumgraph <- matrix(0, n, n)
  for(i in 1:t) {
    cols[,i] <- ifelse(cumadopt[,i], vcols[1], vcols[2])
    cumgraph <- cumgraph + graph[,,i]
  }

  # Getting the coords
  fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")
  coords <- fun(cumgraph, layout.par)

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
  par(mfrow=mfrow.par)
  for(i in 1:t)  {
    set.seed(curseed)
    sna::gplot(graph[,,i],displaylabels =  TRUE, vertex.col = cols[,i], coord=coords,
               main=sprintf(main, i), ...)
  }
  par(oldpar)

  invisible(coords)

}

#' Creates a \code{diffusionnet} class object
#' @param graph Either an adjacency matrix, an array or an edgelist
as_diffusionnet <- function(graph, ...) {

}

#' @rdname as_diffusionnet
as_diffusionnet.matrix <- function(graph, toa, recode=TRUE, ...) {
  # Figuring out if it is an edgelist
  k <- ncol(graph)
  n <- ncol(graph)
  if ((n!=k) & (k>2)) stop("Invalid -graph-. It should be either an edgelist or a square matrix.")
  else if ((k==2))

  t <- length(unique(toa))

  graph <- array(rep(graph, t), dim=c(n, n, t))
  as_diffusionnet.array(graph, toa, recode)
}

#' @rdname as_diffusionnet
as_diffusionnet.array <- function(graph, toa, recode=TRUE, ...) {

  # Getting times of adoption
  adopmats <- toa_mat(toa, recode)

  list(
    graph=graph,
    adopt.mat=adoptmats$adopt,
    cumadopt.mat=adoptmats$cumadopt,
    toa=toa
  )
}

#' Plots threshold
#' @param graph \eqn{n\times n\times T}{n * n * T} array.
#' @param exposure \eqn{n\times T}{n * T} matrix. Esposure to the innovation obtained from \code{\link{exposure}}
#' @param toa Integer vector of size \eqn{n}. Times of exposure
#' @param toa.recode Logical. TRUE when time recoding must be done
#' @param main Character. Title of the plot
#' @param xlab Character. x-axis label
#' @param ylab Character. y-axis label
#' @param vertex.size Numeric vector of size \eqn{n}. Relative size of the vertices
#' @param vertex.col Either a vector of size \eqn{n} or a scalar indicating colors of the vertices
#' @param vertex.label Character vector of size \eqn{n}. Labels of the vertices
#' @param vertex.lab.pos Integer value to be passed to \code{\link[par]{par through pos}}
#' @param edge.width Numeric. Width of the edges
#' @param edge.col Character. Color of the edges
#' @param arrow.length Numeric value to be passed to \code{\link{arrows}}
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn
#' @param bty See \code{\link{par}}
#' @param ... Additional arguments passed to \code{plot} via \code{\link[sna:gplot]{gplot}}
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
#' plot_threshold(graph, expos, toa, vertex.size = indegree)
#'
#' @export
plot_threshold <- function(graph, exposure, toa, times.recode=TRUE,
                           main="Time of Adoption by Network Threshold", xlab="Time", ylab="Threshold",
                           vertex.size=NULL, vertex.col=rep("lightblue", length(toa)), vertex.label=NULL, vertex.lab.pos=3,
                           edge.width = 2, edge.col = "gray", arrow.length=.20,
                           include.grid = TRUE,
                           bty="n", ...) {
  # Getting basic info
  t <- dim(graph)[3]
  n <- dim(graph)[1]

  # Creating the cumulative graph
  cumgraph <- matrix(0, n, n)
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[,,i]
  }

  # Creating the pos vector
  y <- threshold(exposure, toa, times.recode)

  # Jitter to the xaxis and limits
  jit <- jitter(toa, amount = .25)
  xlim <- range(toa) + c(-1,1)
  yran <- range(y)
  ylim <- yran + (yran[2] - yran[1])*.1*c(-1,1)

  # Plotting
  oldpar <- par(no.readonly = TRUE)
  plot(NULL, xlim=xlim, ylim=ylim, bty=bty, xlab=xlab, ylab=ylab, main=main, ...)

  # Rescaling vertex sizes
  if (length(vertex.size)) {

    # First, for x
    vrange <- range(vertex.size)
    vertex.size <- (vertex.size - vrange[2])/(vrange[1] - vrange[2])/4
    v0 <- which(vertex.size==0)
    vertex.size[v0] <- min(vertex.size[-v0])/2
  }
  else vertex.size <- rep(1/(max(toa)-min(toa))/4, length(toa))

  # Now, for y (it should be different)
  xran <- range(xlim)
  yran <- range(ylim)
  vertex.size.y <- vertex.size *(yran[2]-yran[1])/(xran[2]-xran[1])

  # Drawing arrows
  for(i in 1:n)
    for(j in 1:n) {
      if (!cumgraph[i,j] || toa[i]==toa[j]) next

      # Resizing the edge accordingly to the size of the vertex
      d <- dist(
        rbind( c(jit[i],y[i]),c(jit[j],y[j]) )
        )
      alpha <- acos( (jit[j]-jit[i])/d )

      tox <- jit[j] - cos(alpha)*vertex.size[j]

      # For y, you need to know whether is above or below
      if (y[i] < y[j]) toy <- y[j] - sin(alpha)*vertex.size[j]
      else             toy <- y[j] + sin(alpha)*vertex.size[j]

      # sna::gplot.arrow(jit[i], y[i], tox, toy, width=edge.width)
      arrows(jit[i], y[i], tox, toy, lwd = edge.width, col = edge.col,
             length=arrow.length)
    }

  # Drawing the vertices and its labels
  symbols(jit, y, circle=vertex.size, inches=FALSE, bg=vertex.col, add=TRUE)

  # Positioning labels can be harsh, so we try with this algorithm
  if (!length(vertex.label)) vertex.label <- 1:n
  text(x=jit, y=y+vertex.size.y, labels = vertex.label, pos=vertex.lab.pos)

  if (include.grid) grid()

  par(oldpar)

  invisible(data.frame(toa=toa,threshold=y, jit=jit))

}

