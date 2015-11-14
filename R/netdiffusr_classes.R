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
