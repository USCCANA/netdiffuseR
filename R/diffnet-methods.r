#' S3 plotting method for diffnet objects.
#'
#' @param x An object of class \code{\link[=diffnet-class]{diffnet}}
#' @param t Integer scalar indicating the time slice to plot.
#' @param displaylabels Logical scalar. When TRUE, \code{plot} shows vertex labels.
#' @param vertex.col Character scalar/vector. Color of the vertices.
#' @param vertex.cex Numeric scalar/vector. Size of the vertices.
#' @param edge.col Character scalar/vector. Color of the edges.
#' @param mode Character scalar. In the case of \code{plot}, passed to
#' \code{\link[sna:gplot]{gplot}}.
#' @param layout.par Layout parameters (see details).
#' @param gmode Character scalar. Passed to \code{\link[sna:gplot]{sna::gplot}.}
#' @param main Character. A title template to be passed to sprintf.
#' @param ... Further arguments passed to \code{\link[sna:gplot]{gplot}}.
#' @param y Ignored.
#' @export
#'
#' @family diffnet methods
#' @details
#' Plotting is done via the function \code{\link[sna:gplot]{gplot}},
#' and its layout via \code{\link[sna:gplot.layout]{gplot.layout}}, both from
#' the (\pkg{sna}) package.
#'
#' \code{vertex.cex} can either be a numeric scalar, a numeric vector or a character
#' scalar taking any of the following values \code{"degree"}, \code{"indegree"}, or
#' \code{"outdegree"}. The later will be passed to \code{\link{dgr}} to calculate
#' degree of the selected slice and will be normalized as
#'
#' \deqn{
#'   vertex.cex = d/[max(d) - min(d)]\times 2 + .5}{%
#'   vertex.cex = d/[max(d) - min(d)]* 2 + .5 %
#'  }
#'
#' where \code{d=sqrt(dgr(graph))}.
#'
#' @return A matrix with the coordinates of the vertices.
#' @author George G. Vega Yon
#' @examples
#'
#' data(medInnovationsDiffNet)
#' plot(medInnovationsDiffNet)
#'
#'
plot.diffnet <- function(
  x,y=NULL, t=1, displaylabels = FALSE,
  vertex.col = c(adopt="blue", noadopt="grey"),
  gmode=ifelse(x$meta$undirected, "graph", "digraph"),
  vertex.cex = "degree", edge.col = "gray", mode = "fruchtermanreingold",
  layout.par = NULL, main = "Diffusion network in time %d", ...) {

  # Listing arguments
  args <- list(...)

  # Checking that the time period is actually within
  if (!(t %in% 1:x$meta$nper))
    stop("-t- must be an integer within 1 and ",x$meta$nper,".")

  # Extracting the graph to be plotted
  graph <- methods::as(x$graph[[t]], "matrix.csc")

  # Setting the colors
  cols <- with(x, ifelse(cumadopt[,t], vertex.col[1], vertex.col[2]))

  # Calcularing layout

   if ("coord" %in% names(args)) {
     coords <- args[["coord"]]
     args[["coord"]] <- NULL
   } else {
     fun <- getFromNamespace(paste0("gplot.layout.",mode), "sna")
     coords <- fun(graph, layout.par)
   }

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character"))
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(x$graph[[t]], cmode=vertex.cex, undirected = x$meta$undirected)
      vertex.cex <- sqrt(vertex.cex)
      r <- range(vertex.cex)

      # If all the vertices have the same degree
      vertex.cex <- (vertex.cex - r[1]+ .1)/(r[2] - r[1] + .1)*2
    } else {
      stop("Invalid -vertex.cex-")
    }

  do.call(sna::gplot, c(list(dat=graph, displaylabels=displaylabels, vertex.col=cols,
             coord=coords, edge.col=edge.col, label=x$meta$ids,
             main=sprintf(main, x$meta$pers[t]), vertex.cex=vertex.cex, gmode=gmode),
             args))

  invisible(coords)
}

#' @export
#' @rdname diffnet-class
print.diffnet <- function(x, ...) {
  with(x, {
    # Getting attrs
    vsa <- paste0(colnames(vertex.static.attrs), collapse=", ")
    if (nchar(vsa) > 50) vsa <- paste0(strtrim(vsa, 50),"...")
    else if (!nchar(vsa)) vsa <- '-'
    nsa <-ncol(vertex.static.attrs)
    if (nsa) vsa <- paste0(vsa," (",nsa, ")")

    vda <- paste0(colnames(vertex.dyn.attrs[[1]]), collapse=", ")
    if (nchar(vda) > 50) vda <- paste0(strtrim(vda, 50),"...")
    else if (!nchar(vda)) vda <- '-'
    nda <- ncol(vertex.dyn.attrs[[1]])
    if (nda) vda <- paste0(vda," (",nda, ")")

    # Getting nodes labels
    nodesl <- paste0(meta$n," (",
                     paste(head(meta$ids, 8), collapse=", "),
                     ifelse(meta$n>8, ", ...", "") ,")")

    cat(
    "Dynamic network of class -diffnet-",
    paste(" Name               :", meta$name),
    paste(" Behavior           :", meta$behavior),
    paste(" # of nodes         :", nodesl ),
    paste(" # of time periods  :", meta$nper, sprintf("(%d - %d)", meta$pers[1], meta$pers[meta$nper])),
    paste(" Type               :", ifelse(meta$undirected, "undirected", "directed")),
    paste(" Final prevalence   :",
          formatC(sum(cumadopt[,meta$nper])/meta$n, digits = 2, format="f")
          ),
    paste(" Static attributes  :", vsa),
    paste(" Dynamic attributes :", vda),
    sep="\n"
    )
  })
  invisible(x)
}

#' Summary of diffnet objects
#'
#' @export
#' @param object An object of class \code{\link[=as_diffnet]{diffnet}}.
#' @param slices Either an integer or character vector. While integer vectors are used as
#' indexes, character vectors are used jointly with the time period labels.
#' @param valued Logical scalar. When \code{TRUE} weights will be considered.
#' Otherwise non-zero values will be replaced by ones.
#' @param no.print Logical scalar. When TRUE suppress screen messages.
#' @param skip.moran Logical scalar. When TRUE Moran's I is not reported (see details).
#' @param ... Further arguments to be passed to \code{\link{approx_geodesic}}.
#' @details
#' Moran's I is calculated over the
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
#' Moran's I is not reported. This can be useful for both: reducing computing
#' time and saving memory as geodesic distance matrix can become large. Since
#' version \code{1.18.0}, geodesic matrices are approximated using \code{approx_geodesic}
#' which, as a difference from \code{\link[sna:geodist]{geodist}} from the
#' \pkg{sna} package, and \code{\link[igraph:distances]{distances}} from the
#' \pkg{igraph} package returns a matrix of class \code{dgCMatrix} (more
#' details in \code{\link{approx_geodesic}}).
#'
#' @return A data frame with the following columns:
#' \item{adopt}{Integer. Number of adopters at each time point.}
#' \item{cum_adopt}{Integer. Number of cumulative adopters at each time point.}
#' \item{cum_adopt_pcent}{Numeric. Proportion of comulative adopters at each time point.}
#' \item{hazard}{Numeric. Hazard rate at each time point.}
#' \item{density}{Numeric. Density of the network at each time point.}
#' \item{moran_obs}{Numeric. Observed Moran's I.}
#' \item{moran_exp}{Numeric. Expected Moran's I.}
#' \item{moran_sd}{Numeric. Standard error of Moran's I under the null.}
#' \item{moran_pval}{Numeric. P-value for the observed Moran's I.}
#' @author George G. Vega Yon
#'
#' @examples
#' data(medInnovationsDiffNet)
#' summary(medInnovationsDiffNet)
#'
#' @family diffnet methods
#'
summary.diffnet <- function(
  object,
  slices     = NULL,
  no.print   = FALSE,
  skip.moran = FALSE,
  valued     = getOption("diffnet.valued",FALSE),
  ...) {
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
    nlinks(x)/nnodes(x)/(nnodes(x)-1)
    # nelements <- length(x@x)
    # x <-nelements/(meta$n * (meta$n-1))
  }))

  # Computing moran's I
  if (!skip.moran) {

    m <- matrix(NA, nrow=length(slices), ncol=4,
                dimnames = list(NULL, c("moran_obs", "moran_exp", "moran_sd", "moran_pval")))

    for (i in 1:length(slices)) {
      # Computing distances
      g <- approx_geodesic(object$graph[[slices[i]]], ...)

      # Inverting it (only the diagonal may have 0)
      g@x <- 1/g@x

      m[i,] <- unlist(moran(object$cumadopt[,slices[i]], g))
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

  if (!skip.moran) {
    out <- cbind(out, m)
  }

  if (no.print) return(out)

  # Function to print data.frames differently
  header <- c(" Period "," Adopters "," Cum Adopt. (%) ",
              " Hazard Rate "," Density ",
              if (!skip.moran) c(" Moran's I (sd) ") else NULL
              )

  slen   <- nchar(header)
  hline  <- paste(sapply(sapply(slen, rep.int, x="-"), paste0, collapse=""),
                  collapse=" ")
  rule   <- paste0(rep("-", sum(slen) + length(slen) - 1), collapse="")

  # Quick Formatting function
  qf <- function(x, digits=2) sprintf(paste0("%.",digits,"f"), x)

  cat("Diffusion network summary statistics\n",
      "Name     : ", meta$name, "\n",
      "Behavior : ", meta$behavior, "\n",
      rule,"\n",sep="")
  cat(header,"\n")
  cat(hline, "\n")
  for (i in 1:nrow(out)) {
    cat(sprintf(
      paste0("%",slen,"s", collapse=" "),
      qf(meta$pers[slices[i]],0), qf(out[i,1],0),
      sprintf("%s (%s)",
        qf(out$cum_adopt[i],0),
        qf(out$cum_adopt_pcent[i])
        ),
      ifelse(i==1, "-",qf(out$hazard[i])), qf(out$density[i]),
      if (!skip.moran) {
        if (is.nan(out$moran_sd[i]))
          " - "
        else
          sprintf("%s (%s) %-3s",
                  qf(out$moran_obs[i]),
                  qf(out$moran_sd[i]),
                  ifelse(out$moran_pval[i] <= .01, "***",
                         ifelse(out$moran_pval[i] <= .05, "**",
                                ifelse(out$moran_pval[i] <= .10, "*", ""
                              )))
                )
      } else ""
    ), "\n")
  }


  # print(out, digits=2)

  cat(
    rule,
    paste(" Left censoring  :", sprintf("%3.2f (%d)", lc/meta$n, lc)),
    paste(" Right centoring :", sprintf("%3.2f (%d)", rc/meta$n, rc)),
    paste(" # of nodes      :", sprintf("%d",meta$n)),
    "\n Moran's I was computed on contemporaneous autocorrelation using 1/geodesic",
    " values. Significane levels  *** <= .01, ** <= .05, * <= .1.",
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
#' @templateVar dynamic TRUE
#' @templateVar undirected TRUE
#' @template graph_template
#' @param cumadopt \eqn{n\times T}{n*T} matrix.
#' @param slices Integer vector. Indicates what slices to plot. By default all are plotted.
#' @param vertex.col A character vector of size 3 with colors names.
#' @param vertex.shape A character vector of size 3 with shape names.
#' @param vertex.cex Numeric vector of size \eqn{n}. Size of the vertices.
#' @param label Character vector of size \eqn{n}. If no provided, rownames of
#' the graph are used.
#' @param edge.col Character scalar/vector. Color of the edge.
#' @param mode Character scalar. Name of the layout algorithm to implement (see details).
#' @param layout.par Layout parameters (see details).
#' @param mfrow.par Vector of size 2 with number of rows and columns to be passed to \code{\link{par}.}
#' @param main Character scalar. A title template to be passed to \code{\link{sprintf}.}
#' @param gmode Character scalar. See \code{\link[sna:gplot]{gplot}.}
#' @param vertex.frame.color Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}
#' @param edge.arrow.size Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}
#' @param intra.space Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}
#' @param key.height Numeric scalar. Sets the proportion of the plot (y-axis) that the key uses.
#' @param rescale.fun A function to rescale vertex size. By defult it is set to be \code{\link{rescale_vertex_igraph}}
#' @param ... Further arguments to be passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param coords Numeric matrix of size \eqn{n\times 2}{n * 2} with vertices coordinates.
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
#' @author George G. Vega Yon
plot_diffnet <- function(
  graph, cumadopt,
  slices=NULL,
  undirected=TRUE,
  vertex.col=c("white", "red", "blue"),
  vertex.shape=c("square", "circle", "circle"),
  vertex.cex="degree",
  label=NA,
  edge.col="gray",
  mode="fruchtermanreingold", layout.par=NULL,
  mfrow.par=NULL, main="Network in period %d",
  gmode=ifelse(undirected, "graph", "digraph"),
  lgd = list(x="bottom", legend=c("Non-adopters", "New adopters","Adopters"), pch=21,
             bty="n", cex=1.2, horiz=TRUE), coords=NULL,
  vertex.frame.color="gray",
  edge.arrow.size=.25,
  intra.space=c(.15,.15),
  key.height = 0.1,
  rescale.fun = function(x) rescale_vertex_igraph(x, adjust = 100),
  ...
) {
  switch (class(graph),
    array = plot_diffnet.array(
      graph, cumadopt, slices, undirected, vertex.col, vertex.shape, vertex.cex, label,
      edge.col, mode, layout.par, mfrow.par, main, gmode, lgd, coords,
      vertex.frame.color, edge.arrow.size, intra.space, key.height, rescale.fun,...),
    list = plot_diffnet.list(
      graph, cumadopt, slices, undirected, vertex.col, vertex.shape, vertex.cex, label,
      edge.col, mode, layout.par, mfrow.par, main, gmode, lgd, coords,
      vertex.frame.color, edge.arrow.size, intra.space, key.height, rescale.fun,...),
    diffnet = plot_diffnet.list(
      graph$graph, graph$cumadopt, slices, graph$meta$undirected,
      vertex.col, vertex.shape, vertex.cex, label,
      edge.col, mode, layout.par, mfrow.par, main, gmode, lgd, coords,
      vertex.frame.color, edge.arrow.size, intra.space, key.height, rescale.fun,...),
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_diffnet
plot_diffnet.array <- function(graph, ...) {
  graph <- apply(graph, 3, methods::as, Class="dgCMatrix")
  plot_diffnet.list(graph, ...)
}

# @export
# @rdname plot_diffnet
plot_diffnet.list <- function(graph, cumadopt, slices,
                         undirected=TRUE,
                         vertex.col=c("white", "red", "blue"),
                         vertex.shape=c("square", "circle", "circle"),
                         vertex.cex="degree",
                         label=NA,
                         edge.col="gray",
                         mode="fruchtermanreingold", layout.par=NULL,
                         mfrow.par=NULL, main="Network in period %d",
                         gmode=ifelse(undirected, "graph", "digraph"),
                         lgd = list(x="bottom", legend=c("Non adopters", "New adopters","Adopters"), pch=21,
                                    bty="n", cex=1.2, horiz=TRUE),
                         coords=NULL,
                         vertex.frame.color="gray",
                         edge.arrow.size=.25,
                         intra.space=c(.15,.15),
                         key.height = 0.1,
                         rescale.fun = function(x) rescale_vertex_igraph(x, adjust = 100),
                         ...) {

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
  if (!length(coords)) coords <- fun(methods::as(cumgraph, "matrix.csc"), layout.par)

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character"))
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(cumgraph, cmode = vertex.cex, undirected=undirected)
    } else {
      stop("Invalid -vertex.cex-")
    }

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


  # 1. Computing dims ----------------------------------------------------------
  intra.space <- intra.space + 1
  xrange <- range(coords[,1])
  xlim   <- xrange + c(0, (xrange[2]-xrange[1])*intra.space[1]*(mfrow.par[2] - 1))

  yrange <- range(coords[,2])
  ylim   <- yrange - c((yrange[2]-yrange[1])*intra.space[2]*(mfrow.par[1] - 1), 0)

  # Adjustems depending on the number of slice
  irow <- t(matrix(0:(mfrow.par[1]-1), nrow=mfrow.par[1], ncol=mfrow.par[2], byrow = FALSE))
  icol <- t(matrix(0:(mfrow.par[2]-1), nrow=mfrow.par[1], ncol=mfrow.par[2], byrow = TRUE))

  # 2. Set up frame
  plot.new()
  ylim <- grDevices::extendrange(ylim, ylim)
  plot.window(
    grDevices::extendrange(xlim),
    ylim - c((ylim[2]-ylim[1])*key.height, 0))

  # 3. Plotting
  times <- as.integer(names(graph))
  for (i in 1:t) {
    # Colors, new adopters are painted differently
    cols <- ifelse(!cumadopt[,i], vertex.col[1],
                   ifelse(!cumadopt[,i-(i!=1)] | rep(i,n) == 1, vertex.col[2], vertex.col[3]))

    # Shapes
    shapes <- ifelse(!cumadopt[,i], vertex.shape[1],
                     ifelse(!cumadopt[,i-(i!=1)] | rep(i,n) == 1, vertex.shape[2], vertex.shape[3]))

    # Computing coords
    coords_adjs <- cbind(
      coords[,1] + icol[i]*intra.space[1]*(xrange[2]-xrange[1]),
      coords[,2] - irow[i]*intra.space[2]*(yrange[2]-yrange[1]))

    # Title
    if (length(main)) {
      # Creating the title and computing how much space should be
      itit <- sprintf(main, slices[i])
      nlines <- length(strsplit(main, "\n")[[1]])
      tspace <- par("cxy")[2]*(nlines + .5)

      # Adjusting coords
      ystart <- range(coords_adjs[,2])
      coords_adjs[,2] <- (coords_adjs[,2] - ystart[1])/(ystart[2] - ystart[1])*(
        ystart[2] - ystart[1] - tspace) + ystart[1]

      # Including text
      xstart <- range(coords_adjs[,1])
      text(
        x = (xstart[2] + xstart[1])/2,
        y = ystart[2] - tspace,
        labels = itit,
        pos=3
      )

    }

    # Coercing into the right method
    if (!inherits(graph[[i]], "dgCMatrix")) g <- methods::as(graph[[i]], "dgCMatrix")
    else g <- graph[[i]]

    # Checking dimnames
    if (!length(unlist(dimnames(g), recursive = TRUE)))
      dimnames(g) <- list(1:nnodes(g), 1:nnodes(g))

    # Creating igraph object
    ig  <- igraph::graph_from_adjacency_matrix(g)
    ig  <- igraph::permute(ig, match(igraph::V(ig)$name, nodes(g)))


    # Plotting
    igraph::plot.igraph(ig,
                        vertex.color = cols,
                        layout = coords_adjs,
                        edge.color = edge.col,
                        vertex.size = rescale.fun(vertex.cex),
                        vertex.label=label,
                        add=TRUE, rescale=FALSE,
                        edge.arrow.size=edge.arrow.size,
                        vertex.frame.color = vertex.frame.color,
                        vertex.shape=shapes,
                        ...)
  }

  # Legend
  do.call(legend, c(lgd, list(pt.bg=vertex.col)))

  invisible(coords)

}

#' Threshold levels through time
#'
#' Draws a graph where the coordinates are given by time of adoption, x-axis,
#' and threshold level, y-axis.
#'
#' @templateVar dynamic TRUE
#' @templateVar toa TRUE
#' @templateVar undirected TRUE
#' @template graph_template
#' @param expo \eqn{n\times T}{n * T} matrix. Esposure to the innovation obtained from \code{\link{exposure}}
#' @param t0 Integer scalar. Passed to \code{\link{threshold}}.
#' @param include_censored Logical scalar. Passed to \code{\link{threshold}}.
#' @param attrs Passed to \code{\link{exposure}} (via threshold).
#' @param no.contemporary Logical scalar. When TRUE, edges for vertices with the same
#' \code{toa} won't be plotted.
#' @param main Character scalar. Title of the plot.
#' @param xlab Character scalar. x-axis label.
#' @param ylab Character scalar. y-axis label.
#' @param vertex.cex Numeric vector of size \eqn{n}. Relative size of the vertices.
#' @param vertex.col Either a vector of size \eqn{n} or a scalar indicating colors of the vertices.
#' @param vertex.label Character vector of size \eqn{n}. Labels of the vertices.
#' @param vertex.lab.pos Integer value to be passed to \code{\link{text}} via \code{pos}.
#' @param vertex.lab.cex Either a numeric scalar or vector of size \eqn{n}. Passed to \code{text}.
#' @param vertex.lab.adj Passed to \code{\link{text}}.
#' @param vertex.lab.col Passed to \code{\link{text}}.
#' @param jitter.amount Numeric vector of size 2 (for x and y) passed to \code{\link{jitter}}.
#' @param jitter.factor Numeric vector of size 2 (for x and y) passed to \code{\link{jitter}}.
#' @param vertex.bcol Either a vector of size \eqn{n} or a scalar indicating colors of vertices' borders.
#' @param vertex.sides Either a vector of size \eqn{n} or a scalar indicating the
#' number of sides of each vertex (see details).
#' @param vertex.rot Either a vector of size \eqn{n} or a scalar indicating the
#' rotation in radians of each vertex (see details).
#' @param edge.width Numeric. Width of the edges.
#' @param edge.col Character. Color of the edges.
#' @param arrow.length Numeric value to be passed to \code{\link{arrows}}.
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn.
#' @param bty See \code{\link{par}}.
#' @param xlim Passed to \code{\link{plot}}.
#' @param ylim Passed to \code{\link{plot}}.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @family visualizations
#' @seealso Use \code{\link{threshold}} to retrieve the corresponding threshold
#' obtained returned by \code{\link{exposure}}.
#' @keywords hplot
#'
#' @details When \code{vertex.label=NULL} the function uses vertices ids as labels.
#' By default \code{vertex.label=""} plots no labels.
#'
#' Vertices are drawn using an internal function for generating polygons.
#' Polygons are inscribed in a circle of radius \code{vertex.cex}, and can be
#' rotated using \code{vertex.rot}. The number of sides of each polygon
#' is set via \code{vertex.sides}.
#'
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
#' @author George G. Vega Yon
plot_threshold <- function(
  graph, expo, toa,
  include_censored=FALSE,
  t0=min(toa, na.rm = TRUE), attrs=NULL,
  undirected=getOption("diffnet.undirected"), no.contemporary=TRUE,
  main="Time of Adoption by Network Threshold", xlab="Time", ylab="Threshold",
  vertex.cex="degree", vertex.col=rgb(.3,.3,.8,.5),
  vertex.label="", vertex.lab.pos=NULL,  vertex.lab.cex=1,
  vertex.lab.adj = c(.5,.5), vertex.lab.col=rgb(.3,.3,.8,.9),
  vertex.sides = 40L, vertex.rot = 0,
  edge.width = 2, edge.col = rgb(.6,.6,.6,.1), arrow.length=.20,
  include.grid = TRUE, bty="n",
  vertex.bcol=vertex.col, jitter.factor=c(1,0), jitter.amount=c(.25,0),
  xlim=NULL, ylim=NULL, ...
) {

  if (missing(expo))
    if (!inherits(graph, "diffnet")) {
      stop("-expo- should be provided when -graph- is not of class 'diffnet'")
    } else {
      expo <- exposure(graph)
    }

  switch (class(graph),
    array = plot_threshold.array(
      graph, expo, toa, include_censored, t0, attrs, undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label,
      vertex.lab.pos, vertex.lab.cex, vertex.lab.adj,vertex.lab.col,
      vertex.sides, vertex.rot,
      edge.width, edge.col,
      arrow.length, include.grid, bty, vertex.bcol, jitter.factor, jitter.amount,
      xlim, ylim, ...),
    list = plot_threshold.list(
      graph, expo, toa, include_censored, t0, attrs, undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label, vertex.lab.pos, vertex.lab.cex,
      vertex.lab.adj, vertex.lab.col,
      vertex.sides, vertex.rot,
      edge.width, edge.col,
      arrow.length, include.grid, bty, vertex.bcol,jitter.factor, jitter.amount,
      xlim, ylim, ...),
    diffnet = {
      # If graph is diffnet, then we should do something different (because the
      # first toa may not be the firts one as toa may be stacked to the right.
      # see ?as_diffnet)
      # graph$toa <- graph$toa - min(graph$meta$pers) + 1L

      plot_threshold.list(
      graph$graph, expo,
      graph$toa, include_censored, t0=graph$meta$pers[1], attrs, graph$meta$undirected, no.contemporary, main, xlab, ylab,
      vertex.cex, vertex.col, vertex.label, vertex.lab.pos, vertex.lab.cex,
      vertex.lab.adj,vertex.lab.col,
      vertex.sides, vertex.rot,
      edge.width, edge.col,
      arrow.length, include.grid, bty, vertex.bcol, jitter.factor, jitter.amount,
      xlim, ylim, ...)
      },
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_threshold
plot_threshold.array <- function(graph, ...) {
  graph <- apply(graph, 3, methods::as, Class="dgCMatrix")
  plot_threshold.list(graph, ...)
}

# @export
# @rdname plot_threshold
plot_threshold.list <- function(
  graph, expo, toa, include_censored, t0, attrs,
  undirected, no.contemporary,
  main, xlab, ylab,
  vertex.cex, vertex.col, vertex.label, vertex.lab.pos, vertex.lab.cex,
  vertex.lab.adj,vertex.lab.col,
  vertex.sides, vertex.rot,
  edge.width, edge.col, arrow.length,
  include.grid, bty, vertex.bcol,
  jitter.factor, jitter.amount, xlim, ylim, ...) {
  # Step 0: Getting basic info
  t <- length(graph)
  n <- nrow(graph[[1]])


  # Step 1: Creating the cumulative graph
  # Matrix::sparseMatrix(i={}, j={}, dims=c(n, n))
  cumgraph <- methods::new("dgCMatrix", Dim=c(n,n), p=rep(0L, n+1L))
  for(i in 1:t) {
    cumgraph <- cumgraph + graph[[i]]
  }

  # Creating the pos vector
  y <- threshold(expo, toa, t0, attrs=attrs, include_censored=include_censored)
  y <- jitter(y, factor=jitter.factor[2], amount = jitter.amount[2])

  # Jitter to the xaxis and limits
  jit <- jitter(toa, factor=jitter.factor[1], amount = jitter.amount[1])
  xran <- range(toa, na.rm = TRUE)
  if (!length(xlim)) xlim <- xran + c(-1,1)
  yran <- c(0,1)
  if (!length(ylim)) ylim <- yran + (yran[2] - yran[1])*.1*c(-1,1)

  # Step 2: Checking colors and sizes

  # Computing sizes
  if ((length(vertex.cex) == 1) && inherits(vertex.cex, "character")) {
    if (vertex.cex %in% c("degree", "indegree", "outdegree")) {
      vertex.cex <- dgr(cumgraph, cmode=vertex.cex, undirected = undirected)
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

  # Checking sides
  test <- length(vertex.sides)
  if (!inherits(vertex.sides, c("integer", "numeric"))) {
    stop("-vertex.sides- must be integer.")
  } else if (inherits(vertex.sides, "numeric")) {
    warning("-vertex.sides- will be coerced to integer.")
    vertex.sides <- as.integer(vertex.sides)
  }

  if (test == 1) {
    vertex.sides <- rep(vertex.sides, n)
  } else if (test != n) {
    stop("-vertex.sides- must be of the same length as nnodes(graph).")
  }

  # Checking Rotation
  test <- length(vertex.rot)
  if (!inherits(vertex.rot, "integer") & !inherits(vertex.rot, "numeric")) {
    stop("-vertex.rot- must be numeric.")
  } else if (test == 1) {
    vertex.rot <- rep(vertex.rot, n)
  } else if (test != n) {
    stop("-vertex.rot- must be of the same length as nnodes(graph).")
  }

  # Checking colors
  test <- length(vertex.bcol)
  if (test == 1) vertex.bcol <- rep(vertex.bcol, n)
  else if (test != n) stop("-vertex.bcol- must be either of length 1 or nnodes(graph).")

  # Checking colors
  test <- length(vertex.col)
  if (test == 1) vertex.col <- rep(vertex.col, n)
  else if (test != n) stop("-vertex.col- must be either of length 1 or nnodes(graph).")

  # Plotting
  # oldpar <- par(no.readonly = TRUE)
  plot(NULL, xlim=xlim, ylim=ylim, bty=bty, xlab=xlab, ylab=ylab, main=main,
       xaxs="i", yaxs="i",...)

  # Should there be a grid??
  if (include.grid) grid()

  # Now, for y (it should be different)
  xran <- range(xlim, na.rm = TRUE)
  yran <- range(ylim, na.rm = TRUE)

  # Drawing arrows, first we calculate the coordinates of the edges, for this we
  # use the function edges_coords. This considers aspect ratio of the plot.
  edges <- edges_coords(cumgraph, toa, jit, y, vertex.cex, undirected, no.contemporary,
                        dev=par("pin"), ran=c(xlim[2]-xlim[1], ylim[2]-ylim[1]))
  edges <- as.data.frame(edges)

  # with(edges, segments(x0, y0, x1, y1, lwd = edge.width, col = edge.col))

  ran  <- c(xlim[2]-xlim[1], ylim[2]-ylim[1])
  e_arrows <- apply(edges, 1, function(x) {
      y <- edges_arrow(x["x0"], x["y0"], x["x1"], x["y1"],
                   width=arrow.length,
                   height=arrow.length,
                   beta=pi*(2/3),
                   dev=par("pin"), ran=ran)
      polygon(y[,1], y[,2], col=edge.col, border=edge.col)
    })

  # Drawing the vertices and its labels
  # Computing the coordinates
  pol <- vertices_coords(jit, y, vertex.cex, vertex.sides, vertex.rot, par("pin"), ran)

  # Plotting
  lapply(seq_len(length(pol)),
         function(x) {
           polygon(pol[[x]][,1], pol[[x]][,2],
                   border = vertex.bcol[x],
                   col    = vertex.col[x])
           })

  # Positioning labels can be harsh, so we try with this algorithm
  if (!length(vertex.label)) vertex.label <- 1:n
  text(x=jit, y=y, labels = vertex.label,
       pos = vertex.lab.pos,
       cex = vertex.lab.cex,
       col = vertex.lab.col,
       adj = vertex.lab.adj
       )

  # par(oldpar)

  invisible(data.frame(toa=toa,threshold=y, jit=jit))

}

#' Plot distribution of infect/suscep
#'
#' After calculating infectiousness and susceptibility of each individual on the
#' network, it creates an \code{nlevels} by \code{nlevels} matrix indicating the
#' number of individuals that lie within each cell, and draws a heatmap.
#'
#' @templateVar dynamic TRUE
#' @templateVar toa TRUE
#' @template graph_template
#' @param t0 Integer scalar. See \code{\link{toa_mat}}.
#' @param normalize Logical scalar.  Passed to infection/susceptibility.
#' @param K Integer scalar.  Passed to infection/susceptibility.
#' @param r Numeric scalar.  Passed to infection/susceptibility.
#' @param expdiscount Logical scalar.  Passed to infection/susceptibility.
#' @param bins Integer scalar. Size of the grid (\eqn{n}).
#' @param nlevels Integer scalar. Number of levels to plot (see \code{\link{filled.contour}}).
#' @param h Numeric vector of length 2. Passed to \code{\link[MASS:kde2d]{kde2d}} in the \pkg{MASS} package.
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
#' By default the function will try to apply a kernel smooth function via
#' \code{kde2d}. If not possible (because not enought data points), then
#' the user should try changing the parameter \code{h} or set it equal to zero.
#'
#' \code{toa} is passed to \code{infection/susceptibility}.
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
#' # Generating a random graph -------------------------------------------------
#' set.seed(1234)
#' n <- 100
#' nper <- 20
#' graph <- rgraph_er(n,nper, p=.2, undirected = FALSE)
#' toa <- sample(1:(1+nper-1), n, TRUE)
#'
#' # Visualizing distribution of suscep/infect
#' out <- plot_infectsuscep(graph, toa, K=3, logscale = FALSE)
#' @author George G. Vega Yon
plot_infectsuscep <- function(
  graph, toa, t0=NULL,normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE, bins=20,nlevels=round(bins/2), h=NULL,
  logscale=TRUE, main="Distribution of Infectiousness and\nSusceptibility",
  xlab="Infectiousness of ego", ylab="Susceptibility of ego",
  sub=ifelse(logscale, "(in log-scale)", NA),
  color.palette = grDevices::colorRampPalette(grDevices::blues9),
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
      graph, toa, t0, normalize, K, r, expdiscount, bins, nlevels, h, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued, ...),
    list = plot_infectsuscep.list(
      graph, toa, t0, normalize, K, r, expdiscount, bins, nlevels, h, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued,...),
    diffnet = plot_infectsuscep.list(
      graph$graph, graph$toa, t0, normalize, K, r, expdiscount, bins, nlevels, h, logscale, main,
      xlab, ylab, sub, color.palette, include.grid, exclude.zeros, valued,...),
    stopifnot_graph(graph)
  )
}

# @export
# @rdname plot_infectsuscep
plot_infectsuscep.array <- function(graph, ...) {
  graph <- apply(graph, 3, methods::as, Class="dgCMatrix")
  plot_infectsuscep.list(graph, ...)
}

# @export
# @rdname plot_infectsuscep
plot_infectsuscep.list <- function(graph, toa, t0, normalize,
                              K, r, expdiscount,
                              bins,nlevels,
                              h,
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
  complete <- complete.cases(infect, suscep)

  # Performing classification (linear)
  if (logscale) {
    infectp<-log(infect)
    suscepp<-log(suscep)

    # Only keeping complete cases
    complete <- complete & is.finite(infectp) & is.finite(suscepp)

    if (any(!complete)) warning("When applying logscale some observations are missing.")
  }
  else {
    infectp <- infect
    suscepp <- suscep
  }

  infectp <- infectp[complete,]
  suscepp <- suscepp[complete,]

  if ((!length(infectp) | !length(suscepp)) & logscale)
    stop("Can't apply logscale (undefined values).")

  # If excluding zeros
  include <- rep(TRUE,length(infectp))
  if (exclude.zeros) {
    include[!infectp | !suscepp] <- FALSE
  }


  # Computing infect & suscept
  if (length(h) && h==0) {
    coords <- grid_distribution(infectp[include], suscepp[include], bins)
  } else {
    if (!length(h)) h <- c(
      MASS::bandwidth.nrd(infectp[include & infectp!=0]),
      MASS::bandwidth.nrd(suscepp[include & suscepp!=0])
      )

    # Cant use smoother
    if (any((h==0) | is.na(h)))
      stop('Not enought data to perform smooth. Try choosing another value for -h-,',
           ' or set h=0 (no kernel smooth).')
    coords <- MASS::kde2d(infectp[include], suscepp[include], n = bins, h = h)
  }

  # Nice plot



  n <- sum(coords$z)
  with(coords, filled.contour(
    x,y,
    z/n, bty="n", main=main, xlab=xlab, ylab=ylab, sub=sub, color.palette =color.palette,
    xlim=range(x), ylim=range(y),
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

  # Adding some reference
  legend("topleft", legend=
         sprintf('\n%d out of %d obs.\nincluded', sum(include), length(complete)),
         bty="n")

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
#' # Generating a random diffnet -----------------------------------------------
#' set.seed(821)
#' diffnet <- rdiffnet(100, 5, seed.graph="small-world", seed.nodes="central")
#'
#' plot_adopters(diffnet)
#'
#' # Alternatively, we can use a TOA Matrix
#' toa <- sample(c(NA, 2010L,2015L), 20, TRUE)
#' mat <- toa_mat(toa)
#' plot_adopters(mat$cumadopt)
#' @return A matrix as described in \code{\link{cumulative_adopt_count}}.
#' @export
#' @author George G. Vega Yon
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
    cumadopt <- cumulative_adopt_count(obj)
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

  matplot(times, y=mat, ylim=ylim, add=add, type=type,
          lty=lty, col=col, xlab=xlab, ylab=ylab, main=main, pch=pch,
          bg=bg,...)

  # If not been added
  if (!add) {
    if (include.legend)
      legend("topleft", bty="n", pch=pch,
             legend = c("Cumulative adopters", "Adopters")[test], pt.bg = bg, col=col)

    if (include.grid)
      grid()
  }

  invisible(out)
}

# x <- cumulative_adopt_count(diffnet)
# z <- x["num",] - c(0,x["num",1:(ncol(x)-1)])
# cumsum(z)
# x["num",]


#' \code{diffnet} Arithmetic and Logical Operators
#'
#' Addition, subtraction, network power of diffnet and logical operators such as
#' \code{&} and \code{|} as objects
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
#' @export
#' @name diffnet-arithmetic
#' @family diffnet methods
`^.diffnet` <- function(x,y) {

  if (y < 2) return(x)

  for (i in 1:x$meta$nper) {
    g <- x$graph[[i]]
    for (p in 1:(y-1))
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

  if (inherits(x, "diffnet") && (inherits(y, "numeric") | inherits(y, "integer"))) {
    for (i in 1:x$meta$nper)
      x$graph[[i]]@x <- y/(x$graph[[i]]@x)
    return(x)
  } else if (inherits(y, "diffnet") && (inherits(x, "numeric") | inherits(x, "integer"))) {
    for (i in 1:y$meta$nper)
      y$graph[[i]]@x <- x/(y$graph[[i]]@x)
    return(y)
  } else stop("No method for x:", class(x), " and y:", class(y))

}

#' @rdname diffnet-arithmetic
#' @export
#' @examples
#' # Removing the first 3 vertex of medInnovationsDiffnet ----------------------
#' data(medInnovationsDiffNet)
#'
#' # Using a diffnet object
#' first3Diffnet <- medInnovationsDiffNet[1:3,,]
#' medInnovationsDiffNet - first3Diffnet
#'
#' # Using indexes
#' medInnovationsDiffNet - 1:3
#'
#' # Using ids
#' medInnovationsDiffNet - as.character(1001:1003)
`-.diffnet` <- function(x, y) {
  if (inherits(x, "diffnet") & inherits(y, "diffnet")) {

    # Listing the id numbers that wont be removed
    ids.to.remove <- y$meta$ids
    ids.to.remove <- which(x$meta$ids %in% ids.to.remove)
    x[-ids.to.remove, , drop=FALSE]
  } else if (inherits(x, "diffnet") & any(class(y) %in% c("integer", "numeric"))) {

    # Dropping using ids
    x[-y,, drop=FALSE]
  } else if (inherits(x, "diffnet") & inherits(y, "character")) {
    # Checking labels exists
    test <- which(!(y %in% x$meta$ids))
    if (length(test))
      stop("Some elements in -y- (right-hand side of the expression) are not ",
           "in the set of ids of the diffnet object:\n\t",
           paste0(y[test], collapse=", "),".")

    y <- which(x$meta$ids %in% y)
    x[-y,,drop=FALSE]
  } else
    stop("Subtraction between -",class(x),"- and -", class(y), "- not supported.")
}

#' @export
#' @rdname diffnet-arithmetic
`*.diffnet` <- function(x,y) {
  if (inherits(x, "diffnet") & inherits(y, "diffnet")) {

    # Checking dimensions
    test <- all(dim(x) == dim(y))
    if (!test)
      stop('Both -x- and -y- must have the same dimensions.')

    x$graph <- mapply(`*`, x$graph, y$graph)
    return(x)
  } else if (inherits(x, "diffnet") & is.numeric(y)) {
    x$graph <- mapply(`*`, x$graph, y)
    return(x)

  } else
    stop("Multiplication between -",class(x),"- and -", class(y), "- not supported.")
}

#' @export
#' @rdname diffnet-arithmetic
`&.diffnet` <- function(x,y) {
  x$graph <- mapply(function(a,b) methods::as(a & b, "dgCMatrix"), x$graph, y$graph)
  x
}

#' @export
#' @rdname diffnet-arithmetic
`|.diffnet` <- function(x,y) {
  x$graph <- mapply(function(a,b) methods::as(a | b, "dgCMatrix"), x$graph, y$graph)
  x
}

#' Matrix multiplication
#'
#' Matrix multiplication methods, including \code{\link{diffnet}}
#' objects. This function creates a generic method for \code{\link[base:matmult]{\%*\%}}
#' allowing for multiplying diffnet objects.
#'
#' @param x Numeric or complex matrices or vectors, or \code{diffnet} objects.
#' @param y Numeric or complex matrices or vectors, or \code{diffnet} objects.
#'
#' @details This function can be usefult to generate alternative graphs, for
#' example, users could compute the n-steps graph by doing \code{net \%*\% net}
#' (see examples).
#'
#' @return In the case of \code{diffnet} objects performs matrix multiplication
#' via \code{\link{mapply}} using \code{x$graph} and \code{y$graph} as arguments,
#' returnling a \code{diffnet}. Otherwise returns the default according to
#' \code{\link[base:matmult]{\%*\%}}.
#'
#' @examples
#' # Finding the Simmelian Ties network ----------------------------------------
#'
#' # Random diffnet graph
#' set.seed(773)
#' net <- rdiffnet(100, 4, seed.graph='small-world', rgraph.args=list(k=8))
#' netsim <- net
#'
#' # According to Dekker (2006), Simmelian ties can be computed as follows
#' netsim <- net * t(net) # Keeping mutal
#' netsim <- netsim * (netsim %*% netsim)
#'
#' # Checking out differences (netsim should have less)
#' nlinks(net)
#' nlinks(netsim)
#'
#' mapply(`-`, nlinks(net), nlinks(netsim))
#'
#' @export
#' @rdname diffnetmatmult
#' @family diffnet methods
`%*%` <- function(x, y) UseMethod("%*%")

#' @export
#' @rdname diffnetmatmult
`%*%.default` <- function(x, y) {
  if (inherits(y, "diffnet")) `%*%.diffnet`(x,y)
  else base::`%*%`(x=x,y=y)
}

#' @export
#' @rdname diffnetmatmult
`%*%.diffnet` <- function(x, y) {

  mat2dgCList <- function(w,z) {
    w <- lapply(seq_len(nslices(z)), function(u) methods::as(w, "dgCMatrix"))
    names(w) <- dimnames(z)[[3]]
    w
  }

  if (inherits(x, "diffnet") && inherits(y, "diffnet")) {
    x$graph <- mapply(base::`%*%`, x$graph, y$graph)
  } else if (inherits(x, "diffnet") && !inherits(y, "diffnet")) {
    if (identical(rep(dim(x)[1],2), dim(y)))
      x$graph <- mapply(base::`%*%`, x$graph, mat2dgCList(y, x))
    else stop("-y- must have the same dimension as -x-")
  } else if (inherits(y, "diffnet") && !inherits(x, "diffnet")) {
    if (identical(rep(dim(y)[1],2), dim(x))) {
      y$graph <- mapply(base::`%*%`, mat2dgCList(x, y), y$graph)
      return(y)
    }
    else stop("-y- must have the same dimension as -x-")
  }

  x
}



#' Coerce a diffnet graph into an array
#'
#' @param x A diffnet object.
#' @param ... Ignored.
#' @details
#' The function takes the list of sparse matrices stored in \code{x} and creates
#' an array with them. Attributes and other elements from the diffnet object are
#' dropped.
#'
#' \code{dimnames} are obtained from the metadata of the diffnet object.
#'
#' @return A three-dimensional array of \eqn{T} matrices of size \eqn{n\times n}{n * n}.
#' @seealso \code{\link{diffnet}}.
#' @family diffnet methods
#' @examples
#' # Creating a random diffnet object
#' set.seed(84117)
#' mydiffnet <- rdiffnet(30, 5)
#'
#' # Coercing it into an array
#' as.array(mydiffnet)
#' @export
as.array.diffnet <- function(x, ...) {
  # Coercing into matrices
  z <- lapply(x$graph, function(y) {
    as.matrix(y)
  })

  # Creating the array
  out <- with(x$meta, array(dim=c(n, n, nper)))
  for (i in 1:length(z))
    out[,,i] <- z[[i]]

  # Naming dimensions
  dimnames(out) <- with(x$meta, list(ids, ids, pers))
  out
}


#' Count the number of vertices/edges/slices in a graph
#'
#' @template graph_template
#' @return For \code{nvertices} and \code{nslices}, an integer scalar equal to the number
#' of vertices and slices in the graph. Otherwise, from \code{nedges}, either a list
#' of size \eqn{t} with the counts of edges (non-zero elements in the adjacency matrices) at
#' each time period, or, when \code{graph} is static, a single scalar with
#' such number.
#' @details
#' \code{nnodes} and \code{nlinks} are just aliases for \code{nvertices} and
#' \code{nedges} respectively.
#' @export
#' @examples
#' # Creating a dynamic graph (we will use this for all the classes) -----------
#' set.seed(13133)
#' diffnet <- rdiffnet(100, 4)
#'
#' # Lets use the first time period as a static graph
#' graph_mat <- diffnet$graph[[1]]
#' graph_dgCMatrix <- methods::as(graph_mat, "dgCMatrix")
#'
#' # Now lets generate the other dynamic graphs
#' graph_list  <- diffnet$graph
#' graph_array <- as.array(diffnet) # using the as.array method for diffnet objects
#'
#' # Now we can compare vertices counts
#' nvertices(diffnet)
#' nvertices(graph_list)
#' nvertices(graph_array)
#'
#' nvertices(graph_mat)
#' nvertices(graph_dgCMatrix)
#'
#' # ... and edges count
#' nedges(diffnet)
#' nedges(graph_list)
#' nedges(graph_array)
#'
#' nedges(graph_mat)
#' nedges(graph_dgCMatrix)
nvertices <- function(graph) {
  switch(class(graph),
         array     = nrow(graph),
         matrix    = nrow(graph),
         dgCMatrix = nrow(graph),
         list      = nrow(graph[[1]]),
         diffnet   = graph$meta$n,
         stopifnot_graph(graph)
         )
}

#' @rdname nvertices
#' @export
nnodes <- nvertices

#' @export
#' @rdname nvertices
nedges <- function(graph) {
  switch (class(graph),
    array     = {
      # Computing and coercing into a list
      x <- as.list(apply(graph, 3, function(x) sum(x!=0)))

      # Naming
      tnames <- names(x)
      if (!length(tnames)) names(x) <- 1:length(x)
      x
      },
    matrix    = sum(graph != 0),
    dgCMatrix = length(graph@i),
    list      = {
      # Computing
      x <- lapply(graph, function(x) length(x@i))

      # Naming
      tnames <- names(x)
      if (!length(tnames)) names(x) <- 1:length(x)
      x
      },
    diffnet   = lapply(graph$graph, function(x) sum(x@x != 0)),
    stopifnot_graph(graph)
  )
}

#' @export
#' @rdname nvertices
nlinks <- nedges

#' @export
#' @rdname nvertices
nslices <- function(graph) {
  switch (class(graph),
    array     = dim(graph)[3],
    matrix    = 1L,
    dgCMatrix = 1L,
    diffnet   = graph$meta$nper,
    list      = length(graph),
    stopifnot_graph(graph)
  )
}

#' @export
#' @rdname diffnet-class
nodes <- function(graph) {
  cls <- class(graph)
  if ("diffnet" %in% cls)
    return(graph$meta$ids)
  else if ("list" %in% cls) {
    ans <- rownames(graph[[1]])
    if (!length(ans)) stop("There are not names to fetch")
    else return(ans)
  } else if (any(c("matrix", "dgCMatrix", "array") %in% cls)) {
    ans <- rownames(graph)
    if (!length(ans)) stop("There are not names to fetch")
    else return(ans)
  }
  else stopifnot_graph(graph)

}

#' @export
#' @rdname diffnet-class
#' @param FUN a function to be passed to lapply
diffnetLapply <- function(graph, FUN, ...) {
  lapply(seq_len(nslices(graph)), function(x, graph, ...) {
    FUN(x,
        graph               = graph$graph[[x]],
        toa                 = graph$toa,
        vertex.static.attrs = graph$vertex.static.attrs,
        vertex.dyn.attrs    = graph$vertex.dyn.attrs[[x]],
        adopt               = graph$adopt[,x,drop=FALSE],
        cumadopt            = graph$cumadopt[,x,drop=FALSE],
        meta                = graph$meta)
    }, graph=graph,...)
}
# debug(diffnetLapply)
# diffnetLapply(medInnovationsDiffNet, function(x, graph, cumadopt, ...) {
#   sum(cumadopt)
# })

#' @export
#' @rdname diffnet-class
str.diffnet <- function(object, ...) {
  utils::str(unclass(object))
}

#' @export
#' @rdname diffnet-class
dimnames.diffnet <- function(x) {
  with(x, list(
    meta$ids,
    c(colnames(vertex.static.attrs), names(vertex.dyn.attrs[[1]])),
    meta$pers)
  )
}

#' @export
#' @rdname diffnet-class
#' @method t diffnet
t.diffnet <- function(x) {
  x$graph <- lapply(x$graph, getMethod("t", "dgCMatrix"))
  x
}

#' @rdname diffnet-class
#' @export
dim.diffnet <- function(x) {
  k <- length(with(x, c(colnames(vertex.static.attrs), names(vertex.dyn.attrs[[1]]))))
  as.integer(with(x$meta, c(n, k, nper)))
}

