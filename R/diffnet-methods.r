#' S3 plotting method for diffnet objects.
#'
#' @param x An object of class \code{\link[=diffnet-class]{diffnet}}
#' @param t Integer scalar indicating the time slice to plot.
#' @param vertex.color Character scalar/vector. Color of the vertices.
#' @template plotting_template
#' @param main Character. A title template to be passed to sprintf.
#' @param ... Further arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param y Ignored.
#' @export
#'
#' @family diffnet methods
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
  x,y=NULL, t=1,
  vertex.color  = c(adopt="steelblue", noadopt="white"),
  vertex.size = "degree",
  main        = "Diffusion network in time %d",
  minmax.relative.size = getOption("diffnet.minmax.relative.size", c(0.01, 0.04)),
  ...) {

  # Listing arguments
  igraph.args <- list(...)

  # Checking that the time period is actually within
  if (!(t %in% 1:x$meta$nper))
    stop("-t- must be an integer within 1 and ",x$meta$nper,".")

  # Extracting the graph to be plotted
  graph <- diffnet_to_igraph(x)[[t]]

  # Setting the colors
  cols <- with(x, ifelse(cumadopt[,t], vertex.color[1], vertex.color[2]))

  set_igraph_plotting_defaults("igraph.args")

  if (!length(igraph.args$layout))
    igraph.args$layout <- igraph::layout_nicely(graph)

  igraph.args$vertex.color  <- cols

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-1,1),
    ylim = c(-1,1)
    )

  igraph.args$vertex.size <-
    rescale_vertex_igraph(
      compute_vertex_size(x$graph[[t]], vertex.size),
      minmax.relative.size = minmax.relative.size
      )

  do.call(igraph::plot.igraph, c(
    list(
      x           = graph
    ), igraph.args))

  if (length(main))
    graphics::title(main = sprintf(main, x$meta$pers[t]))

  invisible(igraph.args$layout)
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
#' @template graph_template
#' @param cumadopt \eqn{n\times T}{n*T} matrix.
#' @param slices Integer vector. Indicates what slices to plot. By default all are plotted.
#' @param vertex.color A character vector of size 3 with colors names.
#' @param vertex.shape A character vector of size 3 with shape names.
#' @template plotting_template
#' @param mfrow.par Vector of size 2 with number of rows and columns to be passed to \code{\link{par}.}
#' @param main Character scalar. A title template to be passed to \code{\link{sprintf}.}
#' @param ... Further arguments to be passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param legend.args List of arguments to be passed to \code{\link{legend}}.
#' @param background Either a function to be called before plotting each slice, a color
#' to specify the backgroupd color, or \code{NULL} (in which case nothing is done).
#'
#' @details Plotting is done via the function \code{\link[igraph:plot.igraph]{plot.igraph}}.
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
#' The argument \code{vertex.color} contains the colors of non-adopters, new-adopters,
#' and adopters respectively. The new adopters (default color \code{"tomato"}) have a different
#' color that the adopters when the graph is at their time of adoption, hence,
#' when the graph been plotted is in \eqn{t=2} and \eqn{toa=2} the vertex will
#' be plotted in red.
#'
#' \code{legend.args} has the following default parameter:
#' \tabular{ll}{
#'   \code{x} \tab \code{"bottom"} \cr
#'   \code{legend} \tab \code{c("Non adopters", "New adopters","Adopters")} \cr
#'   \code{pch} \tab \code{sapply(vertex.shape, switch, circle = 21, square = 22, 21)} \cr
#'   \code{bty} \tab \code{"n"} \cr
#'   \code{horiz} \tab \code{TRUE} \cr
#' }
#'
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
plot_diffnet <- function(...) UseMethod("plot_diffnet")

#' @export
#' @rdname plot_diffnet
plot_diffnet.diffnet <- function(
  graph, ...
) {

  args <- list(...)

  do.call(
    plot_diffnet.default,
    c(
      list(graph = as_dgCMatrix(graph), cumadopt = graph$cumadopt),
      args
    )
  )


}

#' @rdname plot_diffnet
#' @export
plot_diffnet.default <- function(
  graph, cumadopt,
  slices       = NULL,
  vertex.color = c("white", "tomato", "steelblue"),
  vertex.shape = c("square", "circle", "circle"),
  vertex.size  = "degree",
  mfrow.par    = NULL,
  main         = c("Network in period %s", "Diffusion Network"),
  legend.args  = list(),
  minmax.relative.size = getOption("diffnet.minmax.relative.size", c(0.01, 0.04)),
  background   = NULL,
  ...) {

  set_plotting_defaults("background")

  # Setting parameters
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Setting legend parameters, if specified
  if (length(legend.args) | (!length(legend.args) & is.list(legend.args))) {
    if (!length(legend.args$x)) legend.args$x <- "bottom"
    if (!length(legend.args$legend))
      legend.args$legend <-c("Non adopters", "New adopters","Adopters")
    if (!length(legend.args$pch)) {
      legend.args$pch <- sapply(vertex.shape, switch, circle = 21, square = 22, 21)
    }
    if (!length(legend.args$bty)) legend.args$bty <- "n"
    if (!length(legend.args$horiz)) legend.args$horiz <-TRUE
  }


  igraph.args <- list(...)

  # Coercing into a dgCMatrix list
  graph <- as_dgCMatrix(graph)

  if (!is.list(graph))
    stopifnot_graph(graph)

  # Making sure it has names
  add_dimnames.list(graph)
  colnames(cumadopt) <- names(graph)

  # Checking parameters
  t <- nslices(graph)
  n <- nrow(graph[[1]])

  # Checking slices
  if (!length(slices)) {
    slices <- names(graph)[unique(floor(seq(1, t, length.out = min(t, 4))))]
  } else if (is.numeric(slices)) {
    slices <- names(graph)[slices]
  }

  t <- length(slices)

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


  # Computing legend and main width/height
  legend_height_i <- 0
  if (length(legend.args) && length(legend.args$legend)) {
    legend_height_i <- max(sapply(
      legend.args$legend,
      graphics::strheight,
      units="inches",
      cex = if (length(legend.args$cex)) legend.args$cex else NULL
      ))*2.5
  }

  main_height_i <- graphics::strheight(
    main[2],
    units = "inches",
    cex = if ("cex.main" %in% igraph.args) igraph.args$main.cex else NULL
  )*1.5

  graphics::par(
    mfrow = mfrow.par, mar = rep(.25,4),
    omi = c(legend_height_i, 0, main_height_i, 0),
    xpd = NA, xaxs = "i", yaxs="i"
    )

  # Setting igraph defaults
  set_igraph_plotting_defaults("igraph.args")

  # 3. Plotting ----------------------------------------------------------------
  times <- as.integer(names(graph))

  # Set types:
  # - 1: Non adopter
  # - 2: Adopter in s
  # - 3: Adopter prior to s
  set_type <- function() {

    i <- match(s, colnames(cumadopt))
    j <- match(s, slices)

    # If we are looking at the first of both
    if (i==1 & j ==1)
      return(ifelse(!cumadopt[,s], 1L, 2L))

    # Otherwise, we look at something more complicated
    type <- ifelse(!cumadopt[,s] , 1L, NA)

    if (j > 1) {
      type <- ifelse(!is.na(type), type,
                     ifelse(cumadopt[,slices[j-1]], 3L, 2L))
    } else if (i > 1) {
      type <- ifelse(!is.na(type), type,
                     ifelse(cumadopt[, i-1], 3L, 2L))
    }

    type

  }

  for (s in slices) {
    # Colors, new adopters are painted differently

    # Setting color and shape depending on the type of vertex these are.
    type   <- set_type()
    cols   <- vertex.color[type]
    shapes <- vertex.shape[type]

    # Creating igraph object
    ig  <- igraph::graph_from_adjacency_matrix(graph[[s]], weighted = TRUE)

    # Computing layout
    if (!length(igraph.args$layout)) {
      igraph.args$layout <- igraph::layout_nicely(ig)
    } else if (length(igraph.args$layout) && is.function(igraph.args$layout)) {
      igraph.args$layout <- igraph.args$layout(ig)
    }

    # Computing subtitle height
    graphics::plot.new()
    graphics::plot.window(xlim=c(-1.15,1.15), ylim=c(-1.15,1.15))

    # Should we paint or do something else?
    if (is.function(background)) background()
    else if (length(background))
      graphics::rect(-1.15,-1.15,1.15,1.15, col=background, border=background)

    # Plotting
    do.call(
      igraph::plot.igraph,
      c(
      list(
        ig,
        vertex.color = cols,
        vertex.size  = rescale_vertex_igraph(
          compute_vertex_size(graph, vertex.size, match(s, names(graph))),
          minmax.relative.size = minmax.relative.size
          ),
        vertex.shape = shapes
        ),
      igraph.args)
      )

    # Adding a legend (title)
    if (length(main))
      subtitle(x = sprintf(main[1], names(graph[s])))

  }

  # Legend
  graphics::par(
    mfrow = c(1,1), mai = rep(0,4), new = TRUE, xpd=NA,
    omi = c(0, 0, main_height_i, 0)
  )

  # graphics::par(mfrow=c(1,1), new=TRUE, mar=rep(0,4), oma = rep(0,4), xpd=NA)
  graphics::plot.new()
  graphics::plot.window(c(0,1), c(0,1))
  if (length(main) > 1)
    title(main = main[2], outer=TRUE)

  if (length(legend.args))
    do.call(graphics::legend, c(legend.args, list(pt.bg=vertex.color)))

  invisible(igraph.args$layout)

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
#' @param vertex.size Numeric vector of size \eqn{n}. Relative size of the vertices.
#' @param vertex.color Either a vector of size \eqn{n} or a scalar indicating colors of the vertices.
#' @param vertex.label Character vector of size \eqn{n}. Labels of the vertices.
#' @param vertex.label.pos Integer value to be passed to \code{\link{text}} via \code{pos}.
#' @param vertex.label.cex Either a numeric scalar or vector of size \eqn{n}. Passed to \code{text}.
#' @param vertex.label.adj Passed to \code{\link{text}}.
#' @param vertex.label.color Passed to \code{\link{text}}.
#' @param jitter.amount Numeric vector of size 2 (for x and y) passed to \code{\link{jitter}}.
#' @param jitter.factor Numeric vector of size 2 (for x and y) passed to \code{\link{jitter}}.
#' @param vertex.frame.color Either a vector of size \eqn{n} or a scalar indicating colors of vertices' borders.
#' @param vertex.sides Either a vector of size \eqn{n} or a scalar indicating the
#' number of sides of each vertex (see details).
#' @param vertex.rot Either a vector of size \eqn{n} or a scalar indicating the
#' rotation in radians of each vertex (see details).
#' @param edge.width Numeric. Width of the edges.
#' @param edge.color Character. Color of the edges.
#' @param arrow.width Numeric value to be passed to \code{\link{arrows}}.
#' @param arrow.length Numeric value to be passed to \code{\link{arrows}}.
#' @param arrow.color Color.
#' @param include.grid Logical. When TRUE, the grid of the graph is drawn.
#' @param bty See \code{\link{par}}.
#' @param xlim Passed to \code{\link{plot}}.
#' @param ylim Passed to \code{\link{plot}}.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @param edge.curved Logical scalar. When curved, generates curved edges.
#' @param background TBD
#' @family visualizations
#' @seealso Use \code{\link{threshold}} to retrieve the corresponding threshold
#' obtained returned by \code{\link{exposure}}.
#' @keywords hplot
#'
#' @details When \code{vertex.label=NULL} the function uses vertices ids as labels.
#' By default \code{vertex.label=""} plots no labels.
#'
#' Vertices are drawn using an internal function for generating polygons.
#' Polygons are inscribed in a circle of radius \code{vertex.size}, and can be
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
#' plot_threshold(graph, expos, toa, vertex.size = "indegree")
#'
#' @export
#' @author George G. Vega Yon
plot_threshold <- function(graph, expo, ...) UseMethod("plot_threshold")

#' @export
#' @rdname plot_threshold
plot_threshold.diffnet <- function(graph, expo, ...) {
  # If graph is diffnet, then we should do something different (because the
  # first toa may not be the firts one as toa may be stacked to the right.
  # see ?as_diffnet)
  # graph$toa <- graph$toa - min(graph$meta$pers) + 1L

  if (missing(expo))
    expo <- exposure(graph)

  args <- list(...)

  if (!length(args$undirected)) args$undirected <- graph$meta$undirected
  if (!length(args$t0)) args$t0 <- graph$meta$pers[1]
  if (length(args$toa)) {
    warning("While -graph- has its own toa variable, the user is providing one.")
  } else {
    args$toa <- graph$toa
  }

  do.call(plot_threshold.default, c(list(graph = graph$graph, expo=expo), args))
}

#' @export
#' @rdname plot_threshold
plot_threshold.array <- function(graph, expo, ...) {
  plot_threshold.default(as_dgCMatrix(graph), expo = expo, ...)
}

#' @export
#' @rdname plot_threshold
plot_threshold.default <- function(
  graph,
  expo,
  toa,
  include_censored = FALSE,
  t0               = min(toa, na.rm = TRUE),
  attrs            = NULL,
  undirected       = getOption("diffnet.undirected"),
  no.contemporary  = TRUE,
  main             = "Time of Adoption by\nNetwork Threshold",
  xlab             = "Time",
  ylab             = "Threshold",
  vertex.size      = "degree",
  vertex.color     = NULL,
  vertex.label     = "",
  vertex.label.pos = NULL,
  vertex.label.cex = 1,
  vertex.label.adj = c(.5,.5),
  vertex.label.color = NULL,
  vertex.sides     = 40L,
  vertex.rot       = 0,
  edge.width       = 2,
  edge.color       = NULL,
  arrow.width      = NULL,
  arrow.length     = NULL,
  arrow.color      = NULL,
  include.grid     = FALSE,
  vertex.frame.color = NULL,
  bty              = "n",
  jitter.factor    = c(1,1),
  jitter.amount    = c(.25,.025),
  xlim             = NULL,
  ylim             = NULL,
  edge.curved      = NULL,
  background       = NULL,
  ...
  ) {

  # Setting default parameters
  set_plotting_defaults(c("edge.color", "vertex.frame.color", "vertex.label.color", "edge.curved", "vertex.color", "background", "arrow.color"))

  # # Checking out defaults
  # if (!length(edge.color)) edge.color <- igraph_plotting_defaults$edge.color
  # if (!length(edge.color)) edge.color <- igraph_plotting_defaults$vertex.frame.color

  # Checking if exposure was provided
  if (missing(expo))
    stop("expo should be provided")

  # Checking the type of graph
  graph <- as_dgCMatrix(graph)

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
  y0 <- threshold(expo, toa, t0, attrs=attrs, include_censored=include_censored)
  y <- jitter(y0, factor=jitter.factor[2], amount = jitter.amount[2])

  # Jitter to the xaxis and limits
  jit <- jitter(toa, factor=jitter.factor[1], amount = jitter.amount[1])
  xran <- range(toa, na.rm = TRUE)
  if (!length(xlim)) xlim <- xran + c(-1,1)
  yran <- c(0,1)
  if (!length(ylim)) ylim <- yran + (yran[2] - yran[1])*.1*c(-1,1)

  # Step 2: Checking colors and sizes

  # Computing sizes
  vertex.size <- compute_vertex_size(graph, vertex.size)

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

  # Plotting
  # oldpar <- par(no.readonly = TRUE)
  graphics::plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main,
       xaxs="i", yaxs="i",...)

  # Should we paint or do something else?
  if (is.function(background)) background()
  else if (length(background))
    graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2], col=background, border=background)

  # Checking
  if (!length(arrow.width))
    arrow.width <- with(graphics::par(), (usr[2] - usr[1])/75)

  if (!length(arrow.length))
    arrow.length <- with(graphics::par(), (usr[2] - usr[1])/75)

  # Should there be a grid??
  if (include.grid)
    grid()

  # Now, for y (it should be different)
  xran <- range(xlim, na.rm = TRUE)
  yran <- range(ylim, na.rm = TRUE)

  # Drawing arrows, first we calculate the coordinates of the edges, for this we
  # use the function edges_coords. This considers aspect ratio of the plot.
  vertex.size <- igraph_vertex_rescale(vertex.size, adjust=1)
  edges <- edges_coords(cumgraph, toa, jit, y, vertex.size, undirected, no.contemporary,
                        dev=par("pin"), ran=c(xlim[2]-xlim[1], ylim[2]-ylim[1]))
  edges <- as.data.frame(edges)

  ran  <- c(xlim[2]-xlim[1], ylim[2]-ylim[1])

  # Plotting the edges
  mapply(function(x0, y0, x1, y1, col, edge.curved, arrow.color) {
    y <- edges_arrow(x0, y0, x1, y1, width=arrow.width, height=arrow.length,
                     beta=pi*(2/3), dev=par("pin"), ran=ran, curved = edge.curved)

    # Drawing arrow
    if (edge.curved) {

      # Edge
      graphics::xspline(
        y$edge[,1],y$edge[,2],
        shape = c(0, 1, 0),
        open=TRUE, border = col, lwd=edge.width)

      # Arrow
      graphics::polygon(y$arrow[,1], y$arrow[,2], col = arrow.color, border = arrow.color)
    } else {
      # Edge
      graphics::polygon(y$edge[,1],y$edge[,2], col = col, border = col, lwd=edge.width)

      # Arrow
      graphics::polygon(y$arrow[,1], y$arrow[,2], col = arrow.color, border = arrow.color)
    }


  }, x0 = edges[,"x0"], y0 = edges[,"y0"], x1 = edges[,"x1"], y1 = edges[,"y1"],
  col = edge.color, edge.curved = edge.curved, arrow.color=arrow.color)

  # Drawing the vertices and its labels
  # Computing the coordinates
  pol <- vertices_coords(jit, y, vertex.size, vertex.sides, vertex.rot, par("pin"), ran)

  # Plotting
  mapply(function(coords,border,col)
    graphics::polygon(coords[,1], coords[,2], border = border, col=col),
    coords = pol, border = vertex.frame.color, col=vertex.color)

  # Positioning labels can be harsh, so we try with this algorithm
  if (!length(vertex.label)) vertex.label <- 1:n
  graphics::text(x=jit, y=y, labels = vertex.label,
       pos = vertex.label.pos,
       cex = vertex.label.cex,
       col = vertex.label.color,
       adj = vertex.label.adj
       )

  # par(oldpar)

  invisible(data.frame(toa=toa,threshold=y0, jit=jit))

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
  graph,
  toa,
  t0            = NULL,
  normalize     = TRUE,
  K             = 1L,
  r             = 0.5,
  expdiscount   = FALSE,
  bins          = 20,
  nlevels       = round(bins/2),
  h             = NULL,
  logscale      = TRUE,
  main          = "Distribution of Infectiousness and\nSusceptibility",
  xlab          = "Infectiousness of ego",
  ylab          = "Susceptibility of ego",
  sub           = ifelse(logscale, "(in log-scale)", NA),
  color.palette = function(n) viridisLite::viridis(n),
  include.grid  = TRUE,
  exclude.zeros = FALSE,
  valued        = getOption("diffnet.valued",FALSE),
  ...
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
plot_adopters <- function(
  obj,
  freq           = FALSE,
  what           = c("adopt","cumadopt"),
  add            = FALSE,
  include.legend = TRUE,
  include.grid   = TRUE,
  pch            = c(21,24),
  type           = c("b", "b"),
  ylim           = if (!freq) c(0,1) else NULL,
  lty            = c(1,1),
  col            = c("black","black"),
  bg             = c("tomato","gray"),
  xlab           = "Time",
  ylab           = ifelse(freq, "Frequency", "Proportion"),
  main           = "Adopters and Cumulative Adopters",
  ...
  ) {

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
         igraph    = igraph::vcount(graph),
         network   = network::network.size(graph),
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
    igraph    = igraph::ecount(graph),
    network   = network::network.edgecount(graph),
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

