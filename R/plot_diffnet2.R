#' Takes a numeric vector and maps it into a finite length sequence
#' @param x A numeric or integer vector.
#' @param nlevels Integer scalar. Length of the sequence to be map onto.
#' @param as_factor Logical scalar. When \code{TRUE} the resulting vector is factor.
#' @return A vector of length \code{length(x)} with values mapped to a sequence
#' with \code{nlevels} unique valuess
#' @export
#' @examples
#'
#' x <- rnorm(100)
#' w <- data.frame(as.integer(round_to_seq(x, as_factor = TRUE)),x)
#' plot(w,x)
round_to_seq <- function(x, nlevels=20, as_factor=TRUE) {
  y <- range(x, na.rm = TRUE, finite=TRUE)
  y <- seq(y[1], y[2], length.out = nlevels)

  y <- sapply(x, function(z) {
    if (is.na(z)) return(NA)
    y[which.min(abs(y-z))]
    })

  if (as_factor) as.factor(y)
  else y
}


#' Another way of visualizing diffusion
#' @param graph Either a square matrix or a diffnet object.
#' @param toa Integer vector of length \eqn{n} with the times of adoption.
#' @param pers Integer vector of length \eqn{T} indicating the time periods of the data.
#' @param color.ramp A function as returned by \code{\link[grDevices::colorRamp]{colorRamp}}.
#' @param layout Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param key.width Numeric scalar. Sets the proportion of the plot (x-axis) that the key uses.
#' @param key.title Character scalar. Title of the key (vertex colors).
#' @param main Character scalar. Title of the graph.
#' @param vertex.size Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param vertex.shape Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param vertex.label Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param vertex.frame.color Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param edge.arrow.size Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param edge.curved Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param add.map Character scalar. When \code{"first"} plots a \code{\link{diffusionMap}} before the
#'  graph itself. If \code{"last"} then it adds it at the end. When \code{NULL} adds nothing.
#' @param diffmap.args List. If \code{add.map=TRUE}, arguments passed to \code{diffusionMap}.
#' @param diffmap.alpha Numeric scalar between [0,1]. Alpha level for the map.
#' @param include.white Character scalar. Includes white in the color palette used in the map.
#'  When \code{include.white=NULL} then it won't include it.
#' @param ... Further arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @details If \code{key.width<=0} then no key is created.
#'
#' @return A list with the following elements
#' \item{layout}{}
#' \item{vertex.color}{}
#' \item{vertex.label}{}
#' \item{vertex.shape}{}
#' \item{vertex.size}{}
#' \item{diffmap}{}
#' @export
#' @family visualizations
#' @author George G. Vega Yon
plot_diffnet2 <- function(...) UseMethod("plot_diffnet2")

#' @rdname plot_diffnet2
#' @export
#' @include diffnet-methods.R data.R
plot_diffnet2.diffnet <- function(
  graph, toa=NULL,
  slice=nslices(graph),
  color.ramp = grDevices::colorRamp(c("skyblue","yellow", "red")),
  layout = NULL,
  key.width = .10,
  key.title = "Time of Adoption",
  main = "Diffusion dynamics",
  vertex.size = NULL,
  vertex.shape = NULL,
  vertex.label = "",
  vertex.frame.color="gray",
  edge.arrow.size=.5,
  edge.curved=FALSE,
  add.map = FALSE,
  diffmap.args=list(kde2d.args=list(n=100)),
  diffmap.alpha=.5,
  include.white="first",
  ...
) {
  plot_diffnet2.default(
    graph=graph$graph[[slice]], toa=if (length(toa)) toa else graph$toa, pers=graph$meta$pers,
    color.ramp=color.ramp, layout=layout, key.width=key.width, key.title=key.title,
    main=main,
    vertex.size=vertex.size, vertex.shape=vertex.shape, vertex.label=vertex.label,
    vertex.frame.color=vertex.frame.color, edge.arrow.size=edge.arrow.size,
    edge.curved=edge.curved, add.map=add.map, diffmap.args=diffmap.args,diffmap.alpha,include.white,...)
}

#' @rdname plot_diffnet2
#' @export
plot_diffnet2.default <- function(
  graph,
  toa,
  pers = min(toa, na.rm = TRUE):max(toa, na.rm = TRUE),
  color.ramp = grDevices::colorRamp(c("skyblue","yellow", "red")),
  layout = NULL,
  key.width = .10,
  key.title = "Time of Adoption",
  main = "Diffusion dynamics",
  vertex.size = NULL,
  vertex.shape = NULL,
  vertex.label = "",
  vertex.frame.color="gray",
  edge.arrow.size=.5,
  edge.curved=FALSE,
  add.map=FALSE,
  diffmap.args=list(kde2d.args=list(n=100)),
  diffmap.alpha=.5,
  include.white = "first",
  ...) {

  # Some constants
  nper <- length(pers)

  # Taggin types ---------------------------------------------------------------

  # 1st adopters
  type_1st <- toa == pers[nper]

  # Non Adopters
  type_non <- is.na(toa)

  # Adopters
  type_adopt <- which(!type_1st & !type_non)
  type_1st   <- which(type_1st)
  type_non   <- which(type_non)

  # Colors
  t01 <- pers
  t01 <- c(t01[1], t01[nper])
  col <- color.ramp( (toa - t01[1])/(t01[2] - t01[1]) )
  col[type_non,] <- 255
  col <- rgb(col[,1], col[,2], col[,3], maxColorValue = 255)

  # Shapes
  if (!length(vertex.shape)) {
    vertex.shape <- rep("circle", nnodes(graph))
    vertex.shape[type_non] <- "square"
  }

  # Computing layout -----------------------------------------------------------
  g <- igraph::graph_from_adjacency_matrix(graph, mode="undirected")

  l <- if (!length(layout)) igraph::layout_nicely(g)
  else if (inherits(layout, "function")) layout(g)
  else layout

  xran <- range(l[,1])
  yran <- range(l[,2])

  # Adjusting
  xran[2] <- (xran[2] - xran[1]*(key.width + .1))/(1 - key.width - .1)

  plot.new()
  plot.window(xlim=xran, ylim=yran, xaxs="i", yaxs="i")
  title(main=main)

  # If adding map! -------------------------------------------------------------
  if (length(add.map)) {
    dm <- do.call(diffusionMap.default, c(diffmap.args, list(graph=graph, x=toa,
                                                             layout = l)))
    # Levels
    dmlvls <- pretty(range(dm$map$z), diffmap.args$kde2d.args$n)

    # Colors, in this case we need to extrapolate nper and add white.
    dmcol <- grDevices::rgb(color.ramp(seq(0,1, length.out = nper*2)), maxColorValue = 255)

    # Do we need to include white in the map?
    if (length(include.white))
      if (include.white=="first") dmcol <- c("white", dmcol)
      else if (include.white=="last") dmcol <- c(dmcol, "white")
      else stop('-include.white- should be either NULL, "first" or "last".')

    # Palette
    dmcol <- grDevices::adjustcolor(grDevices::colorRampPalette(dmcol)(length(dmlvls)),
                                    alpha.f=diffmap.alpha)

    # Plot
    if (add.map=="first")
      .filled.contour(dm$map$x, dm$map$y, dm$map$z, levels = dmlvls, col=dmcol)
  } else dm <- NULL

  # Plotting graph -------------------------------------------------------------
  igraph::plot.igraph(
    g, layout=l,
    vertex.color=col,
    vertex.label=vertex.label,
    vertex.shape=vertex.shape,
    vertex.size=vertex.size,
    vertex.frame.color=vertex.frame.color,
    edge.arrow.size=edge.arrow.size,
    edge.curved=edge.curved,
    add=TRUE, rescale=FALSE,
    xlim=xran, ylim=yran,
    ...
  )

  if (length(add.map) && (add.map=="last"))
      .filled.contour(dm$map$x, dm$map$y, dm$map$z, levels = dmlvls, col=dmcol)

  # # Plotting boxes -------------------------------------------------------------
  if (key.width > 0)
    drawColorKey(toa, key.pos = c(1-key.width, 0.975, 0.05, 0.95),
                 nlevels = 100, main = key.title, border="transparent")

  invisible(list(layout=l,vertex.color=col,
                 vertex.label=vertex.label,
                 vertex.shape=vertex.shape,
                 vertex.size=vertex.size,
                 diffmap=dm))
}


#' Creates a heatmap based on a graph layout and times of adoption
#'
#' Basically creates a smooth-scatter plot in which each observation is weighted
#' by \code{x}.
#'
#' @param graph A square matrix of size \eqn{n\times n}{n * n}.
#' @param x An vector of length \eqn{n}. Usually a \code{toa} vector.
#' @param layout Either a \eqn{n\times 2}{n *2} matrix of coordinates or a layout
#'  function applied to \code{graph} (must return coordinates).
#' @param jitter.args A list including arguments to be passed to \code{\link{jitter}}.
#' @param kde2d.args A list including arguments to be passed to \code{\link[MASS:kde2d]{kde2d}}.
#' @details
#' The image is created using the function \code{kde2d} from
#' the \pkg{MASS} package. The complete algorithm follows:
#' \enumerate{
#'  \item \code{x} is coerced into integer and the range is adjusted to start from 1.
#'    \code{NA} are replaced by zero.
#'  \item In no \code{layout} is passed, layout is computed using
#'    \code{\link[igraph:layout_nicely]{layout_nicely}} from \pkg{igraph}
#'  \item Each vertex's coordinates is repeated \code{x} times.
#'  \item The jitter function is applied to the repeated coordinates.
#'  \item 2D kernel is computed using \code{kde2d} over the coordinates.
#' }
#'
#' The resulting matrix can be passed to \code{\link{image}} or similar.
#' @return A list of class \code{diffnet_diffmap}
#' \item{coords}{A matrix of size \eqn{n\times 2}{n*2} of vertices coordinates.}
#' \item{map}{Output from \code{kde2d}. This is a list with 3 elements, vectors
#'  \code{x}, \code{y} and matrix \code{z} of size \eqn{n\times n}{n*n} (passed
#'  via \code{kde2d.args}).}
#' \item{h}{Bandwidth passed to \code{kde2d}.}
#' @export
#' @family visualizations
#' @author George G. Vega Yon
#' @examples
#'
#' # Example with a random graph --------------------------------------------------
#'
#' set.seed(1231)
#'
#' # Random small-world diffusion network
#' x <- netdiffuseR::rdiffnet(300, 20, seed.graph="small-world",
#'  seed.nodes = "random")
#'
#' # Diffusion map (no random toa)
#' dm0 <- diffusionMap(x, kde2d.args=list(n=100, h=4))
#'
#' # Random
#' diffnet.toa(x) <- (x$toa)[order(runif(300))]
#'
#' # Diffusion map (random toa)
#' dm1 <- diffusionMap(x, layout = dm0$coords, kde2d.args=list(n=100, h=4))
#'
#' oldpar <- par(no.readonly = TRUE)
#' col <- adjustcolor(
#'  colorRampPalette(c("white","lightblue", "yellow", "red"))(100),.5)
#' par(mfrow=c(1,2), oma=c(1,0,0,0))
#' image(dm0, col=col, main="Non-random Times of Adoption")
#' image(dm1, col=col, main="Random Times of Adoption")
#' par(mfrow=c(1,1))
#' mtext("Both networks have the same distribution on times of adoption", 1,
#'  outer = TRUE)
#' par(oldpar)
#'
#' # Example with Brazilian Farmers --------------------------------------------
#'
#' dn <- brfarmersDiffNet
#'
#' # Setting last TOA as NA
#' diffnet.toa(dn)[dn$toa == max(dn$toa)] <-
#'   NA
#'
#' # Coordinates
#' coords <- sna::gplot.layout.fruchtermanreingold(
#'   as.matrix(dn$graph[[1]]), layout.par=NULL
#' )
#'
#' # Plotting diffusion
#' plot_diffnet2(dn, layout=coords, vertex.size = 300)
#'
#' # Adding diffusion map
#' out <- diffusionMap(dn, layout=coords, kde2d.args=list(n=100, h=50))
#' col <- adjustcolor(colorRampPalette(c("white","lightblue", "yellow", "red"))(100),.5)
#' with(out$map, .filled.contour(x,y,z,pretty(range(z), 100),col))
#'
diffusionMap <- function(...) UseMethod("diffusionMap")

#' @export
#' @rdname diffusionMap
diffmap <- diffusionMap

#' @export
#' @rdname diffusionMap
diffusionMap.default <- function(
  graph, x, x.adj=round_to_seq, layout=NULL,
  jitter.args = list(),
  kde2d.args  = list(n=100)) {

  # Step 0) Preparing the data
  if (length(x.adj)) {
    if (!is.function(x.adj)) stop('-x.adj- must be a function')
    x <- x.adj(x)
  }

  ntimes <- as.integer(x)
  ntimes <- ntimes - min(ntimes, na.rm = TRUE) + 1L
  ntimes[is.na(ntimes)] <- 0

  # Computing positions
  g <- igraph::graph_from_adjacency_matrix(graph)
  coords <- if (is.function(layout)) layout(g)
  else if (!length(layout)) igraph::layout_nicely(g)
  else if (is.matrix(layout)) layout

  # # Step 1) Expand using toa as weights
  # n <- length(ntimes)
  # i <- seq_len(n)
  # Coords <- cbind(
  #   rep(coords[,1], ntimes),
  #   rep(coords[,2], ntimes)
  # )
  #
  # # Step 2) Jitter
  # Coords[,1] <- do.call(jitter, c(jitter.args, list(x=Coords[,1])))
  # Coords[,2] <- do.call(jitter, c(jitter.args, list(x=Coords[,2])))
  #
  # # Step 3) kde2d
  #
  # if (!length(kde2d.args$h))
  #   kde2d.args$h <- c(MASS::bandwidth.nrd(Coords[,1]), MASS::bandwidth.nrd(Coords[,2]))

  # Step 1) Compute densities per level
  if (!length(kde2d.args$h))
    kde2d.args$h <- c(MASS::bandwidth.nrd(coords[,1]), MASS::bandwidth.nrd(coords[,2]))

  # Mapping limits
  lims  <- c(range(coords[,1]), range(coords[,2]))
  Map   <- with(kde2d.args, list(z=matrix(0, ncol=n, nrow=n)))
  Map$W <- Map$z
  for (i in unique(ntimes)) {
    # Subset and map
    dat <- coords[ntimes==i,,drop=FALSE]
    map <- do.call(MASS::kde2d, c(kde2d.args, list(
      x = dat[,1], y=dat[,2], lims=lims)))

    # Adding up (for weighted average)
    Map$W <- Map$W + map$z
    Map$z <- Map$z + map$z*i

  }

  # Normalizing
  Map$z <- Map$z/(Map$W + 1e-15)
  Map$x <- seq(lims[1], lims[2], length.out = kde2d.args$n)
  Map$y <- seq(lims[3], lims[4], length.out = kde2d.args$n)

  structure(list(
    coords = coords,
    map    = with(Map, list(x=x,y=y,z=z)),
    h       = kde2d.args$h,
    used_x  = x
  ), class="diffnet_diffmap")
}

#' @rdname diffusionMap
#' @export
diffusionMap.diffnet <- function(graph, slice=nslices(graph), ...) {
  with(graph, diffusionMap.default(graph[[slice]], toa, ...))
}

#' @rdname diffusionMap
#' @export
image.diffnet_diffmap <- function(x, ...) {
  image(x$map,...)
}

#' @rdname diffusionMap
#' @export
print.diffnet_diffmap <- function(x, ...) {
  cat("An object of class -diffnet_map-\n")
  cat(str(x))
  cat("Use methods -plot- and -image-.")
}

#' @rdname diffusionMap
#' @export
plot.diffnet_diffmap <- function(x, y=NULL, ...) {
  image.diffnet_diffmap(x, ...)
}
