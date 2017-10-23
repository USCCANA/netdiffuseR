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
#'
#' @seealso Used in \code{\link{diffmap}} and \code{\link{plot_diffnet2}}
round_to_seq <- function(x, nlevels=20, as_factor=FALSE) {
  y <- range(x, na.rm = TRUE, finite=TRUE)
  y <- seq(y[1], y[2], length.out = nlevels)

  y <- sapply(x, function(z) {
    if (is.na(z)) return(NA)
    y[which.min(abs(y-z))]
    })
  # factor(c(1,3), levels = 1:3, labels = letters[1:3])
  if (as_factor) as.factor(y)
  else y
}

#' Another way of visualizing diffusion
#' @templateVar toa TRUE
#' @templateVar slice TRUE
#' @template graph_template
#' @template plotting_template
#' @param pers Integer vector of length \eqn{T} indicating the time periods of the data.
#' @param color.ramp A function as returned by \code{\link[grDevices:colorRamp]{colorRamp}}.
#' @param layout Passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param key.width Numeric scalar. Sets the proportion of the plot (x-axis) that the key uses.
#' @param key.args List. Further arguments to be passed to \code{\link{drawColorKey}}.
#' @param main Character scalar. Title of the graph.
#' @param add.map Character scalar. When \code{"first"} plots a \code{\link{diffusionMap}} before the
#'  graph itself. If \code{"last"} then it adds it at the end. When \code{NULL} adds nothing.
#' @param diffmap.args List. If \code{add.map=TRUE}, arguments passed to \code{diffusionMap}.
#' @param diffmap.alpha Numeric scalar between [0,1]. Alpha level for the map.
#' @param include.white Character scalar. Includes white in the color palette used in the map.
#'  When \code{include.white=NULL} then it won't include it.
#' @param ... Further arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @param no.graph Logical scala. When \code{TRUE} the graph is not drawn. This only makes
#' sense when the option \code{add.map} is active.
#' @details If \code{key.width<=0} then no key is created.
#'
#' By defult, the function passes the following values to \code{plot.igraph}:
#'
#' \itemize{
#' \item{\code{vertex.label} equals to \code{""}}
#' \item{\code{vertex.frame.color} equals to \code{"white"}}
#' \item{\code{add} equals to \code{TRUE}}
#' \item{\code{rescale} equals to \code{FALSE}}
#' \item{\code{vertex.size} equals to \code{rescale.fun(vertex.size)}}
#' }
#'
#' @return A list with the following elements
#' \item{layout}{A numeric matrix with vertex coordinates.}
#' \item{vertex.color}{A character vector with computed colors for each vertex.}
#' \item{vertex.label}{The value passed to \code{plot_diffnet2}.}
#' \item{vertex.shape}{A character vector with assigned shapes.}
#' \item{vertex.size}{A numeric vector with vertices sizes}
#' \item{diffmap}{If \code{add.map=TRUE}, the returned values from \code{\link{diffmap}}}
#' @export
#' @family visualizations
#' @author George G. Vega Yon
plot_diffnet2 <- function(graph, ...) UseMethod("plot_diffnet2")

#' @rdname plot_diffnet2
#' @export
#' @include diffnet-methods.r data.r
plot_diffnet2.diffnet <- function(
  graph,
  toa,
  slice = nslices(graph),
  ...
) {

  if (missing(toa))
    toa <- graph$toa

  plot_diffnet2.default(
    graph         = graph$graph[[slice]],
    toa           = toa,
    pers          = graph$meta$pers,
    ...)
}

#' @rdname plot_diffnet2
#' @export
plot_diffnet2.default <- function(
  graph,
  toa,
  pers          = min(toa, na.rm = TRUE):max(toa, na.rm = TRUE),
  color.ramp    = grDevices::colorRamp(c("steelblue","gray", "tomato")),
  layout        = NULL,
  key.width     = 0.1,
  key.args      = list(),
  main          = "Diffusion dynamics",
  add.map       = NULL,
  diffmap.args  = list(kde2d.args=list(n=100)),
  diffmap.alpha = .5,
  include.white = "first",
  vertex.size   = "degree",
  minmax.relative.size = getOption("diffnet.minmax.relative.size", c(0.01, 0.04)),
  no.graph      = FALSE,
  ...) {

  # Modifying some arguments
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  par(xpd = NA)

  # Collecting arguments
  igraph.args <- list(...)

  # Some constants
  nper <- length(pers)

  if (length(add.map) && !(add.map %in% c("first", "last")))
    stop("When -add.map- is specified it should be either \'first\' or \'last\'.")

  if (!length(add.map) & no.graph)
    stop("If -no.graph=TRUE- then you should specify some value for -add.map-.")

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

  # Adding alpha
  if (ncol(col) < 4)
    col <- cbind(col, 255)

  col[type_non,] <- 255

  col <- rgb(col[,1], col[,2], col[,3], col[,4], maxColorValue = 255)

  # Shapes
  if (!no.graph && !length(igraph.args$vertex.shape)) {
    igraph.args$vertex.shape <- rep("circle", nnodes(graph))
    igraph.args$vertex.shape[type_non] <- "square"
  }

  # Adjmat must have dimnames to make sure sorting in igraph is fine
  add_dimnames.mat(graph)

  # Computing positions
  g <- igraph::graph_from_adjacency_matrix(graph, weighted = TRUE)

  igraph.args$layout <- if (!length(layout)) igraph::layout_nicely(g)
  else if (inherits(layout, "function")) layout(g)
  else layout

  # Keywidth
  key.width <- max(0, key.width)

  graphics::plot.new()
  graphics::plot.window(xlim=c(-1,1 + 5*key.width), ylim=c(-1,1))
  graphics::title(main=main)

  # If adding map! -------------------------------------------------------------
  if (length(add.map)) {

    dm <- do.call(diffusionMap.default, c(diffmap.args, list(graph=graph, x=toa,
                                                             layout = igraph.args$layout)))
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
      graphics::.filled.contour(dm$map$x, dm$map$y, dm$map$z, levels = dmlvls, col=dmcol)
  } else dm <- NULL

  # Plotting graph -------------------------------------------------------------

  # Setting up parameters
  set_igraph_plotting_defaults("igraph.args")

  igraph.args$vertex.size <- rescale_vertex_igraph(
    compute_vertex_size(g, vertex.size),
    minmax.relative.size = minmax.relative.size
    )

  igraph.args$vertex.color <- col

  # Calling igraph
  if (!no.graph)
    do.call(
      what = igraph::plot.igraph,
      args = c(list(g),igraph.args)
    )


  if (length(add.map) && (add.map=="last"))
      graphics::.filled.contour(dm$map$x, dm$map$y, dm$map$z, levels = dmlvls, col=dmcol)

  # # Plotting boxes -------------------------------------------------------------
  if (key.width > 0) {
    # Adjusting the color
    color.palette <- color.ramp(c(0,.5,1))

    if (ncol(color.palette) < 4)
      color.palette <- cbind(color.palette, 255)

    color.palette <- grDevices::rgb(
      color.palette[,1], color.palette[,2], color.palette[,3],
      color.palette[,4],
      maxColorValue = 255)

    color.palette <- grDevices::colorRampPalette(color.palette, TRUE)

    # Filling missings
    if (!length(key.args$main))   key.args$main <- "Time of Adoption"
    if (!length(key.args$na.col)) key.args$na.col <- "transparent"
    if (!length(key.args$na.lab)) key.args$na.lab <- "Non-adopters"
    if (!length(key.args$border)) key.args$border <- "transparent"
    if (!length(key.args$tick.marks)) {
      toaran <- range(toa, na.rm=TRUE)
      key.args$tick.marks <-
        unique(floor(seq(toaran[1], toaran[2], length.out = 5)))
    }



    do.call(
      what = drawColorKey,
      args = c(
        list(toa, key.pos = c(1-key.width, 0.975, 0.05, 0.95), nlevels = 100,
             color.palette = color.palette(100)),
        key.args
      )
    )
  }


  invisible(list(
    layout       = igraph.args$layout,
    vertex.color = col,
    vertex.size  = igraph.args$vertex.size,
    vertex.shape = igraph.args$vertex.shape,
    diffmap      = dm)
    )
}


#' Creates a heatmap based on a graph layout and a vertex attribute
#'
#' Using bi-dimensional kernel smoothers, creates a heatmap based on a graph layout
#' and colored accordingly to \code{x}. This visualization technique is intended
#' to be used with large graphs.
#'
#' @param graph A square matrix of size \eqn{n\times n}{n * n}.
#' @param slice Integer scalar. Slice of the network to be used as baseline for drawing the graph.
#' @param x An vector of length \eqn{n}. Usually a \code{toa} vector.
#' @param layout Either a \eqn{n\times 2}{n *2} matrix of coordinates or a layout
#'  function applied to \code{graph} (must return coordinates).
#' @param jitter.args A list including arguments to be passed to \code{\link{jitter}}.
#' @param kde2d.args A list including arguments to be passed to \code{\link[MASS:kde2d]{kde2d}}.
#' @param sharp.criter A function choose whether to apply a weighted mean for each cell,
#' or randomize over the values present in that cell (see details).
#' @param ... Arguments passed to method.
#' @details
#' The image is created using the function \code{kde2d} from
#' the \pkg{MASS} package. The complete algorithm follows:
#' \enumerate{
#'  \item \code{x} is coerced into integer and the range is adjusted to start from 1.
#'    \code{NA} are replaced by zero.
#'  \item If no \code{layout} is passed, layout is computed using
#'    \code{\link[igraph:layout_nicely]{layout_nicely}} from \pkg{igraph}
#'  \item Then, a \code{kde2d} map is computed for each level of \code{x}. The
#'    resulting matrices are added up as a weighted sum. This only holds if
#'    at the cell level the function \code{sharp.criter} returns \code{FALSE}.
#'  \item The jitter function is applied to the repeated coordinates.
#'  \item 2D kernel is computed using \code{kde2d} over the coordinates.
#' }
#'
#' The function \code{sharp.criter} must take two values, a vector of levels and a
#' vector of weights. It must return a logical scalar with value equal to \code{TRUE}
#' when a randomization at the cell level must be done, in which case the final
#' value of the cell is chosen using \code{sample(x, 1, prob=w)}.
#'
#' The resulting matrix can be passed to \code{\link{image}} or similar.
#'
#' The argument \code{x.adj} uses by default the function \code{\link{round_to_seq}}
#' which basically maps \code{x} to a fix length sequence of numbers such that
#' \code{x.adj(x)} resembles an integer sequence.
#'
#' @return A list of class \code{diffnet_diffmap}
#' \item{coords}{A matrix of size \eqn{n\times 2}{n*2} of vertices coordinates.}
#' \item{map}{Output from \code{kde2d}. This is a list with 3 elements, vectors
#'  \code{x}, \code{y} and matrix \code{z} of size \eqn{n\times n}{n*n} (passed
#'  via \code{kde2d.args}).}
#' \item{h}{Bandwidth passed to \code{kde2d}.}
#' @export
#' @family visualizations
#' @references Vega Yon, George G., and Valente, Thomas W., Visualizing Large Annotated
#' Networks as Heatmaps using Weighted Averages based on Kernel Smoothers (Working paper).
#' @author George G. Vega Yon
#' @examples
#'
#' # Example with a random graph --------------------------------------------------
#'
#' \dontrun{
#' set.seed(1231)
#'
#' # Random scale-free diffusion network
#' x <- rdiffnet(1000, 4, seed.graph="scale-free", seed.p.adopt = .025,
#'                            rewire = FALSE, seed.nodes = "central",
#'                            rgraph.arg=list(self=FALSE, m=4),
#'                            threshold.dist = function(id) runif(1,.2,.4))
#'
#' # Diffusion map (no random toa)
#' dm0 <- diffusionMap(x, kde2d.args=list(n=150, h=.5), layout=igraph::layout_with_fr)
#'
#' # Random
#' diffnet.toa(x) <- sample(x$toa, size = nnodes(x))
#'
#' # Diffusion map (random toa)
#' dm1 <- diffusionMap(x, layout = dm0$coords, kde2d.args=list(n=150, h=.5))
#'
#' oldpar <- par(no.readonly = TRUE)
#' col <- colorRampPalette(blues9)(100)
#' par(mfrow=c(1,2), oma=c(1,0,0,0))
#' image(dm0, col=col, main="Non-random Times of Adoption\nAdoption from the core.")
#' image(dm1, col=col, main="Random Times of Adoption")
#' par(mfrow=c(1,1))
#' mtext("Both networks have the same distribution on times of adoption", 1,
#'       outer = TRUE)
#' par(oldpar)
#' }
#'
#' # Example with Brazilian Farmers --------------------------------------------
#' \dontrun{
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
#' }
#'
diffusionMap <- function(graph, ...) UseMethod("diffusionMap")

#' @export
#' @rdname diffusionMap
diffmap <- diffusionMap

#' Computes weighted variance
#' @param x A numeric vector of length \eqn{n}.
#' @param w A numeric vector of length \eqn{n}.
#' @details \code{weighted_variance} implements weighted variance computation
#' in the following form:
#' \deqn{%
#' \frac{\sum_i w_i'(x_i - \bar x)^2}{(1-n)}
#' }{%
#' sum[w(i)'(x(i) - w.mean(x))^2/(1-n)]
#' }
#'
#' where \eqn{w_i'=w_i/\sum_i w_i}{w(i)' = w(i)/sum(w)}, and
#' \eqn{\bar x = \sum_i w_i'x_i}{w.mean(x)=sum[w(i)'*x(i)]}.
#' @return Numeric scalar with the weighted variance.
#' @export
#' @seealso This function is used in \code{\link{diffmap}}.
weighted_var <- function(x,w) {
  n <- length(x)
  w <- w/sum(w, na.rm=TRUE)*n
  m <- sum(x*w/sum(w, na.rm=TRUE), na.rm=TRUE)
  sum((x - m)^2*w/(n-1+1e-15), na.rm=TRUE)
}

#' @export
#' @rdname weighted_var
wvar <- weighted_var

#' @export
#' @param x.adj Function to adjust \code{x}. If not \code{NULL} then it is applied
#'  to \code{x} at the beginning (see details).
#' @rdname diffusionMap
diffusionMap.default <- function(
  graph, x, x.adj=round_to_seq, layout=NULL,
  jitter.args = list(),
  kde2d.args  = list(n=100),
  sharp.criter=function(x, w) {
    wvar(x,w) > (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))^2/12
    },...) {

  # Step 0) Preparing the data
  if (length(x.adj)) {
    if (!is.function(x.adj)) stop('-x.adj- must be a function')
    x <- x.adj(x)
  }

  # Adjmat must have dimnames to make sure sorting in igraph is fine
  if (!length(unlist(dimnames(graph), recursive = TRUE)))
    dimnames(graph) <- list(1:nnodes(graph), 1:nnodes(graph))


  # Computing positions
  g <- igraph::graph_from_adjacency_matrix(graph, weighted = TRUE)

  coords <- if (is.function(layout)) layout(g)
  else if (!length(layout)) igraph::layout_nicely(g)
  else if (is.matrix(layout)) layout

  # Step 1) Compute densities per level
  if (!length(kde2d.args$h))
    kde2d.args$h <- c(MASS::bandwidth.nrd(coords[,1]), MASS::bandwidth.nrd(coords[,2]))

  # Mapping limits
  lims   <- c(range(coords[,1]), range(coords[,2]))
  lvls   <- unique(x)
  nlvls  <- length(unique(x))
  Maps   <- with(kde2d.args, list(z=array(0, dim=c(n,n,nlvls) )))
  Maps$W <- Maps$z
  for (i in 1:nlvls) {
    # Skip if NA
    if (is.na(lvls[i])) next

    # Subset and map
    dat <- coords[which(x==lvls[i]),,drop=FALSE]
    map <- do.call(MASS::kde2d, c(kde2d.args, list(
      x = dat[,1], y=dat[,2], lims=lims)))

    # Adding up (for weighted average)
    Maps$W[,,i] <- map$z
    Maps$z[,,i] <- map$z*lvls[i]

  }

  # Processing each level
  Map     <- with(kde2d.args, list(z=matrix(0, ncol=n, nrow=n)))
  Map$W   <- Map$z

  for (i in 1:kde2d.args$n)
    for (j in 1:kde2d.args$n) {

      # Computing variance at that level
      if (sharp.criter(lvls,Maps$W[i,j,]) || sum(Maps$W[i,j,]) < 1e-30 ) {
        Map$z[i,j] <- sum(Maps$z[i,j,])/(sum(Maps$W[i,j,]) + 1e-15)
      } else {
        Map$z[i,j] <- sample(lvls, 1, prob=Maps$W[i,j,])
      }
    }

  # Normalizing
  # Map$z <- Map$z/(Map$W + 1e-15)
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
  graphics::image(x$map,...)
}

#' @rdname diffusionMap
#' @export
print.diffnet_diffmap <- function(x, ...) {
  cat("An object of class -diffnet_map-\n")
  cat(utils::str(x))
  cat("Use methods -plot- and -image-.")
}

#' @rdname diffusionMap
#' @param y Ignored.
#' @export
plot.diffnet_diffmap <- function(x, y=NULL, ...) {
  image.diffnet_diffmap(x, ...)
}
