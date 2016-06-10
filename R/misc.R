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
#' # Simple example
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

#' Draw a color key in the current device
#'
#' @param x A numeric vector with the data (it is used to extract the range).
#' @param main Character scalar. Title of the key.
#' @param tick.marks A numeric vector indicating the levels to be included in the axis.
#' @param key.pos A numeric vector of length 4 with relative coordinates of the
#'  key (as \% of the plotting area, see  \code{\link[base:par]{par("usr")}})
#' @param pos Integer scalar. Position of the axis as in \code{\link[graphics:text]{text}}.
#' @param nlevels Integer scalar. Number of levels (colors) to include in the color key.
#' @param color.palette Color palette of \code{length(nlevels)}.
#' @param tick.with Numeric vector of length 2 indicating the length of the inner
#'  and outer tick marks as percentage of the axis.
#' @param add.box Logical scalar. When \code{TRUE} adds a box around the key.
#' @param ... Further arguments to be passed to \code{\link{graphics::rect}{rect}}
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
drawColorKey <- function(
  x, tick.marks = {
    tick <- pretty(x,min.n = 5)
    r    <- range(x, na.rm = TRUE)
    tick[(tick >= r[1]) & (tick <= r[2])]
    }, main=NULL,
  key.pos=c(.925,.975,.05,.95),
  pos = 2, nlevels=length(tick.marks),
  color.palette=grDevices::colorRampPalette(c("lightblue", "yellow", "red"))(nlevels),
  tick.width=c(.01,.0075), add.box=TRUE, ...) {
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
  if (!add.box) graphics::segments(coords[i], coords[3], coords[i], coords[4])
  graphics::segments(coords[i] - tw[1]*sgn, atick.marks, coords[i] + tw[2]*sgn, atick.marks)
  graphics::text(coords[i] - sgn*tw[1], atick.marks, labels = tick.marks, pos=pos)

  # Adding box
  if (add.box) graphics::rect(coords[1], coords[3], coords[2], coords[4])

  if (length(main))
    graphics::text(
      x=(coords[2] + coords[1])/2,
      y=coords[4], pos=3, labels=main
      )

  invisible(NULL)
}

