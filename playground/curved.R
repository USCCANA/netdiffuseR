rotate <- function(mat, p0, alpha) {
  R <- matrix(
    c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)),
    nrow = 2, byrow = TRUE)

  p0 <- matrix(c(p0[1], p0[2]), ncol=2, nrow = nrow(mat), byrow = TRUE)
  t(R %*% t(mat - p0)) + p0
}

#' @param p0,p1 Numeric vector of length 2. Center coordinates
#' @param alpha Numeric scalar. Arc angle in radians.
#' @param n Integer scalar. Number of segments to approximate the arc.
#' @param radii Numeric vector of length 2. Radious
#' @export
#'
arc <- function(p0, p1, alpha = pi/2, n=20, radii = c(0, 0)) {

  # If no curve, nothing to do (old fashioned straight line)
  if (alpha == 0)
    return(rbind(p0, p1))

  alpha0 <- atan2(p1[2]-p0[2], p1[1] - p0[1])

  # Constants
  d <- stats::dist(rbind(p0, p1))

  # if ((d - sum(radii)) < 0) {
  #   alpha <- 2*pi #2*pi - alpha
  # }

  # if (d > 0)
    r <- d/2/sin(alpha/2)
  # else
  #   r <- max(radii)*1.2

  # Center
  M <- cbind(
      p0[1] + d/2,
      p0[2] - cos(alpha/2)*r
    )

  # Angle range
  alpha_start <- 2*asin(radii[1]/2/r)
  alpha_end   <- 2*asin(radii[2]/2/r)

  alpha_i <- seq(
    alpha - alpha_start ,
    alpha_end,
    length.out = n
  ) + (pi - alpha)/2

  ans <- cbind(
    M[1] + cos(alpha_i)*r,
    M[2] + sin(alpha_i)*r
  )


  # Rotation and return
  ans <- rotate(ans, p0, alpha0)
  structure(
    ans,
    alpha0 = atan2(ans[1,2] - p0[2], ans[1,1] - p0[1]),
    alpha1 = atan2(p1[2] - ans[n,2], p1[1] - ans[n,1])
  )

}

#' Edges coordinates
#' @export
edge <- function(p0, p1, radii =c(0.25, 0.25), alpha = pi/4, n = 20) {

  # Creating curve
  arc(
    p0    = p0,
    p1    = p1,
    alpha = alpha,
    n     = n,
    radii = radii
    )

}


#' @param x Numeric vector of length 2. Coordinates of the tip
#'
arrow_fancy <- function(x, a0 = 0, l=.25, a=pi/6, b = pi/1.5) {


  p_top   <- c(x[1], x[2])
  p_left  <- p_top + c(-cos(a), sin(a))*l

  base <- l*sin(a)
  base2 <- base * cos(pi - b)/sin(pi - b)
  p_mid   <- p_left + c(base2, -base)

  p_right <- p_top - c(cos(a), sin(a))*l

  ans <- rbind(p_top, p_left, p_mid, p_right)

  # Rotation
  p0 <- matrix(p_top, nrow = nrow(ans), ncol=2, byrow = TRUE)
  R  <- matrix(c(cos(a0), -sin(a0), sin(a0), cos(a0)), nrow = 2, byrow = TRUE)
  t(R %*% t(ans - p0)) + p0


}

#' Rescale the size of a node to make it relative to the aspect ratio of the device
#' @param size Numeric vector. Size of the node (radious).
#' @param rel Numeric vector of length 2. Relative size for the minimum and maximum
#' of the plot.
#'
#' @details
#' This function is to be called after [plot.new], as it takes the parameter `usr`
#' from the
rescale_node <- function(size, rel=c(.01, .05)) {

  # Rescaling to be between range[1], range[2]
  sran <- range(size, na.rm=TRUE)

  if ((sran[2] - sran[1]) > 1e-10)
    size <- (size - sran[1])/(sran[2] - sran[1]) # 0-1
  else
    size <- size/sran[1]

  size <- size * (rel[2] - rel[1]) + rel[1]

  # Getting coords
  usr <- graphics::par()$usr[1:2]
  size * (usr[2] - usr[1])/2

}

#' @param width Numeric vector. width of the edges
rescale_edge <- function(width, rel=c(1, 3)) {

  ran   <- range(width, na.rm = TRUE)
  if (ran[1] != ran[2])
    width <- (width - ran[1])/(ran[2] - ran[1])
  else
    width <- width/ran[1]

  width*(rel[2] - rel[1]) + rel[1]

}



library(polygons)

#' Adjust coordinates to fit aspect ratio of the device
#' @param coords Two column numeric matrix. Vertices coordinates.
#' @details
#' It first adjusts `coords` to range between `-1,1`, and then, using
#' `graphics::par("pin")`, it rescales the second column of it (`y`) to adjust
#' for the device's aspec ratio.
fit_coords_to_dev <- function(coords, adj = graphics::par("pin")[1:2]) {

  # Making it -1 to 1
  yran <- range(coords[,2], na.rm = TRUE)
  xran <- range(coords[,1], na.rm = TRUE)

  coords[,1] <- (coords[,1] - xran[1])/(xran[2] - xran[1])*2 - 1
  coords[,2] <- (coords[,2] - yran[1])/(yran[2] - yran[1])*2 - 1

  # Adjusting aspect ratio according to the ploting area
  coords[,2] <- coords[,2]*adj[2]/adj[1]

  # Returning new coordinates
  coords

}


# cols <- viridis::cividis(max(igraph::degree(x) + 1))[igraph::degree(x) + 1]
# cols <- heat.colors(max(igraph::degree(x) + 1))[igraph::degree(x) + 1]
# cols <- igraph::V(x)$toa
# cols <- cols - min(cols, na.rm = TRUE) + 1
# cols[is.na(cols)] <- max(cols, na.rm = TRUE) + 1
# cols <- viridis::cividis(max(cols))[cols]
# clus <- igraph::membership(igraph::cluster_walktrap(x))
# cols <- viridis::cividis(max(clus))[clus]
# cols <- colorRampPalette(c("steelblue", "white"), alpha=1)(max(igraph::degree(x)))



#' A wrapper of `rgb(colorRamp)`
#' @param i,j Integer scalar. Indices of ego and alter from 1 through n.
#' @param p Numeric scalar from 0 to 1. Proportion of mixing.
#' @param alpha Numeric scalar from 0 to 1. Passed to [grDevices:colorRamp]
#' @return A color.
edge_color_mixer <- function(i, j, vcols, p = .5, alpha = .15) {

  grDevices::rgb(
    grDevices::colorRamp(vcols[i], alpha = alpha)(p),
    maxColorValue = 255
  )

}

nplot <- function(
  x,
  layout       = igraph::layout_with_fr(x),
  vertex.size  = igraph::degree(x),
  bg.col       = "black",
  vertex.shape = rep(100, igraph::vcount(x)),
  vertex.color = NULL
  ) {

  if (!length(vertex.color)) {
    vertex.color <- length(table(igraph::degree(x)))
    vertex.color <- topo.colors(vertex.color)
    vertex.color <- vertex.color[
      as.factor(igraph::degree(x))
      ]
  }


  # Creating the window
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Setting margin
  par(mai=rep(.25, 4))

  # Adjusting layout to fit the device
  layout <- fit_coords_to_dev(layout)

  # Plotting
  plot(layout, type = "n", bty="n", xaxt="n", yaxt="n", ylab="", xlab="", asp=1)

  # Adding rectangle
  usr <- graphics::par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = bg.col)

  # Rescaling size
  vertex.size <- rescale_node(vertex.size)

  # Arrow size
  l <- (max(layout[,1]) - min(layout[,1]))/150

  # Weights
  W <- igraph::E(x)$weight
  if (!length(W))
    W <- rep(1.0, length(igraph::E(x)))

  W <- rescale_edge(W/max(W, na.rm=TRUE))

  # Computing shapes -----------------------------------------------------------
  E <- igraph::as_edgelist(x, names = FALSE)

  ans <- vector("list", nrow(E))
  for (e in 1:nrow(E)) {

    i <- E[e,1]
    j <- E[e,2]

    ans[[e]] <- arc(layout[i,], layout[j,], radii = vertex.size[c(i,j)])

  }

  # Edges
  for (i in seq_along(ans)) {

    # Not plotting self (for now)
    if (E[i,1] == E[i, 2])
      next

    # Computing edge color
    col <- edge_color_mixer(
      E[i, 1], E[i, 2],
      vertex.color, alpha = W[i])

    # Drawing lines
    lines(ans[[i]], lwd= W[i], col = col)

    # Computing arrow
    arr <- arrow_fancy(
      ans[[i]][nrow(ans[[i]]),1:2],
      a0 = attr(ans[[i]], "alpha1"),
      l = l,
      b = 2
      )

    # Drawing arrows
    polygon(arr, col = col, border=col)


  }

  # Nodes
  for (i in 1:nrow(layout))
    polygon(
      npolygon(
        layout[i,1], layout[i,2],
        n = vertex.shape[i],
        r = vertex.size[i],
        FALSE
        ),
      col    = vertex.color[i],
      border = adjustcolor(vertex.color[i], red.f = .5, blue.f = .5, green.f = .5),
      lwd=4
      )

}

library(netdiffuseR)
nplot(netdiffuseR::diffnet_to_igraph(medInnovationsDiffNet)[[1]])
# nplot(netdiffuseR::diffnet_to_igraph(kfamilyDiffNet)[[1]])

set.seed(1)
x <- igraph::barabasi.game(120, m = 1, power = .95)
x <- igraph::rewire(x, igraph::each_edge(.075))
igraph::E(x)$weight <- runif(igraph::ecount(x))
nplot(x)
