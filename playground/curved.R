rotate <- function(mat, p0, alpha) {
  R <- matrix(
    c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)),
    nrow = 2, byrow = TRUE)

  p0 <- matrix(c(p0[1], p0[2]), ncol=2, nrow = nrow(mat), byrow = TRUE)
  t(R %*% t(mat - p0)) + p0
}

curves <- function(p0, p1, alpha = -pi/2, n=4, sizes = c(0, 0)) {

  # If no curve, nothing to do (old fashioned straight line)
  if (alpha == 0)
    return(rbind(p0, p1))

  alpha0 <- atan2(p1[2]-p0[2], p1[1] - p0[1])

  # Constants
  d <- stats::dist(rbind(p0, p1))
  r <- d/2/sin(alpha/2)

  # Center
  M <- cbind(
    p0[1] + d/2,
    p0[2] - cos(alpha/2)*r
  )

  # Angle range
  alpha_start <- 2*asin(sizes[1]/2/r)
  alpha_end   <- 2*asin(sizes[2]/2/r)

  alpha_i <- seq(
    alpha - alpha_start ,
    alpha_end,
    length.out = n
  ) + (pi - alpha)/2
  # beta <- (pi - alpha)/2.0

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
edge <- function(x0, y0, x1, y1, s0 = 0.25, s1 = 0.25, alpha = pi/4, n = 20) {

  # Creating curve
  curves(
    p0    = c(x0, y0),
    p1    = c(x1, y1),
    alpha = alpha,
    n     = n,
    sizes = c(s0, s1)
    )

}


arrow_fancy <- function(x0, y0, a0 = 0, l=.25, a=pi/6, b = pi/1.5) {


  p_top   <- c(x0, y0)
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

set.seed(1)
x <- igraph::barabasi.game(200, m = 1, power = .95)
x <- igraph::rewire(x, igraph::each_edge(.075))
# data(brfarmersDiffNet, package="netdiffuseR")
# x <- netdiffuseR::diffnet_to_igraph(brfarmersDiffNet)[[1]]
# x <- igraph::erdos.renyi.game(200, .1)
# x <- igraph::sample_smallworld(1, 40, 4, 0.025)
# x <- readr::read_csv("~/Downloads/edges_2008.csv")
# x <- igraph::graph_from_data_frame(x)
N <- cbind(igraph::layout_with_fr(x), igraph::degree(x))
E <- cbind(igraph::as_edgelist(x, names = FALSE), igraph::E(x)$weight)
E <- cbind(E,sample(c(.75, 1, 2), size = nrow(E), replace = TRUE, prob = c(.75,.1,.05)))

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

par(mai=rep(.25, 4))

N[,1:2] <- fit_coords_to_dev(N[,1:2])

plot(N[,1:2], type = "n", bty="n", xaxt="n", yaxt="n", ylab="", xlab="", asp=1)

rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], col = "black")
N[,3] <- rescale_node(N[,3])

ans <- vector("list", nrow(E))
for (e in 1:nrow(E)) {

  i <- E[e,1]
  j <- E[e,2]

  x0 <- N[i, 1]
  y0 <- N[i, 2]
  x1 <- N[j, 1]
  y1 <- N[j, 2]
  s0 <- N[i, 3]
  s1 <- N[j, 3]
  ans[[e]] <- edge(x0, y0, x1, y1, s0 = s0, s1 = s1)

}


# cols <- viridis::cividis(max(igraph::degree(x) + 1))[igraph::degree(x) + 1]
cols <- heat.colors(max(igraph::degree(x) + 1))[igraph::degree(x) + 1]
# cols <- igraph::V(x)$toa
# cols <- cols - min(cols, na.rm = TRUE) + 1
# cols[is.na(cols)] <- max(cols, na.rm = TRUE) + 1
# cols <- viridis::cividis(max(cols))[cols]
# clus <- igraph::membership(igraph::cluster_walktrap(x))
# cols <- viridis::cividis(max(clus))[clus]
# cols <- colorRampPalette(c("steelblue", "white"), alpha=1)(max(igraph::degree(x)))

l <- (max(N[,1]) - min(N[,1]))/150
D <- igraph::degree(x)

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



for (i in seq_along(ans)) {

  col <- edge_color_mixer(E[i,1], E[i, 2], cols,alpha = E[i,3]/max(E[,3]))

  lines(ans[[i]], lwd=1.5, col = col)

  arr <- arrow_fancy(ans[[i]][20,1], ans[[i]][20,2], a0 = attr(ans[[i]], "alpha1"), l = l, b = 2)
  # col <- adjustcolor("darkgray", 1)
  polygon(arr, col = col, border=col)


}

shapes <- rep(100, nrow(N))

for (i in 1:nrow(N))
  polygon(npolygon(N[i,1], N[i,2], n=shapes[i],r= N[i, 3], FALSE),
          col = cols[i],
          border = adjustcolor(cols[i], red.f = .5, blue.f = .5, green.f = .5),
          lwd=4)

# plot(x, vertex.size=sqrt(igraph::degree(x))*4, vertex.label=NA, edge.curved=TRUE, edge.arrow.size=.5, layout = N[,1:2])


