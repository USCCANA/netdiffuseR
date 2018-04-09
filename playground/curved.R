#' Edges coordinates
#' @export
edge <- function(x0, y0, x1, y1, s0 = 0.25, s1 = 0.25, s = 0.25, curved = TRUE) {

  d <- sqrt(sum((c(x0, y0) - c(x1, y1))^2.0)) - s0 - s1

  p0 <- c(x0 + s0, y0)
  p1 <- p0 + c(d, 0)

  pmid <- (p1 + p0)/2.0 - c(0, d*s)

  # Updating points
  beta <- atan2(pmid[2] - p0[2], pmid[1] - p0[1])
  p0   <- c(x0, y0) + c(cos(beta), sin(beta))*s0
  beta <- atan2(p1[2] - pmid[2], p1[1] - pmid[1])
  p1   <- p1 + c(s1, 0) - c(cos(beta), sin(beta))*s1

  p0mid <- (p0 + pmid)/2
  p1mid <- (p1 + pmid)/2

  # # Find points
  alpha <- atan2(y1-y0, x1-x0)


  # Rotating
  ans <- rbind(p0, p0mid, pmid, p1mid, p1)

  R  <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), nrow = 2, byrow = TRUE)
  p0 <- matrix(c(x0, y0), ncol=2, nrow = nrow(ans), byrow = TRUE)
  ans <- t(R %*% t(ans - p0)) + p0

  beta <- atan2(ans["p1mid", 2] - ans["pmid", 2], ans["p1mid", 1] - ans["pmid", 1])

  structure(
    ans[-3,],
    alpha = beta
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

E <- matrix(
  c(1,1,2,2,2,1,3,0,1,3,0,2, 1, 1, 1, 1), byrow = TRUE, ncol=4
)
E<- matrix(runif(4*10, 0, 5), ncol=4)
# set.seed(7)
E <- cbind(E, matrix(rbeta(nrow(E)*2, 2, 15), ncol=2))
# E[4,5:6] <- E[1, 5]
# E <- E[-4,]

N <- rbind(E[,c(1:2, 5)], E[,c(3:4, 6)])

set.seed(1)
# x <- igraph::barabasi.game(20, m = 2, power = 1.25)
# data(brfarmersDiffNet, package="netdiffuseR")
# x <- netdiffuseR::diffnet_to_igraph(brfarmersDiffNet)[[1]]
# x <- netdiffuseR::diffnet_to_igraph(netdiffuseR::rdiffnet(n = 200, t=10, seed.graph = "small-world", seed.nodes = "central"))[[1]]
x <- readr::read_csv("~/Downloads/edges_2008.csv")
x <- igraph::graph_from_data_frame(x)
N <- cbind(igraph::layout_with_fr(x), (igraph::degree(x)+1)^(1/4)/10)
E <- cbind(igraph::as_edgelist(x, names = FALSE), igraph::E(x)$weight)
# E <- cbind(E,sample(c(.75, 1, 2), size = nrow(E), replace = TRUE, prob = c(.75,.1,.05)))

library(polygons)
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


# symbols(N[,1], N[,2], circles = N[,3], inches = FALSE)
# plot(USArrests)
# par(new=TRUE, xpd=NA)
par(mai=rep(.25, 4))
plot(N[,1:2], type = "n", bty="n", asp=1, xaxt="n", yaxt="n", ylab="", xlab="")

rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], col = "black")

# done <- vector(length = nrow(N))

cols <- viridis::cividis(max(igraph::degree(x) + 1))[igraph::degree(x) + 1]
# cols <- igraph::V(x)$toa
# cols <- cols - min(cols, na.rm = TRUE) + 1
# cols[is.na(cols)] <- max(cols, na.rm = TRUE) + 1
# cols <- viridis::cividis(max(cols))[cols]
# clus <- igraph::membership(igraph::cluster_walktrap(x))
# cols <- viridis::cividis(max(clus))[clus]
# cols <- colorRampPalette(c("steelblue", "white"), alpha=1)(max(igraph::degree(x)))

l <- .05 # (max(N[,1]) - min(N[,1]))/160
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

  xspline(ans[[i]], shape = c(0,1,1,0), lwd=E[i,3]/max(E[,3]),
          border = col)

  arr <- arrow_fancy(ans[[i]]["p1",1], ans[[i]]["p1",2], a0 = attr(ans[[i]], "alpha"), l = l, b = 2)
  # col <- adjustcolor("darkgray", 1)
  polygon(arr, col = col, border=col)


}

for (i in 1:nrow(N))
  polygon(npolygon(N[i,1], N[i,2], n=100,r=N[i, 3], FALSE),
          col = cols[i], border = adjustcolor(cols[i], red.f = .8, blue.f = .8, green.f = .8),
          lwd=1.5)

# plot(x, vertex.size=sqrt(igraph::degree(x))*4, vertex.label=NA, edge.curved=TRUE, edge.arrow.size=.5, layout = N[,1:2])


