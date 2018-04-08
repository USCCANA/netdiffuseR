#' Edges coordinates
#' @export
edge <- function(x0, y0, x1, y1, s0 = 0.25, s1 = .25, s = .30, curved = TRUE) {

  # Find points
  alpha <- atan2(y1-y0, x1-x0)

  p0 <- c(x0 + s0*cos(alpha), y0 + s0*sin(alpha))
  p1 <- c(x1 - s1*cos(alpha), y1 - s1*sin(alpha))

  d <- sqrt(sum((p0 - p1)^2.0))
  beta <- alpha - pi/2.0
  pmid <- (p1 + p0)/2.0 + s*d*c(cos(beta), sin(beta))

  # Updating the points
  alpha0m <- atan2(pmid[2] - p0[2], pmid[1] - p0[1])
  p0 <- c(x0 + s0*cos(alpha0m), y0 + s0*sin(alpha0m))
  alpha1m <- atan2(p1[2] - pmid[2], p1[1] - pmid[1])
  p1 <- c(x1 - s1*cos(alpha1m), y1 - s1*sin(alpha1m))

  # Adding midpoints
  # alpha0p0 <- atan2(p0[2] - y0, p0[1] - x0)
  # dp0m     <- sqrt(sum(p0 - pmid)^2)/1000.0
  # p0mid    <- c(p0[1] + cos(alpha0p0)*dp0m, p0[2] + sin(alpha0p0)*dp0m)
  p0mid <- p0
  p1mid <- p1

  # Update midpoint
  d <- sqrt(sum((p0mid - p1mid)^2.0))
  beta <- alpha - pi/2.0
  pmid <- (p1mid + p0mid)/2.0 + s*d*c(cos(beta), sin(beta))

  structure(
    rbind(p0, p0mid, pmid, p1mid, p1),
    class = "matrix",
    alpha = atan2(y1 - p1[2], x1 - p1[1])
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
E<- matrix(runif(4*20, 0, 10), ncol=4)
# set.seed(7)
E <- cbind(E, matrix(rbeta(nrow(E)*2, 2, 15), ncol=2))
# E[4,5:6] <- E[1, 5]
# E <- E[-4,]

N <- rbind(E[,c(1:2, 5)], E[,c(3:4, 6)])

# set.seed(1)
# x <- igraph::barabasi.game(20)
# igraph::V(x)$deg <- igraph::degree(x)
# N <- igraph::layout.auto(x)
# E <- igraph::as_long_data_frame(x)

library(polygons)
ans <- vector("list", nrow(E))
for (i in 1:nrow(E)) {

  ans[[i]] <- edge(E[i, 1], E[i, 2], E[i, 3], E[i, 4], s0 = E[i, 5], s1 = E[i, 6])

}


# symbols(N[,1], N[,2], circles = N[,3], inches = FALSE)
# plot(USArrests)
# par(new=TRUE, xpd=NA)
plot(N[,1:2], type = "n", bty="n", asp=1, xaxt="n", yaxt="n")


# done <- vector(length = nrow(N))
for (i in ans) {


  xspline(i, shape = c(0,1,1,1,0), lwd=2, border = "steelblue")

  arr <- arrow_fancy(i["p1",1], i["p1",2], a0 = attr(i, "alpha"))
  polygon(arr, col = "steelblue", border="steelblue")


}


for (i in 1:nrow(N))
  polygon(circle(N[i,1], N[i,2], N[i, 3], FALSE), col = "tomato")
