// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

arma::vec seq_cpp(double from, double to, int lengthout) {
  arma::vec out(lengthout);
  double step = (to - from)/(lengthout - 1);
  for(int i=0;i<lengthout;i++)
    out(i) = from + i*step;

  return out;

}

//' Distribution over a grid
//'
//' Distribution of pairs over a grid of fix size.
//'
//' @param x Numeric vector of size \eqn{n}
//' @param y Numeric vector of size \eqn{n}
//' @param nlevels Integer scalar. Number of bins to return
//' @details
//'
//' This function ment for internal use only.
//'
//' @export
//' @keywords misc dplot
//' @seealso Used by \code{\link{plot_infectsuscep}}
//' @return Returns a list with three elements
//' \item{x}{Numeric vector of size \code{nlevels} with the class marks for x}
//' \item{y}{Numeric vector of size \code{nlevels} with the class marks for y}
//' \item{z}{Numeric matrix of size \code{nlevels} by \code{nlevels} with the distribution %
//' of the elements in terms of frecuency}
//' @section Examples:
//' \preformatted{
//' # Generating random vectors of size 100
//' x <- rnorm(100)
//' y <- rnorm(100)
//'
//' # Calculating distribution
//' grid_distribution(x,y,20)
//' }
// [[Rcpp::export]]
List grid_distribution(const arma::vec & x, const arma::vec & y, int nlevels=100) {

  // Checking sizes of the vectors
  int m = x.size();
  int s = y.size();

  if (m!=s) stop("x and y don't have the same length.");

  // Crating empty matrix jointly with the sequences
  arma::mat distmat(nlevels,nlevels, arma::fill::zeros);
  double xlim[2], ylim[2];
  xlim[0] = x.min()-1e-10;
  xlim[1] = x.max()+1e-10;
  ylim[0] = y.min()-1e-10;
  ylim[1] = y.max()+1e-10;

  arma::vec xseq = seq_cpp(xlim[0], xlim[1], nlevels + 1);
  arma::vec yseq = seq_cpp(ylim[0], ylim[1], nlevels + 1);

  for(int k=0;k<m;k++)
    for(int i=0;i<nlevels;i++) {
      bool cnt = false;
      for(int j=0;j<nlevels;j++)
        // Testing if x and y are in the range
        if ( ((x(k) <= xseq(i+1)) & (x(k) > xseq(i))) & ((y(k) <= yseq(j+1)) & (y(k) > yseq(j)))) {
          distmat(i,j) += 1;
          cnt=true;
          break;
        }
      if (cnt) break;
    }

  // Output class mark
  NumericVector xmark(nlevels);
  NumericVector ymark(nlevels);

  for(int i=0;i<nlevels;i++)
    xmark[i] = (xseq(i) + yseq(i+1))/2,
      ymark[i] = (yseq(i) + yseq(i+1))/2;

  return List::create(_["x"]=xmark, _["y"]=ymark, _["z"]=distmat);
}


// arma::mat grid_dist(const arma::vec & xran, const arma::vec & yran, const & arma::vec)

/** *R
library(microbenchmark)
microbenchmark(
  seq_cpp(0,1,100),
  seq(0,1,length.out = 100),
  times=1000
)

# Testing
set.seed(123)
x <- rnorm(200000)
y <- rnorm(200000)
z <- grid_distribution(x,y,40)

sum(z$z)
with(z, contour(x,y,z/sum(z), nlevels = 40))
with(z, rgl::persp3d(as.vector(x),as.vector(y),z/sum(z), col="lightblue"))
*/

//' Compute ego/alter edge coordinates considering alter's size and aspect ratio
//'
//' Given a graph, vertices' positions and sizes, calculates the absolute positions
//' of the endpoints of the edges considering the plot's aspect ratio.
//'
//' @param graph A square matrix of size \eqn{n}. Adjacency matrix.
//' @param toa Integer vector of size \eqn{n}. Times of adoption.
//' @param x Numeric vector of size \eqn{n}. x-coordinta of vertices.
//' @param y Numeric vector of size \eqn{n}. y-coordinta of vertices.
//' @param vertex_cex Numeric vector of size \eqn{n}. Vertices' sizes in terms
//' of the x-axis (see \code{\link{symbols}}).
//' @param undirected Logical scalar. Whether the graph is undirected or not.
//' @param no_contemporary Logical scalar. Whether to return (calcular) edges'
//' coordiantes for vertices with the same time of adoption (see details).
//' @param dev Numeric vector of size 2. Height and width of the device (see details).
//' @param ran Numeric vector of size 2. Range of the x and y axis (see details).
//' @return A numeric matrix of size \eqn{m\times 5}{m * 5} with the following
//' columns:
//' \item{x0, y0}{Edge origin}
//' \item{x1, y1}{Edge target}
//' \item{alpha}{Relative angle between \code{(x0,y0)} and \code{(x1,y1)} in terms
//' of radians}
//' With \eqn{m} as the number of resulting edges.
//' @details
//'
//' In order to make the plot's visualization more appealing, this function provides
//' a straight forward way of computing the tips of the edges considering the
//' aspect ratio of the axes range. In particular, the following corrections are
//' made at the moment of calculating the egdes coords:
//'
//' \itemize{
//' \item{Instead of using the actual distance between ego and alter, a relative
//' one is calculated as follows
//' \deqn{d'=\left[(x_0-x_1)^2 + (y_0' - y_1')^2\right]^\frac{1}{2}}{d'=sqrt[(x0-x1)^2 + (y0'-y1')^2]}
//' where \eqn{%
//' y_i'=y_i\times\frac{\max x - \min x}{\max y - \min y} }{%
//' yi' = yi * [max(x) - min(x)]/[max(y) - min(y)]}
//' }
//' \item{Then, for the relative elevation angle, \code{alpha}, the relative distance \eqn{d'}
//' is used, \eqn{\alpha'=\arccos\left( (x_0 - x_1)/d' \right)}{\alpha' = acos[ (x0 - x1)/d' ]}}
//' \item{Finally, the edge's endpoint's (alter) coordinates are computed as follows: %
//' \deqn{%
//'   x_1' = x_1 + \cos(\alpha')\times v_1}{%
//'   x1' = x1 + cos(\alpha') * v1
//' }
//' \deqn{%
//'   y_1' = y_1 -+ \sin(\alpha')\times v_1 \times\frac{\max y - \min y}{\max x - \min x} }{%
//'   y1' = y1 -+ sin(\alpha')*[max(y) - min(y)]/[max(x) - min(x)]
//' }
//' Where \eqn{v_1}{v1} is alter's size in terms of the x-axis, and the sign of
//' the second term in \eqn{y_1'}{y1'} is negative iff \eqn{y_0 < y_1}{y0<y1}.
//' }
//' }
//'
//' The same process (with sign inverted) is applied to the edge starting piont.
//' The resulting values, \eqn{x_1',y_1'}{x1',y1'} can be used with the function
//' \code{\link{arrows}}. This is the workhorse function used in \code{\link{plot_threshold}}.
//'
//' The \code{dev} argument provides a reference to rescale the plot accordingly
//' to the device, and former, considering the size of the margins as well (this
//' can be easily fetched via \code{par("pin")}, plot area in inches).
//'
//' On the other hand, \code{ran} provides a reference for the adjustment
//' according to the range of the data, this is \code{range(x)[2] - range(x)[1]}
//' and \code{range(y)[2] - range(y)[1]} respectively.
//'
//' @keywords misc dplot
//' @examples
//' # --------------------------------------------------------------------------
//' data(medInnovationsDiffNet)
//' library(sna)
//'
//' # Computing coordinates
//' set.seed(79)
//' coords <- sna::gplot(as.matrix(medInnovationsDiffNet$graph[[1]]))
//'
//' # Getting edge coordinates
//' vcex <- rep(1.5, nnodes(medInnovationsDiffNet))
//' ecoords <- edges_coords(
//'   medInnovationsDiffNet$graph[[1]],
//'   diffnet.toa(medInnovationsDiffNet),
//'   x = coords[,1], y = coords[,2],
//'   vertex_cex = vcex,
//'   dev = par("pin")
//'   )
//'
//' ecoords <- as.data.frame(ecoords)
//'
//' # Plotting
//' symbols(coords[,1], coords[,2], circles=vcex,
//'   inches=FALSE, xaxs="i", yaxs="i")
//'
//' with(ecoords, arrows(x0,y0,x1,y1, length=.1))
//' @export
// [[Rcpp::export]]
NumericMatrix edges_coords(
    const arma::sp_mat & graph,
    const arma::colvec & toa,
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & vertex_cex,
    bool undirected=true,
    bool no_contemporary=true,
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
) {

  // The output matrix has the following
  // - x0 and y0
  // - x1 and y1
  // - alpha
  std::vector< double > x0;
  std::vector< double > y0;
  std::vector< double > x1;
  std::vector< double > y1;
  std::vector< double > alpha;

  // Rescaling the vertex sizes
  arma::colvec vertex_size(vertex_cex);

  // If yexpand is too small, just throw an error
  if (ran.length() == 0) {
    ran = NumericVector::create(2);
    ran[0] = x.max() - x.min();
    ran[1] = y.max() - y.min();
  }

  // Expansion factor for y
  double yexpand = 1.0;
  if ( ran[1] > 1e-5 ) yexpand = ran[1]/ran[0];

  // Adjusting for device size
  if (dev.length() == 0)
    dev = NumericVector::create(2,1.0);

  yexpand = yexpand * (dev[0]/dev[1]);

  // The the filled elements of the graph
  arma::umat indexes = sparse_indexes(graph);

  // for(int i=0;i<n;i++) {
  for(unsigned I=0;I<indexes.n_rows;I++) {

    int i = indexes.at(I,0);
    int j = indexes.at(I,1);

    // Checking conditions
    if (undirected && (i < j)) continue;
    if (no_contemporary && (toa(i)==toa(j)) ) continue;

    // Computing angle
    double a = angle(x(i), y(i)/yexpand, x(j), y(j)/yexpand);
    alpha.push_back(a);

    // Adding the xs and the ys
    x0.push_back(x.at(i) + cos(a)*vertex_size.at(i));
    x1.push_back(x.at(j) - cos(a)*vertex_size.at(j));

    // The formula needs an extra help to figure out the ys
    y0.push_back(y.at(i) + sin(a)*vertex_size.at(i)*yexpand);
    y1.push_back(y.at(j) - sin(a)*vertex_size.at(j)*yexpand);
  }

  // Building up the output
  int e = x0.size();
  NumericMatrix out(e,5);
  for(int i=0; i<e; i++) {
    out(i,0) = x0[i];
    out(i,1) = y0[i];
    out(i,2) = x1[i];
    out(i,3) = y1[i];
    out(i,4) = alpha[i];
  }

  colnames(out) = CharacterVector::create("x0", "y0", "x1", "y1", "alpha");

  return out;

}

/** *R
library(netdiffuseR)
set.seed(123)
graph <- rgraph_ba()
toa <- sample(1:5, nnodes(graph), TRUE)
pos <- sna::gplot.layout.random(as.matrix(graph), NULL)
pos[,2] <- pos[,2]*20
cex <- seq(1,3,length.out = nnodes(graph))/40

# Adjusting by device size, mar and mai
arr <- as.data.frame(edges_coords(graph, toa, pos[,1], pos[,2], cex,
                                  dev = par("pin)))

#Fixing sizes and others
xran <- range(pos[,1])
yran <- range(pos[,2])
plot(pos, col="white", xlim= xran, ylim= yran,
     xaxs="i", yaxs="i",xaxt="n", yaxt="n")

yran <- pretty(seq(yran[1], yran[2], length.out = 10), n=10)
xran <- pretty(seq(xran[1], xran[2], length.out = 10), n=10)
axis(1, at =xran)
axis(2, at =yran)

with(arr, arrows(x0, y0, x1, y1))
symbols(pos[,1], pos[,2], circles=cex, add=TRUE, inches = FALSE, bg=rgb(.3,.3,.8,.2),
        xaxs="i", yaxs="i")
text(pos[,1], pos[,2], labels = 1:11)

# Taking a look at one of these
pos[c(4,11),]

a <- atan((-3.4510269 - -0.8881612)/(0.8045981-0.3114116))
*/


// [[Rcpp::export]]
arma::mat edges_arrow(
    const double & x0,
    const double & y0,
    const double & x1,
    const double & y1,
    const double & height,
    const double & width,
    const double beta = 1.5707963267949, // PI/2
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
) {
  // Creating output
  arma::mat coords(6,2);

  // If yexpand is too small, just throw an error ------------------------------
  if (ran.length() == 0) {
    ran = NumericVector::create(2);
    ran[0] = (x1 > x0 ? x1 - x0: x0 -x1);
    ran[1] = (y1 > y0 ? y1 - y0: y0 -y1);
  }

  // Expansion factor for y
  double yexpand = 1.0;
  if ( ran[1] > 1e-5 ) yexpand = ran[1]/ran[0];

  // Adjusting for device size
  if (dev.length() == 0)
    dev = NumericVector::create(2,1.0);

  yexpand = yexpand * (dev[0]/dev[1]);

  // Computing angle and adjusting for sign
  double alpha = angle(x0, y0/yexpand, x1, y1/yexpand);

  // Filling coords ------------------------------------------------------------
  coords.at(0,0) = x1;
  coords.at(0,1) = y1;

  // Left
  coords.at(1,0) = x1 - cos(alpha)*height + cos(beta+alpha)*width;
  coords.at(1,1) = y1 - (sin(alpha)*height - sin(beta+alpha)*width)*yexpand;

  // center
  coords.at(2,0) = x1 - cos(alpha)*height;
  coords.at(2,1) = y1 - sin(alpha)*height*yexpand;

  // Bottom
  coords.at(3,0) = x0;
  coords.at(3,1) = y0;

  // Back to the center
  coords.at(4,0) = coords.at(2,0);
  coords.at(4,1) = coords.at(2,1);

  // Right
  coords.at(5,0) = x1 - cos(alpha)*height + cos(-beta+alpha)*width;
  coords.at(5,1) = y1 - (sin(alpha)*height - sin(-beta+alpha)*width)*yexpand;

  return coords;
}

/** *R
# rm(list = ls())

X <- c(-9,-9)
ran <- X*1.1
h <- 1
w <- .5
beta <- pi/1.5

pol2 <- vector_polygon(0, 0, X[1], X[2], h, w, dev=par("pin"), beta = beta)

plot(pol2[,1], pol2[,2], xlim=ran, ylim=ran, col="white")
polygon(pol2[,1], pol2[,2], col=rgb(.5,.5,.9,.5))

text(pol2[,1], pol2[,2], text=1:3)
segments(0,0,X[1],X[1])

library(netdiffuseR)
data("medInnovationsDiffNet")
x <- plot_threshold(medInnovationsDiffNet)

*/


// [[Rcpp::export]]
List vertices_coords(
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & size,
    const arma::colvec & nsides,
    const arma::colvec & rot,
    NumericVector dev = NumericVector::create(),
    NumericVector ran = NumericVector::create()
) {


  // Checking sizes
  if (x.n_rows != y.n_rows) stop("-x- and -y- lengths do not coincide.");
  if (x.n_rows != size.n_rows) stop("-x- and -size- lengths do not coincide.");
  if (x.n_rows != nsides.n_rows) stop("-x- and -nsides- lengths do not coincide.");
  if (x.n_rows != rot.n_rows) stop("-x- and -rot- lengths do not coincide.");

  List out(x.n_rows);

  // If yexpand is too small, just throw an error
  if (ran.length() == 0) {
    ran = NumericVector::create(2);
    ran[0] = x.max() - x.min();
    ran[1] = y.max() - y.min();
  }

  // Expansion factor for y
  double yexpand = 1.0;
  if ( ran[1] > 1e-5 ) yexpand = ran[1]/ran[0];

  // Adjusting for device size
  if (dev.length() == 0)
    dev = NumericVector::create(2,1.0);

  yexpand = yexpand * (dev[0]/dev[1]);

  for (unsigned i=0;i<x.n_rows;i++) {
    // Getting inner degrees
    double alpha = PI - ((nsides(i) - 2.0)*PI)/nsides(i);
    double beta  = (PI - 2.0*PI/nsides(i))/2.0;

    // Getting step size
    double size_adj = 2.0*cos(beta)*size(i);

    // Suboutput and first coordinate
    arma::mat coords(nsides(i),2);
    coords(0,0) = x(i) + size(i)*cos(beta + PI + rot(i));
    coords(0,1) = y(i) + size(i)*sin(beta + PI + rot(i))*yexpand;

    double ALPHA = rot(i);
    for (int j=1; j<nsides(i); j++) {
      coords(j,0) = coords(j-1,0) + size_adj*cos(ALPHA);
      coords(j,1) = coords(j-1,1) + size_adj*sin(ALPHA)*yexpand;
      ALPHA += alpha;
    }

    // Assigning element
    out[i] = coords;
  }

  return out;
}

/* **R

# Parameters
n   <- c(3,4,5,6)
d   <- rep(.5, 4)
rot <- rep(-pi/6, 4)
y <- x <- c(1, 2, 3, 4)

plot.new()
plot.window(xlim=c(0,5), ylim=c(0,5))
axis(1)
axis(2)

# Computing coordinates
coords <- vertices_coords(x,y,d,n,rot, dev=par()$pin)

# polygon(coords)
invisible(lapply(coords, polygon))

invisible(sapply(x, function(x) symbols(x,x, .5, inches=FALSE, add=TRUE)))

*/
