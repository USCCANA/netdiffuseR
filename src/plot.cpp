// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec seq_cpp(double from, double to, int lengthout) {
  arma::vec out(lengthout);
  double step = (to - from)/(lengthout - 1);
  for(int i=0;i<lengthout;i++)
    out(i) = from + i*step;

  return out;

}

//' Distribution over a grid
//'
//' Distribution of pairs over a grid of fix size (number of elements)
//'
//' @param x Numeric vector of size \eqn{n}
//' @param y Numeric vector of size \eqn{n}
//' @param nlevels Integer scalar. Number of bins to return
//' @details
//'
//' This function ment for internal use only.
//'
//' @export
//' @keywords misc
//' @seealso Used by \code{\link{plot_infectsuscep}}
//' @return Returns a list with three elements
//' \item{x}{Numeric vector of size \code{nlevels} with the class marks for x}
//' \item{y}{Numeric vector of size \code{nlevels} with the class marks for y}
//' \item{z}{Numeric matrix of size \code{nlevels} by \code{nlevels} with the distribution %
//' of the elements in terms of frecuency}
//' @section Examples:
//' \code{# Generating random vectors of size 100}
//'
//' \code{x <- rnorm(100)}
//'
//' \code{y <- rnorm(100)}
//'
//' \code{# Calculating distribution}
//'
//' \code{grid_distribution(x,y,20)}
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

/***R
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

// [[Rcpp::export]]
List edges_coords(
    const arma::mat & graph,
    const arma::colvec & toa,
    const arma::colvec & x,
    const arma::colvec & y,
    const arma::colvec & vertex_cex,
    bool undirected=true,
    bool no_contemporary=true) {

  int n = graph.n_cols;

  // The output matrix has the following
  // - x0 and y0
  // - x1 and y1
  // - alpha
  // - dist
  // - size0 and size1
  // - mutual
  std::vector< double > x0;
  std::vector< double > y0;
  std::vector< double > x1;
  std::vector< double > y1;
  std::vector< double > alpha;
  std::vector< double > dist;
  std::vector< double > size0;
  std::vector< double > size1;
  // std::vector< int > mutual;

  // Rescaling the vertex sizes
  arma::colvec vertex_size(vertex_cex);
  double mins = vertex_size.min();
  double maxs = vertex_size.max();

//   for(int i=0;i<n;i++) {
//     vertex_size(i) = (vertex_size(i) - mins + 1e-5)/(maxs-mins+ 1e-5)/4.0;
//   }

  for(int i=0;i<n;i++) {
    // Verifying undirected or not
    int m=n;
    if (undirected) m=i;

    for(int j=0;j<m;j++) {
      if (!graph(i,j)) continue;

      if (no_contemporary && (toa(i)==toa(j)) ) continue;

      // Euclidean distance
      double d = pow( pow(x(i) - x(j), 2.0) + pow(y(i) - y(j), 2.0) , 0.5 );
      dist.push_back(d);

      // Computing the elevation degree
      double a = acos( (x[i] - x[j])/d );
      alpha.push_back(a);

      // Adding the xs and the ys
      x0.push_back(x(i));
      x1.push_back(x(j) + cos(a)*vertex_size(j));

      y0.push_back(y(i));

      // The formula needs an extra help to figure out the ys
      if (y(i) < y(j)) y1.push_back(y(j) - sin(a)*vertex_size(j));
      else             y1.push_back(y(j) + sin(a)*vertex_size(j));

      // Now the sizes
      size0.push_back(vertex_size(i));
      size1.push_back(vertex_size(j));
    }
  }

  // Building up the output
  int e = x0.size();
  NumericMatrix out(e,8);
  for(int i=0; i<e; i++) {
    out(i,0) = x0[i];
    out(i,1) = y0[i];
    out(i,2) = x1[i];
    out(i,3) = y1[i];
    out(i,4) = size0[i];
    out(i,5) = size1[i];
    out(i,6) = alpha[i];
    out(i,7) = dist[i];
  }

  colnames(out) = CharacterVector::create("x0", "y0", "x1", "y1", "size0",
           "size1", "alpha", "dist");

  return List::create(_["sizes"]=vertex_size, _["edges"]=out);

}

/***R
set.seed(123)
graph <- rand_graph()
toa <- sample(1:5, 10, TRUE)
pos <- sna::gplot.layout.random(graph, NULL)
cex <- seq(1,5,length.out = 10)

arr <- edges_coords(graph, toa, pos[,1], pos[,2], cex/20)

plot(pos, col="white", xlim= c(-2,2), ylim= c(-2,2))
with(arr, arrows(edges[,1], edges[,2], edges[,3], edges[,4]))
symbols(pos, circles=arr$sizes, add=TRUE, inches = FALSE, bg="lightblue")
text(pos[,1], pos[,2], labels = 1:10)
*/
