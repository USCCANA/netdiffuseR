// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
double struct_test_mean(NumericVector & y,
                        std::string funname, bool self=false) {

  int    n =y.size();
  double m = (self? n*n:n*(n-1));
  double ans = 0.0;

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) {
      if (self) ans+= fun(y[i], y[j])/m;
      else if (i!=j) ans+= fun(y[i], y[j])/m;
    }

    return ans;
}


// [[Rcpp::export]]
double struct_test_var(NumericVector & y, std::string funname, bool self=false)  {

  int n      = y.size();
  double m   = (self ? n*n*n : n*(n-1)*(n-1));
  double ans = 0.0;

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j)
      for (int k=0;k<n;++k) {
        if (self) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
        else if (i!=j && i!=k) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
      }

      return ans - pow(struct_test_mean(y, funname, self),2.0);
}

typedef arma::sp_mat::const_iterator spiter;

//' Computes variance of \eqn{Y} at ego level
//' @param graph A matrix of size \eqn{n\times n}{n*n} of class \code{dgCMatrix}.
//' @param Y A numeric vector of length \eqn{n}.
//' @param funname Character scalar. Comparison to make (see \code{\link{vertex_covariate_compare}}).
//' @param all Logical scalar. When \code{FALSE} (default) \eqn{f_i} is mean at
//' ego level. Otherwise is fix for all i (see details).
//' @details
//'
//' For each vertex \eqn{i} the variance is computed as follows
//'
//' \deqn{%
//' (\sum_j a_{ij})^{-1}\sum_j a_{ij} \left[f(y_i,y_j) - f_i\right]^2
//' }{%
//' (sum_j a(ij))^(-1) * \sum_j a(ij) * [f(y(i),y(j)) - f(i)]^2
//' }
//'
//' Where \eqn{a_{ij}}{a(ij)} is the ij-th element of \code{graph}, \eqn{f} is
//' the function specified in \code{funname}, and, if \code{all=FALSE}
//' \eqn{f_i = \sum_j a_{ij}f(y_i,y_j)^2/\sum_ja_{ij}}{f(i)=\sum_j a(ij)f(y(i), y(j))^2/\sum_j a(ij)},
//' otherwise \eqn{f_i = f_j = \frac{1}{n^2}\sum_{i,j}f(y_i,y_j)}{f(i)=f(j)=(1/n^2)\sum_(i,j) f(y_i,y_j)}
//'
//'
//' This is an auxiliary function for \code{\link{struct_test}}. The idea is
//' to compute an adjusted measure of disimilarity between vertices, so the
//' closest in terms of \eqn{f} is \eqn{i} to its neighbors, the smaller the
//' relative variance.
//' @return A numeric vector of length \eqn{n}.
//' @export
//' @seealso \code{\link{struct_test}}
//' @family statistics
// [[Rcpp::export]]
NumericVector ego_variance(const arma::sp_mat & graph, const NumericVector & Y,
                       std::string funname, bool all=false) {

  // Initialization
  int n = Y.length();
  NumericVector ans(n);
  NumericVector fhat(n);
  arma::sp_mat degree = sum(graph,1);

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  // Preparing iterators
  spiter begin = graph.begin();
  spiter end   = graph.end();

  // Which value to use as mean
  if (!all) {
    for (spiter i=begin; i!=end;++i)
      fhat[i.row()] += graph.at(i.row(),i.col())*fun(Y[i.row()], Y[i.col()])/
      (degree[i.row()] + 1e-15);
  } else {
    double val = 0.0;
    for (int i=0;i<n;++i)
      for (int j=0;j<n;++j)
        val += fun(Y[i],Y[j])/n/n;

    fhat.fill(val);
  }

  // Preparing iterators
  begin = graph.begin();
  end   = graph.end();

  // Iterating
  for (spiter i = begin; i!=end; ++i) {
    if (NumericVector::is_na(fhat[i.row()])) {
      ans[i.row()] = NA_REAL;
      continue;
    }
    ans[i.row()] +=  graph.at(i.row(),i.col())*powf(fun(Y[i.row()], Y[i.col()]) - fhat[i.row()],2.0) /
      (degree[i.row()] + 1e-15);
  }

  return ans;

}



/***R
set.seed(1231)
Y <- rnorm(100)

# Greater
struct_test_mean(Y, "greater")
struct_test_var(Y, "greater")

# Diff
struct_test_mean(Y, "distance")
struct_test_var(Y, "distance")

# Equal
struct_test_mean(round(Y,1), "equal")
struct_test_var(round(Y,1), "equal")
#
# L <- letters[1:3]
# for(i in 1:3)
#   for(j in i:3)
#     for (k in i:3) {
#
#     }

library(netdiffuseR)
x     <- kfamilyDiffNet
index <- which(x$toa != max(x$toa, na.rm = TRUE))
x     <- x[index,]
Mobs <- mean(threshold(x), na.rm = TRUE)
# Mobs <- sum(threshold(x)*dgr(x, cmode = "outdegree")[,1]/sum(dgr(x, cmode = "outdegree")[,1]), na.rm = TRUE)
M <- struct_test_mean(x$toa, "smaller")
S <- struct_test_var(x$toa, "smaller")

pval <- pnorm(sqrt(nnodes(x))*(M-Mobs)/sqrt(S))
ifelse(pval >.5,1-pval,pval)*2

d <- dgr(x)[,1]
S <- sapply(sample(1:nnodes(x), 1e4, TRUE), function(i) {
  mean(sample(x$toa, d[i]) < x$toa[i])
})

mean(S, na.rm = TRUE)
var(S, na.rm = TRUE)
pval <- pnorm(sqrt(nnodes(x))*(Mobs-M)/sd(S))
pval <- ifelse(pval > .5, 1-pval, pval)*2
pval
# struct_test(x, function(x) mean(threshold(x), na.rm = TRUE), R=500)

x <- matrix(rnorm(5e5*3), ncol=3)

mean((x[,1] > x[,2])*(x[,1] > x[,3]))
mean(x[,1]>x[,2])

# Med innovations
x     <- kfamilyDiffNet
index <- which(x$toa != max(x$toa, na.rm = TRUE))
x     <- x[index,]

G   <- x$graph[[1]]
ans <- hatf(G, x$toa, "distance")

i <- 73
g   <- which(G[i,] != 0)

toas <- NULL
g
for (j in g) {
  whichg <- which(G[j,] !=0)
  print(whichg)
  for (k in whichg)
    if (G[i,k]!=0) {
      message(sprintf("%03d--%03d", j, k))
      toas <- c(toas, abs(x$toa[j] - x$toa[k]))
    }
}
mean(toas,na.rm=TRUE);ans[i]
stop()

# Another example
set.seed(331)
G <- rdiffnet(1e3, t=10, rewire = FALSE, seed.graph = "small-world",
              rgraph.args = list(k=14,p=.2), threshold.dist = function(x) .2)
ans0 <- hatf(G$graph[[1]],G$toa,"<=")
ans1 <- hatf(G$graph[[1]],sample(G$toa, 1e3, TRUE),"<=")
struct_test_var(G$toa, "distance")
plot(ans0,ans1)
abline(0,1)
cor(ans0,ans1, use="complete.obs")

t.test(ans0,ans1)

x$toa[c(1,6)]

diffnet <- medInnovationsDiffNet
G <- diffnet$graph[[1]]
Y <- diffnet$toa

ans <- struct_test(G, function(g) mean(ego_variance(g,Y,"quaddist",TRUE)), R=1000)
ans
*/
