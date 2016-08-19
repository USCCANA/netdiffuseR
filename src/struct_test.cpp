// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "netdiffuser_extra.h"
using namespace Rcpp;

// [[Rcpp::export]]
double struct_test_mean(NumericVector & y,
                        std::string funname, bool self=false) {

  int    n =y.size();
  double m = (self? n*n/2:n*(n-1)/2);
  double ans = 0.0;

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  for (int i=0;i<n;i++)
    for (int j=i;j<n;j++) {
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

  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++)
      for (int k=0;k<n;k++) {
        if (self) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
        else if (i!=j && i!=k) ans+=fun(y[i],y[j])*fun(y[i],y[k])/m;
      }

      return ans - pow(struct_test_mean(y, funname, self),2.0);
}

typedef arma::sp_mat::const_iterator spiter;

// [[Rcpp::export]]
NumericVector hatf(const arma::sp_mat & G, const  NumericVector & Y,
                   std::string funname) {
  int n = G.n_rows;
  NumericVector ans(n);
  NumericVector d(n);

  // Fetching function
  funcPtr fun;
  st_getfun(funname, fun);

  // Preparing iterators
  spiter begin = G.begin();
  spiter end   = G.end();

  for (spiter i = begin; i!=end; i++) {
    for (spiter j = begin; j!=end; j++) {
      // I and j must be connected
      if (i.col() != j.row()) continue;

      // Are i and k connected
      // To see such we check at G in
      // G(i,k) > 0
      if (G.at(i.row(), j.col()) < 1e-15) continue;

      ans[i.row()] += fun(Y[j.row()], Y[j.col()]);
      d.at(i.row())++;
    }
  }

  for (int i=0;i<n;i++)
    if (d[i]>0) ans[i] /= (d[i] + 1e-15);
    else ans[i] = NA_REAL;

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

i <- 105
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
*/
