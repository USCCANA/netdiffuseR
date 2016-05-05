#' Susceptibility and Infection
#'
#' Calculates infectiousness and susceptibility for each node in the graph
#'
#' @param graph A dynamic graph (see \code{\link{netdiffuseR-graphs}}).
#' @param toa Integer vector with times of adoption (see details)
#' @param t0 Integer scalar. See \code{\link{toa_mat}}.
#' @param normalize Logical. Whether or not to normalize the outcome
#' @param K Integer scalar. Number of time periods to consider
#' @param r Numeric scalar. Discount rate used when \code{expdiscount=TRUE}
#' @param expdiscount Logical scalar. When TRUE, exponential discount rate is used (see details).
#' @param valued Logical scalar. When FALSE non-zero values in the adjmat are set to one.
#' @param outgoing Logical scalar. When \code{TRUE}, computed using outgoing ties.
#' @family statistics
#' @aliases susceptibility
#' @keywords univar
#' @seealso The user can visualize the distribution of both statistics
#' by using the function \code{\link{plot_infectsuscep}}
#' @details
#'
#' Normalization, \code{normalize=TRUE}, is applied by dividing the
#' resulting number from the infectiousness/susceptibility stat
#' by the number of individuals who adopted the innovation at
#' time \eqn{t}.
#'
#' Given that node \eqn{i} adopted the innovation in time \eqn{t}, its
#' Susceptibility is calculated as follows
#'
#' \deqn{S_i = \frac{%
#' \sum_{k=1}^K\sum_{j=1}^n x_{ij(t-k+1)}z_{j(t-k)}\times \frac{1}{w_k}}{%
#' \sum_{k=1}^K\sum_{j=1}^n x_{ij(t-k+1)}z_{j(1\leq t \leq t-k)} \times \frac{1}{w_k} }\qquad \mbox{for }i,j=1,\dots,n\quad i\neq j}{%
#' S(i) = [\sum_k \sum_j (x(ij,t-k+1) * z(j,t-k))/w(k)]/[\sum_k \sum_j (x(ij,t-k+1) * z(j, 1<=t<=t-k))/w(k)]  for  j != i}
#'
#' where \eqn{x_{ij(t-k+1)}}{x(ij,t-k+1)} is 1 whenever there's a link from \eqn{i}
#' to \eqn{j} at time \eqn{t-k+1}, \eqn{z_{j(t-k)}}{z(j,t-k)}
#' is 1 whenever individual \eqn{j} adopted the innovation at time \eqn{t-k},
#' \eqn{z_{j(1\leq t \leq t-k)}}{z(j, 1<=t<=t-k)} is 1 whenever
#' \eqn{j} had adopted the innovation up to \eqn{t-k}, and \eqn{w_k}{w(k)} is
#' the discount rate used (see below).
#'
#' Similarly, infectiousness is calculated as follows
#'
#' \deqn{I_i = \frac{%
#' \sum_{k=1}^K \sum_{j=1}^n x_{ji(t+k-1)}z_{j(t+k)}\times \frac{1}{w_k}}{%
#' \sum_{k=1}^K \sum_{j=1}^n x_{ji(t+k-1)}z_{j(t+k\leq t \leq T)}\times \frac{1}{w_k} }\qquad \mbox{for }i,j=1,\dots,n\quad i\neq j}{%
#' I(i) = [\sum_k \sum_j (x(ji,t) * z(j,t+1))/w(k)]/[\sum_k \sum_j (x(ji,t) * z(j, t+1<=t<=T))/w(k)]   for  j != i}
#'
#' It is worth noticing that, as we can see in the formulas, while susceptibility
#' is from alter to ego, infection is from ego to alter.
#'
#' When \code{outgoing=FALSE} the algorithms are based on incoming edges, this is
#' the adjacency matrices are transposed swapping the indexes \eqn{(i,j)} by
#' \eqn{(j,i)}. This can be useful for some users.
#'
#' Finally, by default both are normalized by the number of individuals who
#' adopted the innovation in time \eqn{t-k}. Thus, the resulting formulas,
#' when \code{normalize=TRUE}, can be rewritten as
#'
#' \deqn{%
#' S_i' = \frac{S_i}{\sum_{k=1}^K\sum_{j=1}^nz_{j(t-k)}\times \frac{1}{w_k}} %
#' \qquad I_i' = \frac{I_i}{\sum_{k=1}^K\sum_{j=1}^nz_{j(t-k)} \times\frac{1}{w_k}}}{%
#' S(i)' = S(i)/[\sum_k \sum_j z(j,t-k)/w(k)]
#'
#' I(i)' = I(i)/[\sum_k \sum_j z(j,t-k)/w(k)]}
#'
#' For more details on these measurements, please refer to the vignette titled
#' \emph{Time Discounted Infection and Susceptibility}.
#'
#' @section Discount rate:
#'
#' Discount rate, \eqn{w_k}{w(k)} in the formulas above, can be either exponential
#' or linear. When \code{expdiscount=TRUE}, \eqn{w_k = (1 + r)^{k-1}}{w(k) = (1+r)^(k-1)}, otherwise
#' it will be \eqn{w_k = k}{w(k)=k}.
#'
#' Note that when \eqn{K=1}, the above formulas are equal to the ones presented
#' in Valente et al. (2015).
#'
#' @references
#' Thomas W. Valente, Stephanie R. Dyal, Kar-Hai Chu, Heather Wipfli, Kayo
#' Fujimoto, Diffusion of innovations theory applied to global tobacco control
#' treaty ratification, Social Science & Medicine, Volume 145, November 2015,
#' Pages 89-97, ISSN 0277-9536
#' (\url{http://dx.doi.org/10.1016/j.socscimed.2015.10.001})
#'
#' Myers, D. J. (2000). "The Diffusion of Collective Violence: Infectiousness,
#' Susceptibility, and Mass Media Networks". American Journal of Sociology, 106(1),
#' 173–208. doi:10.1086/303110
#' @examples
#'
#' # Creating a random dynamic graph
#' set.seed(943)
#' graph <- rgraph_er(n=100, t=10)
#' toa <- sample.int(10, 100, TRUE)
#'
#' # Computing infection and susceptibility (K=1)
#' infection(graph, toa)
#' susceptibility(graph, toa)
#'
#' # Now with K=4
#' infection(graph, toa, K=4)
#' susceptibility(graph, toa, K=4)
#'
#' @export
#' @return A numeric column vector (matrix) of size \eqn{n} with either infection/susceptibility rates.
#' @author George G. Vega Yon
infection <- function(graph, toa, t0=NULL,
                      normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE,
                      valued   = getOption("diffnet.valued", FALSE),
                      outgoing = getOption("diffnet.outgoing", TRUE)) {

  # Checking the times argument
  if (missing(toa))
    if (!inherits(graph, "diffnet")) {
      stop("-toa- should be provided when -graph- is not of class 'diffnet'")
    } else {
      toa <- graph$toa
      t0    <- min(graph$meta$pers)
    }

  # Checking baseline time
  if (!length(t0)) t0 <- min(toa, na.rm=TRUE)

  switch (class(graph),
    array = infection.array(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    list = infection.list(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    diffnet = infection.list(graph$graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    stopifnot_graph(graph)
  )
}

# @rdname infection
# @export
infection.array <- function(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing) {
  toa <- toa - t0 + 1L
  t <- dim(graph)[3]
  n <- nrow(graph)
  ngraph <- vector("list", t)

  for(i in 1:t)
    ngraph[[i]] <- methods::as(graph[,,i], "dgCMatrix")

  out <- infection_cpp(ngraph, toa, normalize, K, r, expdiscount, n, valued, outgoing)

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:n
  dimnames(out) <- list(rn, "infection")

  out
}

# @rdname infection
# @export
infection.list <- function(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  toa <- toa - t0 + 1L

  out <- infection_cpp(graph, toa, normalize, K, r, expdiscount, n, valued, outgoing)

  # Naming
  rn <- rownames(graph[[1]])
  if (!length(rn)) rn <- 1:n
  dimnames(out) <- list(rn, "infection")

  out
}

#' @rdname infection
#' @export
susceptibility <- function(graph, toa, t0=NULL, normalize=TRUE, K=1L, r=0.5,
                           expdiscount=FALSE,
                           valued=getOption("diffnet.valued",FALSE),
                           outgoing=getOption("diffnet.outgoing",TRUE)) {
  # Checking the toa argument
  if (missing(toa))
    if (!inherits(graph, "diffnet")) {
      stop("-toa- should be provided when -graph- is not of class 'diffnet'")
    } else {
      toa <- graph$toa
      t0    <- min(graph$meta$pers)
    }

  # Checking baseline time
  if (!length(t0)) t0 <- min(toa, na.rm=TRUE)

  switch (class(graph),
    array = susceptibility.array(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    list = susceptibility.list(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    diffnet = susceptibility.list(graph$graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing),
    stopifnot_graph(graph)
  )
}

# @rdname infection
# @export
susceptibility.list <- function(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing) {
  t <- length(graph)
  n <- nrow(graph[[1]])
  toa <- toa - t0 + 1L

  out <- susceptibility_cpp(graph, toa, normalize, K, r, expdiscount, n, valued, outgoing)

  # Naming
  rn <- rownames(graph[[1]])
  if (!length(rn)) rn <- 1:n
  dimnames(out) <- list(rn, "susceptibility")

  out
}

# @rdname infection
# @export
susceptibility.array <- function(graph, toa, t0, normalize, K, r, expdiscount, valued, outgoing) {
  toa <- toa - t0 + 1L
  t <- dim(graph)[3]
  n <- nrow(graph)
  ngraph <- vector("list", t)

  for(i in 1:t)
    ngraph[[i]] <- methods::as(graph[,,i], "dgCMatrix")

  out <- susceptibility_cpp(ngraph, toa, normalize, K, r, expdiscount, n, valued, outgoing)

  # Naming
  rn <- rownames(graph)
  if (!length(rn)) rn <- 1:n
  dimnames(out) <- list(rn, "susceptibility")

  out
}
