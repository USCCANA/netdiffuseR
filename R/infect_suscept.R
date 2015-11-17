#' Susceptibility and Infection
#'
#' Calculates infectiousness and susceptibility for each node in the graph
#'
#' @param graph Array of size \eqn{n\times n\times T}{n*n*T} as a dynamic graph
#' @param times Integer vector with times of adoption
#' @param normalize Logical. Whether or not to normalize the outcome
#' @param K Integer. Number of time periods to consider
#' @param r Double. Discount rate used when \code{expdiscount=TRUE}
#' @param expdiscount Logical. When TRUE, exponential discount rate is used (see details)
#' @details
#'
#' Normalization, \code{normalize=TRUE}, is applied by dividing the
#' resulting number from the infectiousness/susceptibility calcuation
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
#' Note that when \eqn{k=1}, the above formulas are equal to the ones presented
#' in the paper.
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
#'
#' @export
#' @return A numeric column vector of size \eqn{n} with either infection/susceptibility rates.
#'
infection <- function(graph, times, normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE) {
  UseMethod("infection")
}

#' @rdname infection
#' @export
infection.array <- function(graph, times, normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE) {
  infection_cpp(graph, times, normalize, K, r, expdiscount)
}

#' @rdname infection
#' @export
susceptibility <- function(graph, times, normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE) {
  UseMethod("susceptibility")
}

#' @rdname infection
#' @export
susceptibility.array <- function(graph, times, normalize=TRUE, K=1L, r=0.5, expdiscount=FALSE) {
  susceptibility_cpp(graph, times, normalize, K, r, expdiscount)
}

