#' Structure dependence test
#'
#' Test whether or not a network estimates can be considered structurally dependent, i.e.
#' a function of the network structure. By rewiring the graph and calculating
#' a particular statistic \eqn{t}, the test compares the observed mean of \eqn{t}
#' against the empirical distribution of it obtained from rewiring the network.
#'
#' @param graph A \code{\link{diffnet}} graph.
#' @param statistic A function that returns either a scalar or a vector.
#' @param R Integer scalar. Number of repetitions.
#' @param rewire.args List. Arguments to be passed to \code{\link{rewire_graph}}
#' @param x A \code{diffnet_boot} class object.
#' @param ... Further arguments passed to the method (see details).
#' @param main Character scalar. Title of the histogram.
#' @param xlab Character scalar. x-axis label.
#' @param breaks Passed to \code{\link{hist}}.
#' @param annotated Logical scalar. When TRUE marks the observed data average and the simulated data average.
#' @return A list of class \code{diffnet_bot} containing the following:
#' \item{graph}{The graph passed to \code{struct_test}.}
#' \item{p.value}{The resulting p-value of the test (see details).}
#' \item{t0}{The observed value of the statistic.}
#' \item{mean_t}{The average value of the statistic applied to the simulated networks.}
#' \item{R}{Number of simulations.}
#' \item{statistic}{The function \code{statistic} passed to \code{struct_test}.}
#' \item{boot}{A \code{boot} class object as return from the call to \code{boot}.}
#'
#' The output from the \code{hist} method is the same as \code{\link{hist.default}}.
#' @details
#' \code{struct_test} is a wrapper for the function \code{\link[boot:boot]{boot}} from the
#' \pkg{boot} package. Instead of resampling data--vertices or edges--in each iteration the function
#' rewires the original graph using \code{\link{rewire_graph}} and applies
#' the function defined by the user in \code{statistic}. In particular, the \code{"swap"} algorithm
#' is used in order to preserve the degree sequence of the graph, in other words,
#' each rewired version of the original graph has the same degree sequence.
#'
#' In \code{struct_test} \code{\dots} are passed to \code{boot}, otherwise are passed
#' to the corresponding method (\code{\link{hist}} for instance).
#'
#' From the \code{print} method, p-value for the null of the statistic been
#' equal between graph and its rewired versions is computed as follows
#'
#' \deqn{%
#' p(\tau)=2\times\min\left(\mbox{Pr}(t\leq\tau), \mbox{Pr}(t\geq\tau)\right) %
#' }{ %
#' p(tau) = 2*min[Pr(t<=tau), Pr(t>=tau)] %
#' }
#'
#' Where \eqn{\mbox{Pr}\{\cdot\}}{Pr(.)} is approximated using the
#' Empirical Distribution Function retrieved from the simulations.
#'
#' The test is actually on development by Vega Yon and Valente. A copy of the
#' working paper can be distributed upon request to \email{g.vegayon@gmail.com}
#' @export
#' @references
#' Vega Yon, George G. and Valente, Thomas W. (On development).
#'
#' Davidson, R., & MacKinnon, J. G. (2004). Econometric Theory and Methods. New York:
#' Oxford University Press.
#' @examples
#' # Creating a random graph
#' set.seed(881)
#' diffnet <- rdiffnet(500, 10, seed.graph="small-world")
#'
#' # Testing structure-dependency of threshold
#' res <- struct_test(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100)
#' res
#' hist(res)
#'
#' # Adding a legend
#' legend("topright", bty="n",
#'  legend=c(
#'    expression(t[0]:~Baseline),
#'    expression(t:~Rewired~average)
#'  )
#'  )
#'
#' # Running in parallel fashion
#' \dontrun{
#' res <- struct_test(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100,
#' ncpus=4, parallel="multicore")
#' res
#' hist(res)
#' }
#' @author George G. Vega Yon
struct_test <- function(
  graph,
  statistic,
  R,
  rewire.args=list(
    p          = c(2000, rep(100, nslices(graph) - 1)),
    undirected = getOption("diffnet.undirected", FALSE),
    copy.first = TRUE,
    algorithm  = "swap"
    ),
  ...
  ) {

  # Checking class
  if (!inherits(graph, "diffnet"))
    stop("-graph- must be of class diffnet.")

  # Preparing the call to boot
  rewire.args$graph <- graph
  statisticpll <- function(d, i, fn, rewire.args, ...) {
    fn(do.call(rewire_graph, rewire.args))
  }

  # Calling boot
  boot_res <- boot::boot(1, statisticpll, R=R, fn=statistic, rewire.args=rewire.args,
                   ...)

  # The t0 must be applied with no rewiring!
  boot_res$t0 <- statistic(graph)

  # Calc pval
  # To be conservative, in a two tail test we use the min of the two
  # So, following davidson & mckinnon Confidence intrval section,
  # p(tau) = 2 * min[F(tau), 1-F(tau)]
  p.value <- 2*min(mean(boot_res$t < boot_res$t0), mean(boot_res$t > boot_res$t0))

  # Creating the object
  out <- list(
    graph       = graph,
    p.value     = p.value,
    t0          = boot_res$t0,
    mean_t      = colMeans(boot_res$t, na.rm = TRUE),
    R           = R,
    statistic   = statistic,
    boot        = boot_res
  )

  return(structure(out, class="diffnet_struct_test"))
}

#' @export
#' @rdname struct_test
print.diffnet_struct_test <- function(x, ...) {
  with(x,  {
    tmean <- colMeans(boot$t, na.rm = TRUE)

    cat("Structure dependence test\n",
        "# Simulations     : ", formatC(nrow(boot$t), digits = 0, format = "f", big.mark = ","),"\n",
        "# nodes           : ", formatC(x$graph$meta$n, digits = 0, format = "f", big.mark = ","),"\n",
        "# of time periods : ", formatC(x$graph$meta$nper, digits = 0, format = "f", big.mark = ","),"\n",
        paste(rep("-",80), collapse=""),"\n",
        " H0: t - t0 = 0 (no structure dependency)\n",
        "   t0 (observed) = ", t0, "\n",
        "   t (simulated) = ", mean_t, "\n",
        "   p-value = ", sprintf("%.5f", p.value), sep="")
  })
  invisible(x)
}

#' @export
#' @param b0 Character scalar. When \code{annotated=TRUE}, label for the value of \code{b0}.
#' @param b Character scalar. When \code{annotated=TRUE}, label for the value of \code{b}.
#' @rdname struct_test
hist.diffnet_struct_test <- function(
  x,
  main="Empirical Distribution of Statistic",
  xlab=expression(Values~of~t),
  breaks=20,
  annotated=TRUE,
  b0=expression(atop(plain("") %up% plain("")), t[0]),
  b =expression(atop(plain("") %up% plain("")), t[]),
  ...) {

  out <- hist(x$boot$t,  breaks=breaks, plot=FALSE,...)

  ran <- range(out$mids)
  if (annotated) {
    mt <- mean(x$boot$t, na.rm=TRUE)
    ran <- range(c(ran, mt, x$boot$t0))
  }

  plot(out, main=main, xlab=xlab, xlim = ran)

  # Adding margin note
  if (annotated) {
    mtext(b0, side = 1, at=x$boot$t0)
    mtext(b, side = 1, at=mt)
  }
  invisible(out)
}
# #' @rdname struct_test
# boot_thr <- struct_test(net, function(x) {
#   t <- threshold(x)
#   cbind(mean=mean(t, na.rm=TRUE), sd=sd(t, na.rm = TRUE))
# }
# , nsim)
#
# ttest<-(boot_thr$t[1,] - mean(threshold(net), na.rm=TRUE))/boot_thr$t[2,]
