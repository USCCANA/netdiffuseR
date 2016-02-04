#' Structure dependence test
#'
#' Test whether or not a network estimates can be considered as structure dependent, i.e.
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
#' \item{graph}{The graph passed to \code{boot_net}.}
#' \item{statistic}{The function \code{statistic} passed to \code{boot_net}.}
#' \item{boot}{A \code{boot} class object as return from the call to \code{boot}.}
#'
#' The output from the \code{hist} method is the same as \code{\link{hist.default}}.
#' @details
#' \code{boot_net} is a wrapper for the function \code{\link[boot:boot]{boot}} from the
#' \pkg{boot} package. As a difference, instead of resampling data the function performs
#' rewiring of the graph in each repetition using \code{\link{rewire_graph}} and applies
#' the function defined by the user in \code{statistic}.
#'
#' In \code{boot_net} \code{\dots} are passed to \code{boot}, otherwise are passed
#' to the corresponding method (\code{\link{hist}} for instance).
#'
#' From the \code{print} method, p-value for the null of the statistic been
#' equal between graph and its rewired versions is computed following Davidson
#' & MacKinnon
#'
#' \deqn{%
#' p(\tau)=2\times\min\left(\mbox{Pr}(t\leq\tau), \mbox{Pr}(t\geq\tau)\right) %
#' }{ %
#' p(tau) = 2*min[Pr(t<=tau), Pr(t>=tau)] %
#' }
#'
#' Where \eqn{\mbox{Pr}\{\cdot\}}{Pr(.)} is approximated using the
#' Empirical Distribution Function retrieved from the simulations.
#' @export
#' @references
#' On development.
#' @examples
#' # Creating a random graph
#' set.seed(881)
#' diffnet <- rdiffnet(100, 10)
#'
#' # Testing structure-dependency of threshold
#' res <- boot_net(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100)
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
#' res <- boot_net(diffnet, function(g) mean(threshold(g), na.rm=TRUE), R=100,
#' ncpus=4, parallel="multicore")
#' res
#' hist(res)
#' }
#' @author Vega Yon
boot_net <- function(
  graph,
  statistic,
  R,
  rewire.args=list(p=1, undirected=getOption("diffnet.undirected"), both.ends=TRUE),
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

  # Creating the object
  out <- list(
    graph = graph,
    statistic = statistic,
    boot=boot_res
  )

  return(structure(out, class="diffnet_boot"))
}

#' @export
#' @rdname boot_net
print.diffnet_boot <- function(x, ...) {
  with(x,  {
    tmean <- colMeans(boot$t, na.rm = TRUE)

    # Calc pval
    # To be conservative, in a two tail test we use the min of the two
    # So, following davidson & mckinnon Confidence intrval section,
    # p(tau) = 2 * min[F(tau), 1-F(tau)]
    test <- 2*min(mean(boot$t < boot$t0), mean(boot$t > boot$t0))
    cat("Network Rewiring graph (",nrow(boot$t)," simulations)\n",
        "# nodes           : ", x$graph$meta$n,"\n",
        "# of time periods : ", x$graph$meta$nper,"\n",
        paste(rep("-",80), collapse=""),"\n",
        " H0: t - t0 = 0\n",
        " t0      = ", boot$t0, "\n",
        " t       = ", tmean, "\n",
        " p-value = ", sprintf("%.5f",test), sep="")
  })
  invisible(x)
}

#' @export
#' @rdname boot_net
hist.diffnet_boot <- function(
  x,
  main="Distribution of Statistic on\nrewired network",
  xlab="Values of t",
  breaks=20,
  annotated=TRUE,
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
    mtext(expression(atop(plain("") %up% plain("")), t[0]),
          side = 1, at=x$boot$t0)
    mtext(expression(atop(plain("") %up% plain("")), t[]),
          side = 1, at=mt)
  }
  invisible(out)
}
# #' @rdname boot_net
# boot_thr <- boot_net(net, function(x) {
#   t <- threshold(x)
#   cbind(mean=mean(t, na.rm=TRUE), sd=sd(t, na.rm = TRUE))
# }
# , nsim)
#
# ttest<-(boot_thr$t[1,] - mean(threshold(net), na.rm=TRUE))/boot_thr$t[2,]
