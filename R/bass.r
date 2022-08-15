#' Bass Model
#'
#' Fits the Bass Diffusion model. In particular, fits an observed curve of
#' proportions of adopters to \eqn{F(t)}, the proportion of adopters at time
#' \eqn{t}, finding the corresponding coefficients \eqn{p}, Innovation rate,
#' and \eqn{q}, imitation rate.
#'
#' @param Time Integer vector with values greater than 0. The \eqn{t} parameter.
#' @param p Numeric scalar. Coefficient of innovation.
#' @param q Numeric scalar. Coefficient of imitation.
#' @param dat Either a diffnet object, or a numeric vector.
#' Observed cumulative proportion of adopters.
#' @param x An object of class \code{diffnet_bass}.
#' @param ... Further arguments passed to the method.
#'
#' @details The function fits the bass model with parameters \eqn{[p, q]} for
#' values \eqn{t = 1, 2, \dots, T}, in particular, it fits the following function:
#'
#' \deqn{
#'   F(t) = \frac{1 - \exp{-(p+q)t}}{1 + \frac{q}{p}\exp{-(p+q)t}}
#' }{
#'   F(t) = [1 - exp(-(p + q)*t)]/[1 + exp(-(p + q)*t)*(q/p)]
#' }
#'
#' Which is implemented in the \code{bass_F} function. The proportion of adopters
#' at time \eqn{t}, \eqn{f(t)} is:
#'
#' \deqn{
#'   f(t) = \left\{\begin{array}{ll}
#'   F(t), & t = 1 \\
#'   F(t) - F(t-1), & t > 1
#'   \end{array}\right.
#' }{
#'   f(t) = ifelse(t == 1, F(t), F(t) - F(t-1))
#' }
#'
#' and it's implemented in the \code{bass_f} function.
#'
#' For testing purposes only, the gradient of \eqn{F} with respect to \eqn{p}
#' and \eqn{q} is implemented in \code{bass_dF}.
#'
#' The estimation is done using \code{\link[stats:nls]{nls}}.
#'
#'
#' @return An object of class \code{nls} and \code{diffnet_bass}. For more
#' details, see \code{\link[stats:nls]{nls}} in the \pkg{stats} package.
#'
#' @examples
#' # Fitting the model for the Brazilian Farmers Data --------------------------
#' data(brfarmersDiffNet)
#' ans <- fitbass(brfarmersDiffNet)
#'
#' # All the methods that work for the -nls- object work here
#' ans
#' summary(ans)
#' coef(ans)
#' vcov(ans)
#'
#' # And the plot method returns both, fitted and observed curve
#' plot(ans)
#'
#' @references
#' Bass's Basement Institute Institute. The Bass Model. (2010).
#' Available at: \url{http://www.bassbasement.org/BassModel/Default.aspx}. (Accessed: 29th March 2017)
#' @name bass
#' @author George G. Vega Yon
#' @family statistics
NULL

#' @rdname bass
#' @export
fitbass <- function(dat, ...) UseMethod("fitbass")

#' @export
#' @rdname bass
fitbass.diffnet <- function(dat, ...) {
  .fitbass(cumulative_adopt_count(dat$cumadopt)["prop",], ...)
}

#' @export
#' @rdname bass
fitbass.default <- function(dat, ...) {
  .fitbass(as.vector(dat), ...)
}

.fitbass <- function(dat, ...) {

  # Constants
  Time <- seq_along(dat)

  # # Optimization (Fit nonlinear regression)
  ans <- stats::nls(dat ~ bass_F(Time, p, q), start=list(p = dat[1], q = .5 ),...)

  structure(c(
    ans,
    list(
      q      = stats::coef(ans)[2],
      p      = stats::coef(ans)[1],
      Fvals  = bass_F(Time, stats::coef(ans)[1], stats::coef(ans)[2]),
      fvals  = bass_f(Time, stats::coef(ans)[1], stats::coef(ans)[2]),
      nper   = length(Time),
      dat    = dat,
      Time   = Time
    )
  ), class = c("diffnet_bass", class(ans)))

}

#' @rdname bass
#' @param y Integer vector. Time (label).
#' @param pch Passed to \code{\link[graphics:plot]{matplot}}.
#' @param main Passed to \code{\link[graphics:plot]{matplot}}.
#' @param xlab Character scalar. Label of the \code{x} axis.
#' @param ylab Character scalar. Label of the \code{y} axis.
#' @param type Passed to \code{\link[graphics:plot]{matplot}}.
#' @param lty Passed to \code{\link[graphics:plot]{matplot}}.
#' @param col Passed to \code{\link[graphics:plot]{matplot}}.
#' @param bg Passed to \code{\link[graphics:plot]{matplot}}.
#' @param include.legend Logical scalar. When \code{TRUE}, draws a legend.
#' @param add Passed to \code{\link[graphics:plot]{matplot}}.
#' @export
plot.diffnet_bass <- function(
  x, y=1:length(x$m$lhs()), add=FALSE,
  pch  = c(21,24),
  main = "Bass Diffusion Model",
  ylab = "Proportion of adopters",
  xlab = "Time",
  type = c("b", "b"),
  lty  = c(2,1),
  col  = c("black","black"),
  bg   = c("lightblue","gray"),
  include.legend = TRUE,
  ...) {

  if (length(type) == 1) type <- rep(type,2)
  if (length(pch) == 1)  pch  <- rep(pch,2)
  if (length(lty) == 1)  lty  <- rep(lty,2)
  if (length(col) == 1)  lty  <- rep(col,2)
  if (length(bg) == 1)   bg   <- rep(bg,2)

  mat <- with(x, cbind(m$lhs(), m$fitted()))

  matplot(mat, xlab = xlab, ylab = ylab, main=main,
          ylim = c(0,1), type=type, lty=lty, bg=bg, col=col,
          pch=pch, add=add, ...)

  # Adding legend
  if (!add && include.legend)
    legend(
      "topleft", bty="n",
      legend = c("Observed Cumulative Adopters", "Predicted Cumulative Adopters"),
      col=col, pt.bg=bg, pch=pch
      )


  invisible(x)

}

#' @rdname bass
#' @export
bass_F <- function(Time, p, q) {
  (1 - exp(-(p + q)*Time))/(1 + (q/p)*exp(-(p+q)*Time))
}

#' @rdname bass
#' @export
bass_dF <- function(p, q, Time) {
  expv <- exp((p+q)*Time)

  rbind(
    expv * p^2*Time + q*(expv*(1 + p*Time) - 1)/
      (expv*p + q)^2,
    p * (1 + expv*(p*Time + q*Time - 1)) /
      (expv*p + q)^2
  )
}

#' @rdname bass
#' @export
bass_f <- function(Time, p, q) {
  ifelse(Time==1, bass_F(Time, p, q), bass_F(Time, p, q) - bass_F(Time-1, p, q))
}

