#' Calculate the number of adoption changes between ego and alter.
#'
#' This function calculates the 16 possible configurations between ego and alter
#' over two time points in terms of their behavior and tie changes.  From time
#' one to time two, given a binary state of behavior, ego and alter can be
#' related in 16 different ways. The function \code{adopt_changes} is just an
#' alias for \code{select_egoalter}.
#'
#' @templateVar dynamic TRUE
#' @template graph_template
#' @param adopt \eqn{n\times T}{n*T} matrix. Cumulative adoption matrix obtained from \code{\link{toa_mat}}.
#' @param period Integer scalar. Optional to make the count for a particular period of time.
#' @details
#' The 16 possibilities are summarized in this matrix:
#'
#' \tabular{rrrcccc}{
#'       \tab       \tab       \tab Alter \tab     \tab     \tab     \cr
#'       \tab       \tab \eqn{t-1}   \tab  No   \tab     \tab Yes \tab     \cr
#'       \tab \eqn{t-1}   \tab \eqn{t}     \tab  No   \tab Yes \tab No  \tab Yes \cr
#'   Ego \tab No    \tab No    \tab   1   \tab  2  \tab   9 \tab  10 \cr
#'       \tab       \tab Yes   \tab   3   \tab  4  \tab 11  \tab  12 \cr
#'       \tab Yes   \tab No    \tab   5   \tab  6  \tab 13  \tab  14 \cr
#'       \tab       \tab Yes   \tab   7   \tab  8  \tab 15  \tab  16
#' }
#'
#' The
#' first two Yes/No columns represent Ego's adoption of the innovation in \eqn{t-1}
#' and \eqn{t}; while the first two Yes/No rows represent Alter's adoption of the
#' innovation in \eqn{t-1} and t respectively. So for example, number 4 means that
#' while neither of the two had addopted the innovation in \eqn{t-1}, both have in \eqn{t}.
#' At the same time, number 12 means that ego adopted the innovation in \eqn{t}, but
#' alter had already adopted in \eqn{t-1} (so it has it in both, \eqn{t} and \eqn{t-1}).
#'
#' @return An object of class \code{diffnet_adoptChanges} and \code{data.frame}
#' with \eqn{n\times (T-1)}{n * (T-1)} rows and \eqn{2 + 16\times 3}{2 + 16 * 3}
#' columns. The column names are:
#' \item{\code{time}}{Integer represting the time period}
#' \item{\code{id}}{Node id}
#' \item{\code{select_a_01}, \dots, \code{select_a_16}}{Number of new links classified
#' between categories 1 to 16.}
#' \item{\code{select_d_01}, \dots, \code{select_d_16}}{Number of remove links classified
#' between categories 1 to 16.}
#' \item{\code{select_s_01}, \dots, \code{select_s_16}}{Number of unchanged links
#' classified between categories 1 to 16.}
#' @references
#' Thomas W. Valente, Stephanie R. Dyal, Kar-Hai Chu, Heather Wipfli, Kayo
#' Fujimoto, \emph{Diffusion of innovations theory applied to global tobacco control
#' treaty ratification}, Social Science & Medicine, Volume 145, November 2015,
#' Pages 89-97, ISSN 0277-9536
#' \doi{10.1016/j.socscimed.2015.10.001}
#' @export
#' @author George G. Vega Yon & Thomas W. Valente
#' @examples
#' # Simple example ------------------------------------------------------------
#' set.seed(1312)
#' dn <- rdiffnet(20, 5, seed.graph="small-world")
#'
#' ans <- adopt_changes(dn)
#' str(ans)
#' summary(ans)
select_egoalter <- function(graph, adopt, period=NULL) {

  if (missing(adopt))
    if (!inherits(graph, "diffnet"))
      stop("-adopt- should be provided when -graph- is not of class 'diffnet'")

  cls <- class(graph)

  if ("array" %in% cls) {
    select_egoalter.array(graph, adopt, period)
  } else if ("list" %in% cls) {
    select_egoalter.list(graph, adopt, period)
  } else if ("diffnet" %in% cls) {
    select_egoalter.list(graph$graph, graph$adopt, period)
  } else
    stopifnot_graph(graph)

}

#' @rdname select_egoalter
#' @export
adopt_changes <- select_egoalter

# @rdname select_egoalter
# @export
select_egoalter.array <- function(graph, adopt, period=NULL) {
  graph <- apply(graph, 3, methods::as, Class="dgCMatrix")
  select_egoalter.list(graph, adopt, period)
}

# @rdname select_egoalter
# @export
select_egoalter.list <- function(graph, adopt, period=NULL) {

  # Computing selection mat and coercing into a single matrix
  n <- dim(graph[[1]])[1]
  t <- length(graph)
  if (length(period) && !(period %in% 2:t))
    stop('Invalid period selected. Should be between -', 2,'- and -', t,'-')

  # Creating column, and row names
  cn <- c('time', 'id',
          sprintf('select_a_%02d', 1:16), sprintf('select_d_%02d', 1:16),
          sprintf('select_s_%02d', 1:16))

  # Output parameters
  ids <- rownames(graph[[1]])
  if (!length(ids)) ids <- 1:n

  tn <- names(graph)
  if (!length(tn)) tn <- 1:t

  if (length(period)) {
    out <- data.frame(time=tn[period], ids, select_egoalter_cpp(
      graph[[period-1]], graph[[period]],
      adopt[,period-1], adopt[,period]))

    # Assigning names
    colnames(out) <- cn

    return(out)
  }

  # Looping over periods
  out <- vector("list", t-1)
  for (i in 2:t) {
    out[[i-1]] <- data.frame(time=tn[i], ids, select_egoalter_cpp(
      graph[[i-1]], graph[[i]],adopt[,i-1], adopt[,i]))

    colnames(out[[i-1]]) <- cn
  }

  out <- do.call(rbind, out)
  class(out) <- c("diffnet_adoptChanges", class(out))
  return(out)
}

#' @export
#' @param object An object of class \code{diffnet_adoptChanges}.
#' @param ... Ignored.
#' @rdname select_egoalter
summary.diffnet_adoptChanges <- function(object, ...) {
  timevar <- object[["time"]]

  # Labels
  rn <- paste0(c(
    paste0("Ego|", rep("No (t-1)|", 2)),
    paste0("Ego|", rep("Yes (t-1)|", 2))
  ), c("No (t)","Yes (t)"))

  cn <- paste0(c(
    paste0("Alter|", rep("No (t-1)|", 2)),
    paste0("Alter|", rep("Yes (t-1)|", 2))
  ), c("No (t)","Yes (t)"))

  # Function for tabulating
  egoaltertab <- function(expr) {

    ans <- grepl(expr, colnames(object))
    ans <- object[,ans]
    colnames(ans) <- as.integer(
      gsub(".+_(?=[0-9]+)", "", colnames(ans), perl = TRUE))

    ans <- unclass(by(ans, timevar, colSums))

    ans <- lapply(ans, function(x) {
      x <- matrix(x, byrow = TRUE, ncol=2)
      x <- cbind(x[1:4,], x[5:8,])
      dimnames(x) <- list(rn, cn)
      x
    })

    names(ans) <- paste("period",names(ans))
    ans

  }

  # Added Links
  # debug(egoaltertab)
  links_added   <- egoaltertab("select_a")
  links_deleted <- egoaltertab("select_d")
  links_static  <- egoaltertab("select_s")

  list(
    `added (a)`     = links_added,
    `removed (d)`   = links_deleted,
    `unchanged (s)` = links_static
  )
}
