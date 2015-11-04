#' Calculate the number of adoption changes between ego and alter.
#'
#' The changes are based on 16 categories combining (ego, alter) x (adopt in t) x
#' (adopt in t-1) (see details).
#'
#' @param graph Either an adjacency matrix or an array.
#' @param adopt nxT matrix. Cumulative addoption matrix.
#' @param period Integer. Optional to make the count for a particular period of time.
#' @param ... Further arguments to be passed to the method
#' @details The 16 categories are classified using the table that follows. The
#' first two Yes/No columns represent Ego's adoption of the innovation in t-1
#' and t; while the first two Yes/No rows represent Alter's adoption of the
#' innovation in t-1 and t respectively. So for example, number 4 means that
#' while neither of the two had addopted the innovation in t-1, both have in t.
#' At the same time, number 12 means that ego adopted the innovation in t, but
#' alter had already adopted in t-1 (so it has it in both, t and t-1).
#'
#' \tabular{rrrcccc}{
#'       \tab       \tab       \tab Alter \tab     \tab     \tab     \cr
#'       \tab       \tab t-1   \tab  No   \tab     \tab Yes \tab     \cr
#'       \tab t-1   \tab t     \tab  No   \tab Yes \tab No  \tab Yes \cr
#'   Ego \tab No    \tab No    \tab   1   \tab  2  \tab   9 \tab  10 \cr
#'       \tab       \tab Yes   \tab   3   \tab  4  \tab 11  \tab  12 \cr
#'       \tab Yes   \tab No    \tab   5   \tab  6  \tab 13  \tab  14 \cr
#'       \tab       \tab Yes   \tab   7   \tab  8  \tab 15  \tab  16
#' }
#'
#' @return An array of the count of selection changes for the 16 categories by node.
#' @references
#' Thomas W. Valente, Stephanie R. Dyal, Kar-Hai Chu, Heather Wipfli, Kayo
#' Fujimoto, Diffusion of innovations theory applied to global tobacco control
#' treaty ratification, Social Science & Medicine, Volume 145, November 2015,
#' Pages 89-97, ISSN 0277-9536
#' (\url{http://dx.doi.org/10.1016/j.socscimed.2015.10.001})
#' @export
select_egoalter <- function(graph, ...) UseMethod("select_egoalter")

#' @describeIn select_egoalter Method for arrays
select_egoalter.array <- function(graph, adopt, period=NULL) {

  # Computing selection mat and coersing into a single matrix
  nper <- dim(graph)[3]
  if (length(period) && !(period %in% 2:nper))
    stop('Invalid period selected. Should be between -', 2,'- and -', nper,'-')

  # Creating column names
  cn <- c('time', 'id',
          sprintf('select_a_%02d', 1:16), sprintf('select_d_%02d', 1:16),
          sprintf('select_s_%02d', 1:16))

  # Output parameters
  n   <- dim(graph)[1]
  ids <- 1:n

  if (length(period)) {
    out <- cbind(per=period, ids, do.call(cbind, select_egoalter_cpp(
      graph[,,period-1], graph[,,period],
      adopt[,period-1], adopt[,period])))

    # Assigning colnames
    colnames(out) <- cn

    return(out)
  }

  # Looping over periods
  out <- vector("list", nper-1)
  for (i in 2:nper) {
    out[[i-1]] <- cbind(time=i, ids, do.call(cbind, select_egoalter_cpp(
      graph[,,i-1], graph[,,i],adopt[,i-1], adopt[,i])))

    # Assigning colnames
    colnames(out[[i-1]]) <- cn
  }

  return(array(unlist(out), dim=c(n,length(cn),nper)))
}
