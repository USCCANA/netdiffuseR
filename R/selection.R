#' Selection function
#' @param adjmat Adjacency matrix
#' @param adopt Adoption matrix
#' @param period Integer time period
#' @export
select_egoalter <- function(...) UseMethod("select_egoalter")

#' @describeIn select_egoalter Method for array
select_egoalter.array <- function(adjmat, adopt, period=NULL) {

  # Computing selection mat and coersing into a single matrix
  nper <- dim(adjmat)[3]
  if (length(period) && !(period %in% 2:nper))
    stop('Invalid period selected. Should be between -', 2,'- and -', nper,'-')

  # Creating column names
  cn <- c('time', 'id',
          sprintf('select_a_%02d', 1:16), sprintf('select_d_%02d', 1:16),
          sprintf('select_s_%02d', 1:16))

  # Output parameters
  n   <- dim(adjmat)[1]
  ids <- 1:n

  if (length(period)) {
    out <- cbind(per=period, ids, do.call(cbind, select_egoalter_cpp(
      adjmat[,,period-1], adjmat[,,period],
      adopt[,period-1], adopt[,period])))

    # Assigning colnames
    colnames(out) <- cn

    return(out)
  }

  # Looping over periods
  out <- vector("list", nper-1)
  for (i in 2:nper) {
    out[[i-1]] <- cbind(time=i, ids, do.call(cbind, select_egoalter_cpp(
      adjmat[,,i-1], adjmat[,,i],adopt[,i-1], adopt[,i])))

    # Assigning colnames
    colnames(out[[i-1]]) <- cn
  }

  return(array(unlist(out), dim=c(n,length(cn),nper)))
}
