#' Convert survey-like data to a diffnet object
#'
#' This convenient function turns network nomination data into diffnet obejects.
#' It works as a wrapper of \code{\link{edgelist_to_adjmat}}
#' and \code{\link{as_diffnet}}, and as a preprocessing routine of survey-like data,
#' creating the respective edgelist to be imported as a network.
#'
#' @inheritParams edgelist_to_adjmat
#' @param dat A data frame.
#' @param idvar Character scalar. Name of the id variable.
#' @param netvars Character vector. Names of the network nomination variables.
#' @param toavar Character scalar. Name of the time of adoption variable.
#' @param groupvar Character scalar. Name of cohort variable (e.g. city).
#' @param no.unsurveyed Logical scalar. When TRUE the nominated individuals
#' that do not show in \code{idvar} are set to \code{NA} (see details).
#' @param ... Further arguments to be passed to \code{\link{as_diffnet}}
#' @details
#'
#' All of \code{idvar}, \code{netvars}, \code{toavar} and \code{groupvar}
#' must be integers. Were these numeric they are coerced into integers, otherwise,
#' when neither of both, the function returns with error.
#'
#' In field work it is not unusual that some respondents nominate unsurveyed
#' individuals. In such case, in order to exclude them from the analysis,
#' the user can set \code{no.unsurveyed=TRUE} (the default), telling the
#' function to exclude such individuals from the adjacency matrix. This is
#' done by setting variables in \code{netvars} equal to \code{NA} when the
#' nominated id can't be found in \code{idvar}.
#'
#' If the network nomination process was done in different groups (location
#' for example) the survey id numbers may be define uniquely within each group
#' but not across groups (there may be many individuals with \code{id=1},
#' for example). To encompass this issue, the user can tell the function what
#' variable can be used to distinguish between groups through the \code{groupvar}
#' argument. When \code{groupvar} is provided, function redifines \code{idvar}
#' and the variables in \code{netvars} as follows:
#'
#' \preformatted{
#'    dat[[idvar]] <- dat[[idvar]] + dat[[groupvar]]*z
#' }
#'
#' Where \code{z = 10^nchar(max(dat[[idvar]]))}.
#'
#' @export
#' @return A \code{\link{diffnet}} object.
#' @seealso \code{\link{fakesurvey}}
#' @author
#' Vega Yon
#' @examples
#' # Loading a fake survey (data frame)
#' data(fakesurvey)
#'
#' # Diffnet object keeping isolated vertices ----------------------------------
#' dn1 <- survey_to_diffnet(fakesurvey, "id", c("net1", "net2", "net3"), "toa",
#'    "group", keep.isolates=TRUE)
#'
#' # Diffnet object NOT keeping isolated vertices
#' dn2 <- survey_to_diffnet(fakesurvey, "id", c("net1", "net2", "net3"), "toa",
#'    "group", keep.isolates=FALSE)
#'
#' # dn1 has an extra vertex than dn2
#' dn1
#' dn2
#'
#' # Reproducing medInnovationsDiffNet object ----------------------------------
#' data(medInnovations)
#'
#' # What are the netvars
#' netvars <- names(medInnovations)[grepl("^net", names(medInnovations))]
#'
#' medInnovationsDiffNet2 <- survey_to_diffnet(
#'    medInnovations,
#'    "id", netvars, "toa", "city")
#'
#' medInnovationsDiffNet2
#'
#' # Comparing with the package's version
#' all(diffnet.toa(medInnovationsDiffNet2) == diffnet.toa(medInnovationsDiffNet)) #TRUE
#' all(
#'    diffnet.attrs(medInnovationsDiffNet2, as.df = TRUE) ==
#'    diffnet.attrs(medInnovationsDiffNet, as.df = TRUE),
#'    na.rm=TRUE) #TRUE
#'
#'
survey_to_diffnet <- function(
  dat, idvar, netvars, toavar,
  groupvar=NULL,
  no.unsurveyed=TRUE,
  t = NULL,
  undirected = getOption("diffnet.undirected", FALSE),
  self = getOption("diffnet.self", FALSE),
  multiple=getOption("diffnet.multiple", FALSE),
  keep.isolates=TRUE, recode.ids=TRUE,
  ...) {

  # Coercing data into numeric variables
  for (x in c(idvar, groupvar, netvars, toavar)) {
    cl <- class(dat[[x]])
    if (cl == "integer") next
    else if (cl == "numeric") {
      warning("Coercing -",x,"- into integer.")
      dat[[x]] <- as.integer(dat[[x]])
    } else {
      stop("The variable ",x, "is not integer. All variables listed on ",
           "-idvar-, -toavar-, -groupvar- and -netvars- should be integer.")
    }
  }

  # Changing ids
  if (length(groupvar)) {
    idord <- nchar(max(dat[[idvar]]))
    for (x in c(netvars, idvar))
      dat[[x]] <- dat[[x]] + dat[[groupvar]]*(10^(idord))
  }

  # Removing unsurveyed
  if (no.unsurveyed) {
    surveyed <- unique(dat[[idvar]])
    for (x in netvars)
      dat[[x]][which(!(dat[[x]] %in% surveyed))] <- NA
  }

  # Reshaping data (so we have an edgelist)
  dat.long <- reshape(
    dat[,c(idvar, netvars)], v.names= "net",
    varying = netvars,
    idvar="id", direction="long")

  # Computing the times of adoption
  rtoa <- range(dat[[toavar]], na.rm = TRUE)
  t    <- rtoa[2] - rtoa[1] + 1

  graph <- with(
    dat.long,
    edgelist_to_adjmat(edgelist = cbind(id, net), t = t,
                       undirected=undirected, self=self, multiple = multiple,
                       keep.isolates = keep.isolates, recode.ids = recode.ids)
  )

  used.vertex <- rownames(graph[[1]])

  colstoexport <- which(!(colnames(dat) %in% c(idvar, toavar)))
  used <- which(dat[[idvar]] %in% used.vertex)
  as_diffnet(
    graph=graph, toa=dat[[toavar]][used],
    vertex.static.attrs = dat[used,colstoexport],
    ...
  )
}
