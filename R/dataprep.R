# Data preparation functions
# Includes: netprep, adjmatbuild, adjByTime, adoptMat, toaMat, ExposureCalc
#

# network data preparation ------------------------------------------------


#netprep function and documentation
#'
#'take a stata dataset and turn it into a linklist within a 4-column dataframe
#'@param data a dataset in stata format with columns for year, ego_number, alter_number, and dichot (tie or no tie)
#'@param varnames variable names you want the output dataframe to have, default is year, id, alter, dichot
#'@return a linklist in form of dataframe with columns for year and dichot (1 or 0 indicating tie presence)
netprep<-function(data, varnames=c("year", "id", "alter", "dichot")){
  rawdata <- read.dta(data)
  dyad1 <- cbind(rawdata$year, rawdata$ego_number, rawdata$alter_number, rawdata$dichot)
  dyad1 <- data.frame(dyad1)
  colnames(dyad1) <- varnames
  return(dyad1)
}


# adjacency matrix --------------------------------------------------------


#adjacency matrix building function and documentation
#
#working on generalizing this to a linklist with a valued column
#work on getting this better with the directed stuff
#
#'
#' turns a linklist into an adjacency matrix
#' @param dyad a linklist with 3 columns (ego, alter, tie value (or 1 and 0 if tie is only present or not))
#' @param n n is the number of nodes in the network. Default value is n.
#' @param directed a true/false statement indicating if a network is directed or undirected.
#' @return an adjacency matrix (#nodes by #nodes, value for ties or 1 & 0 throughout)
adjmatbuild<- function(dyad, n=n, nodeids, directed=TRUE){
  adjmat_mat <- matrix(0,n,n)
  if (ncol(dyad)==3) {
    for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,1], dyad[m,2] ] <- dyad[m,3]
  } else {
    for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,1], dyad[m,2] ] <- 1
  }
  if (directed==FALSE) {
    if (ncol(dyad)==3) {
      for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,2], dyad[m,1] ] <- dyad[m,3]
    } else {
      for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,2], dyad[m,1] ] <- 1
    }
  }
  adjmat_mat[is.na(adjmat_mat)]<-0
  diag(adjmat_mat) <- 0
  adjmat_mat<-adjmat_mat[nodeids,nodeids]
  return(adjmat_mat)
}

#function to build an array of adjacency matrices representing a dynamic network
#
#'
#'create an array of adjacency matrices from dynamic data
#'
#'@param data 4-column data frame with year, ego, alter, and tie value (or 0 & 1)
#'@param n number of nodes in network. Default is n.
#'@param Time Time is the number of timepoints you have data for
#'@param type kind of network data---static, dynamic, or cumulative (tie at current timepoint if it was present at any previous timepoint)
#'@param directed a true/false statement indicating if a network is directed or undirected.
#'@return Returns an #nodes x #nodes x #timepoints array of adjacency matrices
adjByTime <- function(data, n = n, Time = 10, type="dynamic", directed=TRUE){
  all_nets <-array(0, dim=c(n,n,Time))
  nodeids<-data[, 1]
if (type=="dynamic"){
    if (mean(data$year) != 1){
     for(i in 1:Time){
       dyad  <- data[which(data$year==i),]
       dyad$year <- NULL
       adjmat_mat<-adjmatbuild(dyad=dyad, n=n, directed=directed)
       all_nets[,,i] <- adjmat_mat
     }
      }else if (mean(data$year) == 1){ dyad = data
       dyad$year <- NULL
       adjmat_mat<-adjmatbuild(dyad=dyad, n=n, directed=directed)
       all_nets[,,] <- adjmat_mat
    }
  }
  if (type=="cumulative"){
    for(i in 1:Time){
      dyad<-data[which(data$year<=i),]
      dyad$year <- NULL
      adjmat_mat<-adjmatbuild(dyad=dyad, n=n, directed=directed)
      all_nets[,,i] <- adjmat_mat
    }
  }
  if (type=="static"){
    dyad<-data
    dyad$year <- NULL
    adjmat_mat<-adjmatbuild(dyad=dyad, n=n, directed=directed)
    all_nets[,,] <- adjmat_mat
  }
  res <- all_nets
  return(res)
}


# adoption matrices ------------------------------for(int j=0;j<n;j++)
for(int k=0;k<n;k++)
  NUMERATOR(j) = dynadjmat(j,k,i)/sedist(j,k)*adopt(k,i);-------------------------



#adoptMat function to create adoption matrix, and documentation
#'
#'creates matrices indicating the adoption of innovation by year
#' @return A list of matrices of dim #nodes x # years indicating for each node what year they adopted
#' an innovation by. Adopt_mat has 1's for all years in which a node indicates having the innovation.
#' Adopt_mat1 has a zero only for the year of adoption and zeros for all other years.
#' @param yearObject vector of year of adoption in order of node ID. Must have all nodes in it, whether adopted or not.
adoptMat <- function(yearObject, id=FALSE){
  maxYear <- max(yearObject)
  n <- length(yearObject)
  Adopt_mat  <- matrix(data = 0, nrow = n, ncol = (maxYear), byrow = TRUE)
  Adopt_mat1  <- matrix(data = 0, nrow = n, ncol = (maxYear), byrow = TRUE)     # for infect & suscep

  for (i in 1:n) {
    Adopt_mat[i,yearObject[i]:maxYear] <- 1
    Adopt_mat1[i,yearObject[i]] <- 1
  }
  if (id[1] !=FALSE){
    row.names(Adopt_mat)<-id
    row.names(Adopt_mat1)<-id
  }
  res <- list(Adopt_mat = Adopt_mat, Adopt_mat1 = Adopt_mat1)
  return(res)
}


# time to adoption matrix -------------------------------------------------


#Function to make time to adoption matrix
#'
#'Creates #nodes x #nodes matrix indicating the difference in times of adoption between each pair of nodes
#'
#'@param timevec Vector indicating the timepoint in which each node adopted the innovation
#'@return #nodes x #nodes matrix indicating the difference in times of adoption between each pair of nodes
## TODO: should this be changed into a matrix with node names and times instead of just a vector of times?
toaMat <- function(timevec){
  toa <- timevec
  n <- NROW(toa)
  toamat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {next}
      toamat[i,j] <- toa[j] - toa[i]   #indicates how sooner or later i adopted relative to j
    }
  }
  return(toamat)
}


# exposure ----------------------------------------------------------------



#exposurecalc function and documentation.
#
# Add in exposure0 at time, not time+1
###also....if (time < (T)) {Exposure1[,time+1]   <-((adjmat_mat %*% Adopt_mat1[,time]) / (rowSums(adjmat_mat)+.0001))}
###do I want to have this output all 3 exposure matrices or separatly specify each one?
#'
#'exposurecalc uses network and innovation adoption data to calculate exposure to the innovation for each node
#'
#'@param all_nets Adjacency matrices for the dynamic network.
#'@param Adopt_mat Matrix indicating if a node has adopted an innovation
#'@param Time The timepoint that you are calculating exposure for
#'@return A list of 2 exposure matrices, one through structural equivalence and the other through adjacency
ExposureCalc <- function(all_nets, Adopt_mat){
  n <- dim(all_nets)[1]
  maxTime <- ncol(Adopt_mat)
  Exposure   <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
  ExposureSE <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
  ExposureC  <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)

  for(Time in 1:maxTime){
    #Anything that is multiplied by adoption is specific by convention
    if(Time < maxTime){
      adjmat_mat <- all_nets[,,Time]
      semat_mat <- (1 /(sedist(adjmat_mat, method="euclidean")))
      semat_mat[is.infinite(semat_mat[,])]<- 0
      diag(semat_mat) <- 0
      in_deg <- (degree(as.network(adjmat_mat), cmode="indegree"))
      ExposureC[,Time+1]    <-((adjmat_mat %*% (Adopt_mat[,Time] * in_deg)) / (rowSums(adjmat_mat)+.0001))
      Exposure[,Time+1]   <-((adjmat_mat %*% Adopt_mat[,Time]) / (rowSums(adjmat_mat)+.0001))
      ExposureSE[,Time+1] <-((semat_mat %*% Adopt_mat[,Time]) / (rowSums(semat_mat)+.0001))
    }
  }

  res <- list(Exposure = Exposure, ExposureSE = ExposureSE, ExposureC = ExposureC)
  return(res)
}
