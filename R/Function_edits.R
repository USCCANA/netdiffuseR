
# cumulative adopters function --------------------------------------------


#cumulative adopters function
#'
#'calculates the percent of adopters by time period
#'@param Adopt_mat matrix #nodes x # timeperiods indicating if node adopted or not
#'@return vector indicating the adoption rate and percent adopted by time period
#'
cumulativeAdopterse<- function(Adopt_mat){
  n <- nrow(Adopt_mat)
  maxTime <- ncol(Adopt_mat)
  diff<- (colSums(Adopt_mat) / n)
  rates  <- ((diff[3:(maxTime-1)] - diff[2:(maxTime-2)]) /diff[2:(maxTime-2)])
  rate   <- mean(rates)
  if (is.infinite(rate)) {rate <- NA}
  diff1 <- cbind(rate, t(diff))
  colnames(diff1) <- c("rate", paste("cumadoptperc",(1:(dim(t(diff))[2])),sep="_"))
  return(diff1)
}


# hazard rate function ----------------------------------------------------
##numbering on hazard rate column names, is this considered correct?
##compare new adopts to rates from above cumadopt
#hazard rate function
#'
#'report the hazard rate
#'
#'@param Adopt_mat matrix #nodes x # timeperiods indicating if node adopted or not
#'@return vector of hazard rates
hazardratee<- function(Adopt_mat){
  maxTime <- ncol(Adopt_mat)
  n <- nrow(Adopt_mat)
  cum_num    <- colSums(Adopt_mat)
  new_adopts <- cum_num[3:(maxTime-1)] - cum_num[2:(maxTime-2)]
  pool       <- (n - (colSums(Adopt_mat)))
  pool       <- pool[2:(maxTime-2)]
  hazard     <- (new_adopts / pool)
  hazard1 <- cbind(t(hazard))
  colnames(hazard1) <- c(paste("hazard",(1:(dim(t(hazard))[2])),sep="_"))
  return(hazard1)
}


# Threshold ---------------------------------------------------------------
###threshold currently is only using exposure not exposurese

#Threshold
#'
#'calculates threshold
#'
#'@param exposureObject exposure matrices from ExposureCalc
#'@param toa

thresholde<- function(exposureObject, toa){

  Exposure <- exposureObject$Exposure
  ExposureSE <- exposureObject$ExposureSE
  n <- length(toa)
  maxTime <- ncol(Exposure)
  threshold <- matrix(0, n, 1)
  for (thresh_id in 1:n) {
    threshold[thresh_id]  <-  Exposure[thresh_id, toa[thresh_id]]
  }
  res <-cbind(as.numeric(row.names(Exposure)), threshold)
  colnames(res)<-c("id", "threshold")
  return(res)
}


# event history -----------------------------------------------------------


#event history dataset
#'
#'makes event history database
#'
#'@param exposureObject exposure matrices from ExposureCalc
#'@param toa vector indicating times of adoption
#'@param Adopt_mat
eventhistorye<- function(exposureObject, toa, Adopt_mat){
  Exposure <- exposureObject$Exposure
  ExposureSE <- exposureObject$ExposureSE
  n <- length(toa)
  maxTime <- ncol(Exposure)
  cn <- c("id", "ratify_year", paste("time",(1:maxTime),sep="_"), paste("exp",(1:maxTime),sep="_"),
          paste("expSE",(1:maxTime),sep="_"), paste("adopt",(1:maxTime),sep="_"))
  event <- NULL
  for (i in 1:n) {
    ry <- toa[i]
    padd <- t(rep(NA, maxTime-ry))
    event2 <- cbind(as.numeric(row.names(Adopt_mat)[i]), ry, t(1:ry), padd, t(Exposure[i,1:ry]), padd,
                    t(ExposureSE[i, 1:ry]), padd, t(Adopt_mat[i,1:ry]),padd)
    colnames(event2) <- cn
    if(is.null(event)){
      event  <- event2
    }else{
      event  <- rbind(event, event2)
    }
    event  <- round(event, digits=3)
  }
  return(event)
}


# infect & suscept --------------------------------------------------------
###maybe remove the time column from the final output?
#
#Function to calculate susceptibility and infectiousness
#
#'
#'Calculates susceptibility and infectiousness
#'
#'@param toa timevector
#'@param toamat difference in adoption times between all nodes in network
#'@param Adopt_mat1
#'@param all_nets adjacency matrices
#'@param maxTime number of timepoints
susceptInfecte <- function(toa, toamat, Adopt_mat1, all_nets){
  maxTime<-ncol(Adopt_mat1)
  res <- NULL
  n <- length(toa)
  cum_ados1 <- colSums(Adopt_mat1)
  rowadopt<-matrix(rep(toa,n), n, n)
  for(i in 1:maxTime){
    #Susceptibility

    adjmat_mat <- all_nets[ , , i]

    sus_mat <- as.integer((adjmat_mat + (toamat==-1) + (rowadopt==i)) ==3)
    sus_mat <- matrix(sus_mat, n, n)

    sus_num <- rowSums(sus_mat)
    sus_mat2 <- as.integer((adjmat_mat + (toamat<=-1) + (rowadopt==i)) ==3)
    sus_mat2 <- matrix(sus_mat2, n, n)

    out_deg_before <- rowSums(sus_mat2)
    suscep_1n <- (sus_num / (out_deg_before+.00001))
    suscep_1n_n <- rep(0, n)
    for (u in 1:n) {
      if (toa[u]==1)  {(suscep_1n_n[u] <- suscep_1n[u]) #just be zero?
      }else{
        suscep_1n_n[u] <- (suscep_1n[u] / (cum_ados1[toa[u]-1]+.0001))
      }
    }
    inf_mat <- as.integer((t(adjmat_mat) + (toamat==1) + (rowadopt==i)) ==3)
    inf_mat <- matrix(inf_mat, n, n)

    inf_num <- rowSums(inf_mat)
    inf_mat2 <- as.integer((t(adjmat_mat) + (toamat>=1) + (rowadopt==i)) ==3) # Nominated ego and adopted later
    inf_mat2 <- matrix(inf_mat2, n, n)

    in_deg_after <- rowSums(inf_mat2)
    infect_1n <- (inf_num / (in_deg_after+.00001))
    infect_1n_n <- rep(0, n)
    for (u in 1:n) {
      if (toa[u]==maxTime) {(infect_1n_n[u] <- infect_1n[u])
      }else{
        infect_1n_n[u] <- (infect_1n[u] / (cum_ados1[toa[u]+1]+.0001)) #just be zero?
      }
    }
    suscep_infect <- cbind(i, toa, as.numeric(row.names(Adopt_mat1)), suscep_1n, infect_1n,  suscep_1n_n,
                          infect_1n_n, out_deg_before, in_deg_after)
    colnames(suscep_infect)<-c("time", "toa", "id", "suscep_1n", "infect_1n",  "suscep_1n_n", "infect_1n_n",
                               "out_deg_before", "in_deg_after")
    suscep_infect <- round(suscep_infect, digits=3)
    if(is.null(res)){
      res = suscep_infect
    }else{
      res <- rbind(res, suscep_infect)
    }
  }
  res<-subset(res, subset=(res[,1]==res[,2]))
  res<-res[order(res[,3]),]
  return(res)
}


# Selection Function ------------------------------------------------------------------------

#Selection function and documentation
#
#' Calculate the number of times the 16 selection changes occurred.
#'
#'
#' @param all_nets An array of networks stored as adjacency matrices.
#' @param Adopt_mat Matrix of row=nodes, column=year, Binary if adopted by year.
#' @param Time Timepoint in which calculating selection for.
#' @return Matrix of of the count of selection changes for the 16 categories by node.
selectionFunctionEgoAltere <- function(all_nets, Adopt_mat, Time){
  n <- nrow(all_nets)
  if (Time>1) {
    change_mat <- (all_nets[,,(Time-1)] + all_nets[,,(Time)])
    change_mat_s <- (all_nets[,,(Time-1)] - all_nets[,,(Time)])
    net_stabl  <- matrix(as.integer(change_mat ==2), n, n)
    net_added  <- matrix(as.integer(change_mat_s ==-1), n, n)
    net_dropd  <- matrix(as.integer(change_mat_s ==1), n, n)

    adopts <- Adopt_mat[,Time]
    adoptstm1 <- Adopt_mat[,(Time-1)]
    select_mat <- matrix(0, n, n)

    for (i in 1:n) {
      for (j in 1:n) {
        if (i==j) {next}

        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-1}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-2}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-3}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-4}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-5}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-6}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==0)) {select_mat[i,j] <-7}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==0) & (adopts[j]==1)) {select_mat[i,j] <-8}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-9}
        if ((adoptstm1[i]==0) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-10}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-11}
        if ((adoptstm1[i]==0) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-12}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-13}
        if ((adoptstm1[i]==1) & (adopts[i]==0) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-14}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==0)) {select_mat[i,j] <-15}
        if ((adoptstm1[i]==1) & (adopts[i]==1) & (adoptstm1[j]==1) & (adopts[j]==1)) {select_mat[i,j] <-16}
      }
    }
    select_mat_d  <- (select_mat * net_dropd)
    select_mat_s  <- (select_mat * net_stabl)
    select_mat <- (select_mat * net_added)

    #apply(select_mat, 1, table (x, responseName=select))

    select1 <- rowSums((select_mat) == 1)
    select2 <- rowSums((select_mat) == 2)
    select3 <- rowSums((select_mat) == 3)
    select4 <- rowSums((select_mat) == 4)
    select5 <- rowSums((select_mat) == 5)
    select6 <- rowSums((select_mat) == 6)
    select7 <- rowSums((select_mat) == 7)
    select8 <- rowSums((select_mat) == 8)
    select9 <- rowSums((select_mat) == 9)
    select10 <- rowSums((select_mat) == 10)
    select11 <- rowSums((select_mat) == 11)
    select12 <- rowSums((select_mat) == 12)
    select13 <- rowSums((select_mat) == 13)
    select14 <- rowSums((select_mat) == 14)
    select15 <- rowSums((select_mat) == 15)
    select16 <- rowSums((select_mat) == 16)

    selectd1 <- rowSums((select_mat_d) == 1)
    selectd2 <- rowSums((select_mat_d) == 2)
    selectd3 <- rowSums((select_mat_d) == 3)
    selectd4 <- rowSums((select_mat_d) == 4)
    selectd5 <- rowSums((select_mat_d) == 5)
    selectd6 <- rowSums((select_mat_d) == 6)
    selectd7 <- rowSums((select_mat_d) == 7)
    selectd8 <- rowSums((select_mat_d) == 8)
    selectd9 <- rowSums((select_mat_d) == 9)
    selectd10 <- rowSums((select_mat_d) == 10)
    selectd11 <- rowSums((select_mat_d) == 11)
    selectd12 <- rowSums((select_mat_d) == 12)
    selectd13 <- rowSums((select_mat_d) == 13)
    selectd14 <- rowSums((select_mat_d) == 14)
    selectd15 <- rowSums((select_mat_d) == 15)
    selectd16 <- rowSums((select_mat_d) == 16)

    selects1 <- rowSums((select_mat_s) == 1)
    selects2 <- rowSums((select_mat_s) == 2)
    selects3 <- rowSums((select_mat_s) == 3)
    selects4 <- rowSums((select_mat_s) == 4)
    selects5 <- rowSums((select_mat_s) == 5)
    selects6 <- rowSums((select_mat_s) == 6)
    selects7 <- rowSums((select_mat_s) == 7)
    selects8 <- rowSums((select_mat_s) == 8)
    selects9 <- rowSums((select_mat_s) == 9)
    selects10 <- rowSums((select_mat_s) == 10)
    selects11 <- rowSums((select_mat_s) == 11)
    selects12 <- rowSums((select_mat_s) == 12)
    selects13 <- rowSums((select_mat_s) == 13)
    selects14 <- rowSums((select_mat_s) == 14)
    selects15 <- rowSums((select_mat_s) == 15)
    selects16 <- rowSums((select_mat_s) == 16)

    select <- cbind(time=Time, id=as.numeric(row.names(Adopt_mat)), select1_=select1, select2_=select2,
                    select3_=select3, select4_=select4, select5_=select5, select6_=select6, select7_=select7,
                    select8_=select8, select9_=select9, select10_=select10, select11_=select11,
                    select12_=select12, select13_=select13, select14_=select14, select15_=select15,
                    select16_=select16, selectd1_=selectd1, selectd2_=selectd2, selectd3_=selectd3,
                    selectd4_=selectd4, selectd5_=selectd5, selectd6_=selectd6, selectd7_=selectd7,
                    selectd8_=selectd8, selectd9_=selectd9, selectd10_=selectd10, selectd11_=selectd11,
                    selectd12_=selectd12, selectd13_=selectd13, selectd14_=selectd14, selectd15_=selectd15,
                    selectd16_=selectd16, selects1_=selects1, selects2_=selects2, selects3_=selects3,
                    selects4_=selects4, selects5_=selects5, selects6_=selects6, selects7_=selects7,
                    selects8_=selects8, selects9_=selects9, selects10_=selects10, selects11_=selects11,
                    selects12_=selects12, selects13_=selects13, selects14_=selects14, selects15_=selects15,
                    selects16_=selects16)
    return(select)
  }
}

# Data preparation functions
# Includes: netprep, adjmatbuild, adjByTime, adoptMat, toaMat, ExposureCalc
#

# network data preparation ------------------------------------------------


#netprep function and documentation
#'
#'take a stata dataset and turn it into a linklist within a 4-column dataframe
#'@param data in order rawdata$year, rawdata$ego_number, rawdata$alter_number, rawdata$dicho
#'@param varnames variable names you want the output dataframe to have, default is year, id, alter, dichot
#'@return a linklist in form of dataframe with columns for year and dichot (1 or 0 indicating tie presence)
netprepe<-function(data, varnames=c("year", "id", "alter", "dichot")){
  dyad1 <- data.frame(data)
  colnames(dyad1) <- varnames
  return(dyad1)
}

# adjacency matrix --------------------------------------------------------


#adjacency matrix building function and documentation
#
#working on generalizing this to a linklist with a valued column
#working on making this work for undirected networks
#
#'
#' turns a linklist into an adjacency matrix
#' @param dyad a linklist with 2-3 columns (ego, alter, and optional tie value)
#' @param n
#' @param nodeids the id #s retained in the matrix. These will be indicated by the row and column names
#' @return an adjacency matrix (#nodes by #nodes, value for ties or 1 & 0 throughout)
adjmatbuilde<- function(dyad, n=max(dyad[,1:2]), nodeids=unique(dyad[, 1])){
  if (n < max(nodeids)) {
    stop("The highest node id number cannot be less than the highest node number to be retained.")
  }
  adjmat_mat <- matrix(0,n,n)
  if (ncol(dyad)==3) {
    for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,1], dyad[m,2] ] <- dyad[m,3]
  } else {
    for(m in 1:nrow(dyad)) adjmat_mat[dyad[m,1], dyad[m,2] ] <- 1
  }
  adjmat_mat[is.na(adjmat_mat)]<-0
  diag(adjmat_mat) <- 0
  adjmat_mat<-adjmat_mat[nodeids,nodeids]
  colnames(adjmat_mat)<-nodeids
  row.names(adjmat_mat)<-nodeids
  return(adjmat_mat)
}

#function to build an array of adjacency matrices representing a dynamic network
#
# Do we want to have this be a datafame/ matrix thing or have people specify the vectors for everything?
#Like, specify a vector for id, alter, time, and tie?
#
#'
#'create an array of adjacency matrices from dynamic data
#'
#'@param data a 2 - 4 column data frame or matrix with columns ordered with for timepoint, id, alter, and tie
#'value (optional). If data is static, a 2-3 column dataset is expected with columns ordered id, alter, and
#'tie value (optional).
#'@param n The highest node id number in the network (numeric). Default is the maximum id present in the id or
#'alter column (columns 2 & 3).
#'@param Time Time is the number of timepoints you have data for (numeric, whole numbers). Default is the
#'maximum of column 1 of your data
#'@param type kind of network data (string) ---static, dynamic, or cumulative (tie at current timepoint if it
#'was present at any previous timepoint)
#'@param nodeids Numeric vector of id numbers to be retained
#'@return Returns a #nodes x #nodes x #timepoints array of adjacency matrices
adjByTimee <- function(data, n = max(data[,2:3]), Time = max(data[,1]), type="dynamic", nodeids=unique(data[, 2])){
  if (n < max(nodeids)) {
    stop("The highest node id number cannot be less than the highest node number to be retained.")
  }
  all_nets <-array(0, dim=c(length(nodeids),length(nodeids),Time))
  if (type=="dynamic"){
    if (mean(data[,1]) != 1){
      for(i in 1:Time){
        dyad  <- data[which(data[,1]==i),2:ncol(data)]
        adjmat_mat<-adjmatbuild(dyad=dyad, n=n, nodeids=nodeids)
        all_nets[,,i] <- adjmat_mat
      }
    }else if (mean(data[,1]) == 1){ dyad = data[,2:ncol(data)]
      adjmat_mat<-adjmatbuild(dyad=dyad, n=n, nodeids=nodeids)
      all_nets[,,] <- adjmat_mat
    }
  }
  if (type=="cumulative"){
    for(i in 1:Time){
      dyad<-data[which(data[,1]<=i), 2:ncol(data)]
      adjmat_mat<-adjmatbuild(dyad=dyad, n=n, nodeids=nodeids)
      all_nets[,,i] <- adjmat_mat
    }
  }
  if (type=="static"){
    dyad<-data
    adjmat_mat<-adjmatbuild(dyad=dyad, n=n, nodeids=nodeids)
    all_nets[,,] <- adjmat_mat
  }
  res <- all_nets
  dimnames(res)<-list(id=nodeids, alter=nodeids, wave=paste("Time", 1:Time, sep="_"))
  return(res)
}


# adoption matrices -------------------------------------------------------



#adoptMat function to create adoption matrix, and documentation
#'
#'creates matrices indicating the adoption of innovation by year
#' @return A list of matrices of dim #nodes x # years indicating for each node what year they adopted
#' an innovation by. Adopt_mat has 1's for all years in which a node indicates having the innovation.
#' Adopt_mat1 has a zero only for the year of adoption and zeros for all other years.
#' @param timevec Numeric vector of year of adoption in order of node ID. Must have all nodes in it, whether
#' adopted or not.
adoptMate <- function(timevec, id=FALSE){
  maxYear <- max(timevec)
  n <- length(timevec)
  Adopt_mat  <- matrix(data = 0, nrow = n, ncol = (maxYear), byrow = TRUE)
  Adopt_mat1  <- matrix(data = 0, nrow = n, ncol = (maxYear), byrow = TRUE)     # for infect & suscep

  for (i in 1:n) {
    Adopt_mat[i,timevec[i]:maxYear] <- 1
    Adopt_mat1[i,timevec[i]] <- 1
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
toaMate <- function(timevec, id){
  n <- NROW(timevec)
  toamat <- matrix(0, n, n)
  row.names(toamat)<-id
  colnames(toamat)<-id
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {next}
      toamat[i,j] <- timevec[j] - timevec[i]   #indicates how sooner or later i adopted relative to j
    }
  }
  return(toamat)
}


# exposure ----------------------------------------------------------------

#might we ever only want to calculate exposure for one timepoint?
#exposurecalc function and documentation.
#'
#'exposurecalc uses network and innovation adoption data to calculate exposure to the innovation for each node
#'
#'@param all_nets Adjacency matrices for the dynamic network.
#'@param Adopt_mat Matrix indicating if a node has adopted an innovation
#'@param method Method to use for calculating exposure, options "direct", "twostep", "threestep", or "se"
#'(structural equivalence)
#'@return A matrix with number of rows = number of nodes and number of columns = number of timepoints
#'indicating exposure for each node at each timepoint.
ExposureCalce <- function(all_nets, Adopt_mat, method="direct"){
  n <- dim(all_nets)[1]
  maxTime <- ncol(Adopt_mat)
  Exposure   <- matrix(data = 0, nrow = n, ncol = (maxTime), byrow = TRUE)
  row.names(Exposure)=row.names(Adopt_mat)
  if (method=="direct"){
    for(Time in 1:(maxTime-1)){
      adjmat_mat <- all_nets[,,Time]
      Exposure[,Time]   <-((adjmat_mat %*% Adopt_mat[,Time]) / (rowSums(adjmat_mat)+.0001))
    }
  }
  if (method=="twostep"){
    for(Time in 1:(maxTime-1)){
      adjmat_mat <- all_nets[,,Time]
      twostep <- matrix(NA,nrow(adjmat_mat),nrow(adjmat_mat))
      for (row in 1:nrow(adjmat_mat)){
        x<-adjmat_mat[row,]
        y<-names(x)[x==1]
        if (length(y)>1) twostep[row,]<- x + colSums(adjmat_mat[y,])
        if (length(y)==1) twostep[row,]<- x + adjmat_mat[y,]
        if (length(y)==0) twostep[row,]<- x
      }
      twostep[twostep>1]<-1
      diag(twostep) <- 0
      Exposure[,Time]   <-((twostep %*% Adopt_mat[,Time]) / (rowSums(twostep)+.0001))
    }
  }
  if (method=="threestep"){
    for(Time in 1:(maxTime-1)){
      adjmat_mat <- all_nets[,,Time]
      threestep <- matrix(NA,nrow(adjmat_mat),nrow(adjmat_mat))
      for (row in 1:nrow(adjmat_mat)){
        x<-adjmat_mat[row,]
        y<-names(x)[x==1]
        if (length(y)>1) twosteprow<- x + colSums(adjmat_mat[y,])
        if (length(y)==1) twosteprow<- x + adjmat_mat[y,]
        if (length(y)==0) twosteprow<- x
        y2<-names(x)[twosteprow>=1]
        if (length(y2)>1) threestep[row,]<- twosteprow + colSums(adjmat_mat[y2,])
        if (length(y2)==1) threestep[row,]<- x + adjmat_mat[y2,]
        if (length(y2)==0) threestep[row,]<- x
      }
      threestep[threestep>1]<-1
      diag(threestep) <- 0
      Exposure[,Time]   <-((threestep %*% Adopt_mat[,Time]) / (rowSums(threestep)+.0001))
    }
  }
  if (method=="se"){
    for(Time in 1:(maxTime-1)){
      adjmat_mat <- all_nets[,,Time]
      semat_mat <- (1 /(sedist(adjmat_mat, method="euclidean")))
      semat_mat[is.infinite(semat_mat[,])]<- 0
      diag(semat_mat) <- 0
      Exposure[,Time] <-((semat_mat %*% Adopt_mat[,Time]) / (rowSums(semat_mat)+.0001))
    }
  }
  return(Exposure)
}
