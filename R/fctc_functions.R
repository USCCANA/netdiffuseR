
# cumulative adopters function --------------------------------------------


#cumulative adopters function
#'
#'calculates the percent of adopters by time period
#'@param Adopt_mat matrix #nodes x # timeperiods indicating if node adopted or not
#'@return vector indicating the adoption rate and percent adopted by time period
#'
cumulativeAdopters<- function(Adopt_mat){
  n <- nrow(Adopt_mat)
  maxTime <- ncol(Adopt_mat)
  diff<- (colSums(Adopt_mat) / n)
  rates  <- ((diff[3:(maxTime-1)] - diff[2:(maxTime-2)]) /diff[2:(maxTime-2)]) # Calcualte rates on 2 - 10 interval
  rate   <- mean(rates)
  if (is.infinite(rate)) {rate <- NA}
  diff1 <- cbind(rate, t(diff))
  colnames(diff1) <- c("rate", paste("cumadoptperc",(1:(dim(t(diff))[2])),sep="_"))
  return(diff1)
}


# hazard rate function ----------------------------------------------------


#hazard rate function
#'
#'report the hazard rate
#'
#'@param Adopt_mat matrix #nodes x # timeperiods indicating if node adopted or not
#'@return vector of hazard rates
hazardrate<- function(Adopt_mat){
  maxTime <- ncol(Adopt_mat)
  n <- nrow(Adopt_mat)
  cum_num    <- colSums(Adopt_mat)
  new_adopts <- cum_num[3:(maxTime-1)] - cum_num[2:(maxTime-2)]  # Calcualte hazard rates on 2 - 10 interval
  pool       <- (n - (colSums(Adopt_mat)))
  pool       <- pool[2:(maxTime-2)]
  hazard     <- (new_adopts / pool)
  hazard1 <- cbind(t(hazard))
  colnames(hazard1) <- c(paste("hazard",(1:(dim(t(hazard))[2])),sep="_"))
  return(hazard1)
}


# Threshold ---------------------------------------------------------------


#Threshold
#'
#'calculates threshold
#'
#'@param exposureObject exposure matrices from ExposureCalc
#'@param toa
#'@param net_no number assigned to network (might remove)
#'@param conv convention number (might remove)
threshold<- function(exposureObject, toa, net_no, conv){

  Exposure <- exposureObject$Exposure
  ExposureSE <- exposureObject$ExposureSE
  n <- length(toa)
  maxTime <- ncol(Exposure)

  threshold_netconv <- matrix(0, n, 1)
  for (thresh_id in 1:n) {
    threshold_netconv[thresh_id]  <-  Exposure[thresh_id, toa[thresh_id]]
  }
  threshold <-  cbind(1:n, threshold_netconv)

  return(threshold)
}


# event history -----------------------------------------------------------


#event history dataset
#'
#'makes event history database
#'
#'@param exposureObject exposure matrices from ExposureCalc
#'@param toa vector indicating times of adoption
#'@param Adopt_mat
#'@param net_no number assigned to network (might remove)
#'@param conv convention number (might remove)
eventhistory<- function(exposureObject, toa, Adopt_mat, net_no, conv){

  Exposure <- exposureObject$Exposure
  ExposureSE <- exposureObject$ExposureSE
  ExposureC <- exposureObject$ExposureC
  n <- length(toa)
  maxTime <- ncol(Exposure)

  # Make Event History Dataset
  cn <- c("conv", "net_no", "id", "ratify_year", paste("time",(1:maxTime),sep="_"), paste("exp",(1:maxTime),sep="_"), paste("expSE",(1:maxTime),sep="_"), paste("expC",(1:maxTime),sep="_"), paste("adopt",(1:maxTime),sep="_"))

  event <- NULL

  for (i in 1:n) {
    ry <- toa[i]
    padd <- t(rep(NA, maxTime-ry))
    event2 <- cbind(conv, net_no, i, ry, t(1:ry), padd, t(Exposure[i,1:ry]), padd,
                    t(ExposureSE[i, 1:ry]), padd, t(ExposureC[i,1:ry]), padd, t(Adopt_mat[i,1:ry]),padd)
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
susceptInfect <- function(toa, toamat, Adopt_mat1, all_nets, maxTime){
  res <- NULL
  n <- length(toa)
  cum_ados1 <- colSums(Adopt_mat1)
  rowadopt<-matrix(rep(toa,191), n, n)
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
    suscep_infect <- cbind(i, toa, 1:191, suscep_1n, infect_1n,  suscep_1n_n, infect_1n_n, out_deg_before, in_deg_after)
    colnames(suscep_infect)<-c("time", "toa", "id", "suscep_1n", "infect_1n",  "suscep_1n_n", "infect_1n_n", "out_deg_before", "in_deg_after")
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
#' @param Time Year in which calculating selection for.
#' @return Matrix of of the count of selection changes for the 16 categories by node.
selectionFunctionEgoAlter <- function(all_nets, Adopt_mat, Time){
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

    select <- cbind(time=Time, id=1:n, select1_=select1, select2_=select2, select3_=select3, select4_=select4, select5_=select5, select6_=select6, select7_=select7, select8_=select8, select9_=select9, select10_=select10,
                    select11_=select11, select12_=select12, select13_=select13, select14_=select14, select15_=select15, select16_=select16,
                    selectd1_=selectd1, selectd2_=selectd2, selectd3_=selectd3, selectd4_=selectd4, selectd5_=selectd5, selectd6_=selectd6, selectd7_=selectd7, selectd8_=selectd8, selectd9_=selectd9, selectd10_=selectd10,
                    selectd11_=selectd11, selectd12_=selectd12, selectd13_=selectd13, selectd14_=selectd14, selectd15_=selectd15, selectd16_=selectd16,
                    selects1_=selects1, selects2_=selects2, selects3_=selects3, selects4_=selects4, selects5_=selects5, selects6_=selects6, selects7_=selects7, selects8_=selects8, selects9_=selects9, selects10_=selects10,
                    selects11_=selects11, selects12_=selects12, selects13_=selects13, selects14_=selects14, selects15_=selects15, selects16_=selects16)
    return(select)
  }
}


