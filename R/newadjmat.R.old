##this was just a try at something different
adjByTime2 <- function(data, n = n, Time = 10, type="dynamic"){
  all_nets <-array(0, dim=c(n,n,Time))
  if (type=="dynamic"){
    if (mean(data$year) != 1){
      listadj<-lapply(split(data, data$year), function(x) {
        x$year<-NULL
        adjmatbuild(dyad=x, n=n)
      })
      all_nets<-array(as.vector(listadj), dim=c(191,191,10))
    }else if (mean(data$year) == 1){
    data$year <- NULL
    adjmat_mat<-adjmatbuild(dyad=data, n=n)
    all_nets[,,] <- adjmat_mat
  }
}
if (type=="cumulative"){
  for(i in 1:Time){
    dyad<-data[which(data$year<=i),]
    dyad$year <- NULL
    adjmat_mat<-adjmatbuild(dyad=dyad, n=n)
    all_nets[,,i] <- adjmat_mat
  }
}
if (type=="static"){
  dyad<-data
  dyad$year <- NULL
  adjmat_mat<-adjmatbuild(dyad=dyad, n=n)
  all_nets[,,] <- adjmat_mat
}
res <- all_nets
return(res)
}
