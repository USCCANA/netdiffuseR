sedist<-function(dat,g=c(1:dim(dat)[1]),method="hamming",joint.analysis=FALSE,mode="digraph",diag=FALSE,code.diss=FALSE){
   #First, prepare the data
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("sedist requires input graphs to be of identical order.")
   if(length(dim(dat))>2){
      n<-dim(dat)[2]
      m<-length(g)
      d<-dat[g,,,drop=FALSE]
   }else{
      n<-dim(dat)[2]
      m<-1
      d<-array(dim=c(m,n,n))
      d[1,,]<-dat
   }
   if(!diag)
      d<-diag.remove(d)
   #Are we conducting a joint analysis?
   if(joint.analysis){
      o<-array(dim=c(1,n,n))
      #Build the data matrix
      v<-vector()
      for(i in 1:n)
         v<-cbind(v,c(as.vector(d[,i,]),as.vector(d[,,i])))
      #Proceed by method
      if(method=="correlation"){
         o[1,,]<-cor(v,use="pairwise")
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="euclidean"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
      }else if(method=="hamming"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
      }else if(method=="gamma"){
         for(i in 1:n)
            for(j in 1:n){
               concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
               discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
               o[1,i,j]<-(concord-discord)/(concord+discord)
            }                  
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="exact"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
      }
   }else{  #Analyze each graph seperately
      o<-array(dim=c(m,n,n))
      for(k in 1:m){
         #Build the data matrix
         v<-vector()
         for(i in 1:n)
            v<-cbind(v,c(as.vector(d[k,i,]),as.vector(d[k,,i])))
         #Proceed by method
         if(method=="correlation"){
            o[k,,]<-cor(v,use="pairwise")
            o[k,,][is.na(o[k,,])]<-0
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="euclidean"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
         }else if(method=="hamming"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
         }else if(method=="gamma"){
            for(i in 1:n)
               for(j in 1:n){
                  concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
                  discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
                  o[k,i,j]<-(concord-discord)/(concord+discord)
               }                  
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="exact"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
         }
      }
   }
   #Dump the output
   if(dim(o)[1]==1)
      as.matrix(o[1,,])
   else
      o
}
