#initialize libraries for creating data
library(foreign)
library(network)
library(cshapes)
library(statnet)

setwd("C:\\Users\\stepharp\\Desktop\\FCTC\\diffusiontest")

netprep<-function(data, varnames=c("year", "id", "alter", "dichot")){
  rawdata <- read.dta(data)
  dyad1 <- cbind(rawdata$year, rawdata$ego_number, rawdata$alter_number, rawdata$dichot)
  dyad1 <- data.frame(dyad1)
  colnames(dyad1) <- varnames
  return(dyad1)
}

ratify2  <- read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\ratify_new.dta")    # here
ratify <- ratify2[order(ratify2$ego_number),]
attributes(ratify)<-NULL
ratify<-data.frame(ratify, stringsAsFactors=FALSE)
names(ratify)<-c("id", "countrynameun", "treaty", "year")
devtools::use_data(ratify, ratify, overwrite = TRUE)

#get distance matrix; cshapes:a GIS dataset of country boundaries<=distance computations on country polygons for specific points in time
worlddate<-as.Date("2003-1-1") #set a bit later to make sure all countries are there
distance_cap<-distmatrix(date=worlddate,type="mindist", useGW=FALSE) #COW has all our countries
#get our codes
codematch<-read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\COWandMatrix.dta") #pulls in our coding scheme

#put our codes on the matrix and drop the country not in our country set
colrownames<-dimnames(distance_cap) #pull out the current names/codes
unlistedcolnames<-matrix(data=as.numeric(unlist(x=colrownames)), nrow=192, ncol=2, byrow=FALSE) #turn current column names into matrix
unlistedcolnames<-data.frame(unlistedcolnames) #make that matrix into dataframe
mergecodematch<-merge(x=unlistedcolnames, y=codematch, by.x= "X1", by.y="COWcode", sort=FALSE, all.x=TRUE) # merge our codes to theirs, adds an NA to our codes
mergereorder<-(merge(x=unlistedcolnames, y=mergecodematch, by.x="X1",by.y="X1", sort=FALSE, all.x=TRUE)) #merge a second time, to reorder due to that NA
newdimnames<-list(mergereorder$matrix_id, mergereorder$matrix_id) #turn our newly ordered codes into a list
i<-1:192
#next 2 lines commented out, they don't need to run everytime
#whatcol<-is.na(newdimnames[[1]][i]) #mark true where that NA is to find it
#print(whatcol) #visually find the NA, 177 is location of NA (this is taiwan)
distance_cap_ourcodes<-matrix(data=distance_cap, dimnames=newdimnames, nrow=192, ncol=192) #apply the list as the new dimnames
distance_cap_ourcodes<-distance_cap_ourcodes[-177,i] #remove NA column from row
j<-1:191
distance_cap_ourcodes<-distance_cap_ourcodes[j,-177] #remove NA from column

#Now turn the matrix into an edgelist
distance_edgelist <-matrix(data=0,36481,3) #make empty matrix to put the edgelist in
row<-1 #start a counter for what row we're filling in
for ( i in 1:191 ){ #for the rows
  for ( j in 1:191) {   #for the columns
    distance_edgelist[row,1]<-as.numeric(dimnames(distance_cap_ourcodes)[[1]][i]) #puts in x row, first column the code for the row
    distance_edgelist[row,2]<-as.numeric(dimnames(distance_cap_ourcodes)[[2]][j]) #puts in x row, second column code for column
    distance_edgelist[row,3]<-distance_cap_ourcodes[i,j] #puts in the number in the matrix for the cell corresponding to the row and column we want, in the x row, third column
    row<-row+1 #go to the next row
  }#loop back for all 191 country matches with i country
}# loop back for all 191 countries
#bind that edgelist to the year 2003 column, add names, call it dyad1
dyad1 <- cbind(1, distance_edgelist)
dyad1 <-data.frame(dyad1)
colnames(dyad1) <- c("year", "id", "alter", "dichot2")
#   write.table(dyad1, file ="dyad_net_1_debug.txt", append=FALSE,
#              quote=FALSE, sep="\t", row.names=F, col.names=TRUE, na=".")

dyad1$dichot[dyad1$dichot2 >= 6529.994] <- 0
dyad1$dichot[dyad1$dichot2 < 6529.994] <- 1
#KF: remove extra column, Dec 2013
dyad1$dichot2 <- NULL
distance<-dyad1
devtools::use_data(distance, distance, overwrite = TRUE)

###network2
trade2  <- read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\trade_2002_2011.dta")
trade   <- trade2[which(trade2$network==2),]
trade   <- trade[which(trade$year>0),]
dyad1 <- cbind(trade$year, trade$ego_number, trade$alter_number, trade$dichot)
dyad1 <-data.frame(dyad1)
colnames(dyad1) <- c("year", "id", "alter", "dichot")
trade<-dyad1
devtools::use_data(trade, trade, overwrite = TRUE)

###network3
trade2  <- read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\trade_tobacco_2002_2011.dta")
trade   <- trade2[which(trade2$year>0),]
dyad1 <- cbind(trade$year, trade$ego_number, trade$alter_number, trade$dichot)
dyad1 <-data.frame(dyad1)
colnames(dyad1) <- c("year", "id", "alter", "dichot")
tobtrade<-dyad1
devtools::use_data(tobtrade, tobtrade, overwrite = TRUE)

#Phone net 4
phone<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\phone_2002_2011.dta")
devtools::use_data(phone, phone, overwrite = TRUE)

# colony.dta net 5
colony<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\colony.dta")
devtools::use_data(colony, colony, overwrite = TRUE)

# INB Network net 6
inbfctc<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\inb_ll.dta")
devtools::use_data(inbfctc, inbfctc, overwrite = TRUE)

# COP Participation Network net 7
copfctc<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\cop_ll.dta")
devtools::use_data(copfctc, copfctc, overwrite = TRUE)

# inb_cac net 8
adhocuncac<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\adhoc_uncac_II.dta")
devtools::use_data(adhocuncac, adhocuncac, overwrite = TRUE)

# cop_corrupt_II.dta net 9
copuncac<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\cop_corrupt_II.dta")
devtools::use_data(copuncac, copuncac, overwrite = TRUE)

# INC_pollute.dta net10
incpollute<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\INC_pollute.dta")
devtools::use_data(incpollute, incpollute, overwrite = TRUE)

# cop_pollute_II.dta net11
coppollute<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\cop_pollute_II.dta")
devtools::use_data(coppollute, coppollute, overwrite = TRUE)

#Referrals net12
rawdata <- read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\refer_full.dta")
dyad1 <- cbind(rawdata$year, rawdata$ego_number, rawdata$alter_number, rawdata$dichot)
dyad1 <- data.frame(dyad1)
colnames(dyad1) <- c("year", "id", "alter", "dichot")
glreferral<-dyad1
devtools::use_data(glreferral, glreferral, overwrite = TRUE)

# Posts Network net13
glposts<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\posts_full.dta")
devtools::use_data(glposts, glposts, overwrite = TRUE)

# Subscription (interest group)Network <= run "two_mode.R" script first to compute C matrix, Affiliation, year constant, year and dichot all ones
# Reminder: copy output file to current directory
#network 14
dyad1<-read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\sub_all.dta")
dyad1 <- data.frame(dyad1)
colnames(dyad1) <- c("year", "id", "alter", "dichot")
glsubscription<- dyad1
devtools::use_data(glsubscription, glsubscription, overwrite = TRUE)

# BIS  (Bank of International Stettlements) network 15
internatsettle<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\bis_dichot.dta")
devtools::use_data(internatsettle, internatsettle, overwrite = TRUE)

# Bilateral Aid Flow net16
bilataid<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\bilat_aid_flow_clean.dta")
devtools::use_data(bilataid, bilataid, overwrite = TRUE)

# Bilateral Investment Treaties network 17
bilatinvesttreaty<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\bit_2.19.14.dta")
devtools::use_data(bilatinvesttreaty, bilatinvesttreaty, overwrite = TRUE)

# IMF Assets network 18
imfassets<-netprep(data="C:\\Users\\stepharp\\Desktop\\FCTC\\Data 9.8.14\\imf.dta")
devtools::use_data(imfassets, imfassets, overwrite = TRUE)

# Random network net 19
rand_mat <- matrix(rbinom((191*191),1,.20),191,191)
diag(rand_mat) <- 0
rand_net <- as.network(rand_mat)
dyad  <- cbind(1, as.edgelist(rand_net), 1)
dyad1 <- data.frame(dyad)
colnames(dyad1) <- c("year", "id", "alter", "dichot")
random<-dyad1
devtools::use_data(random, random, overwrite = TRUE)


#medical innovation

mi.raw<-read.dta("C:\\Users\\stepharp\\Desktop\\FCTC\\classicdiffusiondatasets\\mi_v2.dta")
devtools::use_data(mi.raw, mi.raw, overwrite = TRUE)
