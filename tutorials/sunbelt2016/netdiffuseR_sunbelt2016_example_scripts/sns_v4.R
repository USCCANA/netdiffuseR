# Reading it into netdiffuseR --------------------------------------------------

rm(list=ls())
library(netdiffuseR)
library(foreign)

# Step 0: Reloading the data -

mainset <- read.dta("playground/20160323_sns/SNS datamerged_HS030316.dta",
                   convert.factors = FALSE, convert.underscore = TRUE)
# mainset <- read.dta("c:/misc/rc/SNS datamerged_HS030316.dta",
#                     convert.factors = FALSE, convert.underscore = TRUE)

# Fixing the colname eversmk -> t[1:4]_j1
cn <- colnames(mainset)
test <- which(grepl("^eversmk[0-9]$",cn))
colnames(mainset)[test] <-
  sprintf("t%s.j1",gsub("[a-z]+", "",cn[test]))
colnames(mainset)[colnames(mainset) == "t4.j1"] <- "t4.nj1"

# Step 1: Putting it into long format ------------------------------------------
varying <- c(
  # Network variables
  lapply(seq_len(7), function(x) sprintf("sch.friend%d%d", 1:4, x)), # School friend
  lapply(seq_len(7), function(x) sprintf("sch.admire%d%d", 1:4, x)),  # School admire
  lapply(seq_len(7), function(x) sprintf("sch.succeed%d%d", 1:4, x)), # School succed
  lapply(seq_len(7), function(x) sprintf("sch.popular%d%d", 1:4, x)), # School popular
  # Tobacco and Alcohol: The following variables are not present in the dataset
  # "t4.nj11" "t4.nj23" "t1.j24"  "t2.j24"  "t3.j24"
  lapply(c(1:10,12:22), function(x) c(sprintf("t%d.j%d", 1:3 ,x), sprintf("t4.nj%d",x)))
)

v.names <- c(
  list(sprintf("sch.friend%02d",1:7)),
  list(sprintf("sch.admire%02d",1:7)),
  list(sprintf("sch.succeed%02d",1:7)),
  list(sprintf("sch.popular%02d",1:7)),
  list(sprintf("tj%02d",c(1:10,12:22)))
)

# Checking which are not in
test <- unlist(varying, recursive = TRUE)
test <- test[which(!(test %in% colnames(mainset)))]

# Stop if no match!
stopifnot(!length(test))

# Step 2: Creating covariates --------------------------------------------------

# Adoption

# J1. Have you ever tried cigarette smoking,
# even one or two puffs?
mainset$toa_smoke <- NA
mainset$toa_smoke[with(mainset, is.na(toa_smoke) & t1.j1==1)] <- 2010
mainset$toa_smoke[with(mainset, is.na(toa_smoke) & t2.j1==1)] <- 2011
mainset$toa_smoke[with(mainset, is.na(toa_smoke) & t3.j1==1)] <- 2012
mainset$toa_smoke[with(mainset, is.na(toa_smoke) & t4.nj1==1)] <- 2013

# J10. During your life, on how many days have you had at least one drink of alcohol?
# (please do not count drinking alcohol for religious purposes like communion wine)
# 1 = 0 days (never drink)
mainset$toa_drink <- NA
mainset$toa_drink[with(mainset, is.na(toa_drink) & t1.j10!=1)] <- 2010
mainset$toa_drink[with(mainset, is.na(toa_drink) & t2.j10!=1)] <- 2011
mainset$toa_drink[with(mainset, is.na(toa_drink) & t3.j10!=1)] <- 2012
mainset$toa_drink[with(mainset, is.na(toa_drink) & t4.nj10!=1)] <- 2013

# Age
# library(lubridate)
# mainset$birthday[mainset$birthday==""] <- NA
# mainset$bday <- lubridate::mdy(mainset$birthday)

# One failed:
#      bday   birthday school photoid
# 1702   NA 21/30/1994    113      48
#View(subset(mainset, !is.na(birthday) & is.na(bday), select=c(bday, birthday, school, photoid)))

# Recoding age
mainset$age        <- (12:18)[mainset$t1.i1]
# mainset$age_l1   <- (12:18)[mainset$t1.i1]
# mainset$age_bday <- floor((ymd(20100615) - mainset$bday)/365)
# mainset$age <- with(mainset, ifelse(is.na(age_l1), age_bday, age_l1))
# with(mainset, View(data.frame(age_l1, age_bday, age)))

# Gender: Should be coded as i3 but it shows as i2
mainset$female <- mainset$t1.i2 == 1

# Lunch
lunch <- c("t1.h1", "t2.h2", "t3.h1", "t4.h1")

#Academic Grades
grades <- c("t1.i11", "t2.i11", "t3.i11", "t4.ni11")

for (i in grades) {
  mainset[[i]][is.na(mainset[[i]])] <- 10L
  mainset[[i]] <- 10 - mainset[[i]]
}

# w1.grade <- mainset$t1.i11
# #replace missing with "I wasn't in school last year (10)"
# w1.grade[is.na(w1.grade) | w1.grade == 10] <- mainset$t1.i11[is.na(w1.grade) | w1.grade == 10]
# w1.grade <- 10-w1.grade #re-coded so higher score means higher grade
#
# w4.grade <- mainset$t4.ni11
# #replace missing with "I wasn't in school last year (10)"
# w4.grade[is.na(w4.grade) | w4.grade == 10] <- mainset$t4.ni11[is.na(w4.grade) | w4.grade == 10]
# w4.grade <- 10-w4.grade #re-coded so higher score means higher grade

nrooms <- c("t1.i10", "t2.i10", "t3.i10", "t4.ni10")
# ####Room per HH member
# w1.rooms <- mainset$t1.i10
# w4.rooms <- mainset$t4.ni10

health <- c("t1.i20", "t2.i20", "t3.i20", "t4.ni13")
mainset$t1.i20[mainset$t1.i20==6] <- NA

# ####Self-reported health
# w1.health <- mainset$t1.i20
# w1.health [w1.health == 6] <- NA #one out of range (range=1~5)
# w4.health <- mainset$t4.ni13

parent_smoke <- c("t1.j6", "t2.j6", "t3.j6", "t4.nj7")
# ####Parent smoke
# w1.parsmoke <- mainset$t1.j6
# w4.parsmoke <- mainset$t4.nj7


sibling_smoke <- c("t1.j7", "t2.j7", "t3.j7", "t4.nj8")
# ####Sibling smoke
# w1.sibsmoke <- mainset$t1.j7
# w4.sibsmoke <- mainset$t4.nj8

parent_drink <- c("t1.j14", "t2.j14", "t3.j14", "t4.nj18")

# ####Parent drink
# w1.pardrink <- mainset$t1.j14
# w4.pardrink <- mainset$t4.nj18

sibling_drink <- c("t1.j19", "t2.j19", "t3.j19", "t4.nj19")
# Recoding accordingly
mainset$t4.nj19 <- with(
  mainset, ifelse(t4.nj19 > 3, 1, ifelse(t4.nj19 < 3, 2, 3)))

varying <- c(varying, list(lunch), list(grades), list(nrooms), list(health),
             list(parent_smoke), list(sibling_smoke),
             list(parent_drink), list(sibling_drink))

v.names <- c(v.names, "lunch", "grades", "nrooms", "health", "parent_smoke",
             "sibling_smoke", "parent_drink", "sibling_drink")

# ####Sibling drink!!!!!!!!!!NOTE:answer options were different between wave1 & wave4-re-coding needed
# w1.sibdrink <- mainset$t1.j19
# w4.sibdrink <- mainset$t4.nj19

vars_static  <- c("photoid", "school", "birthday",
                  "asian", "black", "amerind", "hispanic", "white", "multi", "noethnic", "pacific",
                  "race", "age", "female",
                  "toa_smoke", "toa_drink")

# Reshaping
mainset_long <- reshape(
  subset(mainset, select=c(vars_static, unlist(varying))),
  varying = varying, v.names=unlist(v.names),
  sep="", direction = "long", times=2010:2013
)

# Step 3: Creating the diffnet objects -----------------------------------------

# Friends smoking
diffnet_smoke_friends <- survey_to_diffnet(
  mainset_long,
  idvar    = "photoid",
  netvars  = v.names[[1]],
  toavar   = "toa_smoke",
  timevar  = "time", groupvar = "school",
  warn.coercion = FALSE
)

# Friends drinking
diffnet_drink_friends <- survey_to_diffnet(
  mainset_long,
  idvar    = "photoid",
  netvars  = v.names[[1]],
  toavar   = "toa_drink",
  timevar  = "time", groupvar = "school",
  warn.coercion = FALSE
)



# Step 4: Models ---------------------------------------------------------------
# Here I'm only using: friends network

# Computing exposures
# indexes <- sample.int(nnodes(diffnet_drink_friends_admire), 1e3, TRUE)
# diffnet_drink_friends_admire <- diffnet_drink_friends_admire[indexes]

schools <- diffnet_drink_friends[["school"]][[1]]

# # SE exposure
# diffnet_drink_friends[["seexp"]] <- exposure(
#   diffnet_drink_friends_admire, alt.graph = "se", groupvar=schools,
#   valued=TRUE
# )

# Cohesive exposure
diffnet_smoke_friends[["cohexp_smoke"]] <- exposure(
  diffnet_smoke_friends, valued=FALSE
)
# Degree weighted exposure
diffnet_smoke_friends[["cohexp_degree_smoke"]] <- exposure(
  diffnet_smoke_friends, attrs=dgr(diffnet_smoke_friends), valued=FALSE
)

# Cohesive exposure
diffnet_drink_friends[["cohexp_drink"]] <- exposure(
  diffnet_drink_friends, valued=FALSE
)
# Degree weighted exposure
diffnet_drink_friends[["cohexp_degree_drink"]] <- exposure(
  diffnet_drink_friends, attrs=dgr(diffnet_drink_friends), valued=FALSE
)

# Building dataset from scratch
#dat <- diffnet.attrs(diffnet_drink_friends)
dat <- diffnet.attrs(diffnet_smoke_friends)
sub.dat <- data.frame(
  ever_smk4   = dat[[4]]$tj01 == 1,
  ever_smk1   = dat[[1]]$tj01 == 1,
  ever_drink4 = dat[[4]]$tj10 != 1,
  ever_drink1 = dat[[1]]$tj10 != 1,
  age4      = dat[[4]]$age,
  female4   = dat[[4]]$female,
  lunch4    = dat[[4]]$lunch,
  grades4   = dat[[4]]$grades,
  health4   = dat[[4]]$health,
  rooms4    = dat[[4]]$nrooms,
  parent_smoke4  = dat[[4]]$parent_smoke,
  sibling_smoke4 = dat[[4]]$sibling_smoke,
#  seexp4     = dat[[4]]$seexp,
  school      = dat[[4]]$school,
  cohexp_smoke4 = dat[[4]]$cohexp_smoke,
  cohexp_smoke1 = dat[[1]]$cohexp_smoke,
  cohexp_deg_smoke4 = dat[[4]]$cohexp_degree_smoke
#   cohexp_drink4 = dat[[4]]$cohexp_drink,
#   cohexp_drink1 = dat[[1]]$cohexp_drink
  )

dat <- diffnet.attrs(diffnet_drink_friends)
sub.dat.d <- data.frame(
  ever_smk4   = dat[[4]]$tj01 == 1,
  ever_smk1   = dat[[1]]$tj01 == 1,
  ever_drink4 = dat[[4]]$tj10 != 1,
  ever_drink1 = dat[[1]]$tj10 != 1,
  age4      = dat[[4]]$age,
  female4   = dat[[4]]$female,
  lunch4    = dat[[4]]$lunch,
  grades4   = dat[[4]]$grades,
  health4   = dat[[4]]$health,
  rooms4    = dat[[4]]$nrooms,
  parent_smoke4  = dat[[4]]$parent_smoke,
  sibling_smoke4 = dat[[4]]$sibling_smoke,
  #  seexp4     = dat[[4]]$seexp,
  school      = dat[[4]]$school,
  cohexp_drink4 = dat[[4]]$cohexp_drink,
  cohexp_drink1 = dat[[1]]$cohexp_drink,
  cohexp_deg_drink4 = dat[[4]]$cohexp_degree_drink
)
# Checking complete... it turns out the drops rooms4 b/c singularities... it
# means that is perfectly correlated with some other variable!
which_complete <- complete.cases(sub.dat)
table(sub.dat$rooms4[which_complete])

# Smoking behavior
mod_smoke <- formula(ever_smk4 ~ factor(school) + age4 + female4 + lunch4 + grades4 + health4 + rooms4 +
                 parent_smoke4 + sibling_smoke4 + ever_smk1 + ever_drink4 +
                 cohexp_smoke1 + cohexp_smoke4)
out_smoke <- glm(mod_smoke, data=sub.dat, family = binomial(link="probit"))
summary(out_smoke)

# Smoking behavior degree weighted cohesion
mod_smoke <- formula(ever_smk4 ~ factor(school) + age4 + female4 + lunch4 + grades4 + health4 + rooms4 +
                       parent_smoke4 + sibling_smoke4 + ever_smk1 + ever_drink4 +
                       cohexp_deg_smoke4)
out_smoke <- glm(mod_smoke, data=sub.dat, family = binomial(link="probit"))
summary(out_smoke)

# Drinking behavior
mod_drink <- formula(ever_drink4 ~ factor(school) + age4 + female4 + lunch4 + grades4 + health4 + rooms4 +
                 parent_smoke4 + sibling_smoke4 + ever_drink1 + ever_smk4 +
                 cohexp_drink1 + cohexp_drink4)
mod_drink <- glm(mod_drink, data=sub.dat.d, family = binomial(link="probit"))
summary(mod_drink)

# Drinking behavior with degree weighted exposure
mod_drink <- formula(ever_drink4 ~ factor(school) + age4 + female4 + lunch4 + grades4 + health4 + rooms4 +
                       parent_smoke4 + sibling_smoke4 + ever_drink1 + ever_smk4 +
                       cohexp_deg_drink4)
mod_drink <- glm(mod_drink, data=sub.dat.d, family = binomial(link="probit"))
summary(mod_drink)

