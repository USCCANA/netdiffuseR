rm(list=ls())

fakesurvey <- read.csv("data-raw/fakesurvey.csv")
fakesurvey$note <- as.character(fakesurvey$note)

dat <- rbind(fakesurvey)
dat$toa <- sample(1990:1991, nrow(dat), TRUE)
dat <- rbind(dat, dat)

dat$time <- c(
  rep(1990, nrow(fakesurvey)),
  rep(1991, nrow(fakesurvey))
  )

# Aging people
dat$age[dat$time==1991] <-
  dat$age[dat$time==1991] + 1

# Changing some links (isolated wont be such anymore!)
dat$net1[with(dat, which(id==2 & group == 2 & time == 1991))] <- 1
dat$note[with(dat, which(id==2 & group == 2 & time == 1991))] <- "Now is not isolated!"
dat$note[dat$time==1990] <- paste("First wave:",dat$note[dat$time==1990])
dat$note[dat$time==1991] <- paste("First wave:",dat$note[dat$time==1991])

fakesurveyDyn <- dat

save(fakesurveyDyn, file="data/fakesurveyDyn.rdata")

# Playing with the data
library(netdiffuseR)

idvar <- "id"
netvars <- paste0("net", 1:3)
groupvar <- "group"

x <- survey_to_diffnet(
  dat, "id", c("net1", "net2", "net3"), "toa", "group" ,
  timevar = "time", keep.isolates = TRUE)

plot_diffnet(x, vertex.cex = 1.5, displaylabels = TRUE)
