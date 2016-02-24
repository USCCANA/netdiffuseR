rm(list=ls())
fakesurvey <- read.csv("data-raw/fakesurvey.csv")
save(fakesurvey, file = "data/fakesurvey.rdata")
