rm(list=ls())
fakesurvey <- read.csv("data-raw/fakesurvey.csv")
fakesurvey$note <- as.character(fakesurvey$note)
save(fakesurvey, file = "data/fakesurvey.rdata")
