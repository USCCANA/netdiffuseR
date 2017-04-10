# This R script uses netdifuseR to analyze the diffusion of the
# Framework Convention for Tobacco Control (FCTC) replicating th
# analyses reported in Valente et al (2015) The Global Diffusion
# of Tobacco Control, Social Science & Medicine


rm(list = ls())
library(foreign)
library(netdiffuseR)

# fctc_att    <- read.dta("c:/misc/fctc/data/attributes_v3.dta")
# fctc_net    <- read.dta("c:/misc/fctc/data/allnets.dta")
fctc_att    <- read.dta("playground/20160317_fctc/attributes_v3.dta")
fctc_net    <- read.dta("playground/20160317_fctc/allnets.dta")

# Create lists to be store the diffnet file and the outputs from the model.
fctc_diffnets <- vector("list", 6)   # Because we read 6 networks
fctc_logits   <- vector("list", 6)

relations <- c(1, 2, 3, 12, 13, 14)

# Create loop to cycle through the networks
for (i in 1:length(relations)) {
  sub_fctc_net <- subset(fctc_net, relation==relations[i]) # take a subset

  # [2016-03-17]: Both the edgelist and the data must span the same time period
  #  otherwise you get an error.
  tran <- range(sub_fctc_net$year)

  # Network 1 (distance) is static so create this condition on t1 so it only varies for
  #  the networks other than the one inferred by relation == 1
  t1 <- if (relations[i] %in%  c(1,12:14)) NULL else sub_fctc_net$year

  fctcdiffnet2 <- edgelist_to_diffnet(
    edgelist = sub_fctc_net[,c("id","nom")],
    t0       = sub_fctc_net$year,
    t1       = t1,
    # Be sure the data are the same as the spanned time in the network!
    dat      = subset(fctc_att, year %in% tran[1]:tran[2]),
    idvar = "id", toavar = "toa_year_fctc", timevar="year",
    warn.coercion = FALSE,
    fill.missing = "both"                         # fill missing in either network or adoption
  )
  # [2006-03-27]: There is no 'fctcdiffnet' object in the data.
  # # Make sure it looks good
  # summary(fctcdiffnet)

  # Computing exposure
  fctcdiffnet2[["cohexp"]]  <- exposure(fctcdiffnet2)
  #fctcdiffnet2[["seexp"]]   <- exposure(fctcdiffnet2, alt.graph = "se", valued = TRUE)

  # Calculate degree scores
  fctcdiffnet2[["indegree"]]  <- dgr(fctcdiffnet2, cmode= "indegree")
  fctcdiffnet2[["outdegree"]]  <- dgr(fctcdiffnet2, cmode= "outdegree")

  # [2016-03-27]: Going back to 'test <- diffnet.attrs(' otherwise a bunch of
  #   errors. You were storing the data.frame in the list of diffnets so the code
  #   wont run after this loop. That's why I changed this
  # fctcdiffnet2 <- diffnet.attrs(fctcdiffnet2, as.df = TRUE)
  # fctcdiffnet2$adopted <- as.integer(with(fctcdiffnet2, toa == per))
  # #drop post adoption cases
  # fctcdiffnet2 <- subset(fctcdiffnet2, per <=  toa)
  # # Estimate regression for each network separately
  # mod_all <- as.formula(paste("adopted ~ factor(per) + cohexp + factor(continent) + population "))
  #
  # # [2016-03-27]: There is no 'test' object
  # # out_all <- glm(mod_all, data=test, family = binomial(link="logit"))
  # out_all <- glm(mod_all, data=fctcdiffnet2, family = binomial(link="logit"))
  #
  # print(summary(out_all))
  #
  # #  Storing outputs
  # fctc_diffnets[[i]] <- fctcdiffnet2
  # fctc_logits[[i]]   <- out_all
  test <- diffnet.attrs(fctcdiffnet2, as.df = TRUE)
  test$adopted <- as.integer(with(test, toa == per))
  #drop post adoption cases
  test <- subset(test, per <=  toa)
  # Estimate regression for each network separately
  mod_all <- as.formula(paste("adopted ~ factor(per) + cohexp + factor(continent) + population "))

  # [2016-03-27]: There is no 'test' object
  # out_all <- glm(mod_all, data=test, family = binomial(link="logit"))
  out_all <- glm(mod_all, data=test, family = binomial(link="logit"))

  print(summary(out_all))

  #  Storing outputs
  fctc_diffnets[[i]] <- fctcdiffnet2
  fctc_logits[[i]]   <- out_all
}

# Plotting
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,3))
for (i in 1:length(fctc_diffnets))
  plot_threshold(fctc_diffnets[[i]], main=paste("TOA and Threshold for relation",i))

par(oldpar)

# Merging exposure data
for (i in 1:length(fctc_diffnets)) {

  # Getting the data
  exposurei <- diffnet.attrs(fctc_diffnets[[i]], as.df = TRUE)
  exposurei <- subset(exposurei, select=c(id, per, cohexp))

  degreei <- diffnet.attrs(fctc_diffnets[[i]], as.df = TRUE)
  degreei <- subset(degreei, select=c(id, per, indegree))

  # Renaming the last column (the cohesive)
  colnames(exposurei)[3] <- sprintf("cohexp_%02d", relations[i])
  colnames(degreei)[3] <- sprintf("degree_%02d", relations[i])

  # Merging
  if (i == 1) exposures <- exposurei
  else exposures <- merge(exposures, exposurei, by=c("id", "per"))
  if (i == 1) degrees <- degreei
  else degrees <- merge(degrees, degreei, by=c("id", "per"))
}

# Changing names of the id and per columns to merge it with the data
colnames(exposures)[1:2] <- c("id", "year")
colnames(degrees)[1:2] <- c("id", "year")

# Merging with the original data
fctc_att <- merge(fctc_att, exposures, by=c("id", "year"))
fctc_att <- merge(fctc_att, degrees, by=c("id", "year"))

# Running the model
mod_all <- as.formula(paste(
  # Here we generate the variable 'adopted' in the formula
  "I(as.integer(toa_year_fctc==year)) ~ factor(year) + factor(continent) +
  population + democracy + GL_total + tobac_prod +
  degree_01 + degree_02 + degree_03 + degree_12 + degree_13 + degree_14 +
  cohexp_01 + cohexp_02 + cohexp_03 + cohexp_12 + cohexp_13 + cohexp_14 "
  )
)

out_all <- glm(mod_all, data=fctc_att, family = binomial(link="logit"),
               subset = year <=  toa_year_fctc)
summary(out_all)

# Regressions on Thresholds, infection, & susceptibility

# Adding time varying co-variates


# Implementation Regression
imp_mod_all <- as.formula(paste(
  "meanall ~ factor(continent) + toa_year_fctc +
  population + democracy + GL_total + tobac_prod +
  degree_01 + degree_02 + degree_03 + degree_12 + degree_13 + degree_14 +
  cohexp_01 + cohexp_02 + cohexp_03 + cohexp_12 + cohexp_13 + cohexp_14 "
)
)

# [2016-03-27]: Logit is only for dichotomous response
# imp_out_all <- glm(imp_mod_all, data=fctc_att, family = binomial(link="logit"),
#                subset = year ==  10)

# This is more appropiate
imp_out_all <- glm(imp_mod_all, data=fctc_att, family = gaussian(link="identity"),
               subset = year ==  10)

# Which is equivalent to a OLS
summary(lm(imp_mod_all, data=fctc_att, subset = year ==  10))
summary(imp_out_all)
