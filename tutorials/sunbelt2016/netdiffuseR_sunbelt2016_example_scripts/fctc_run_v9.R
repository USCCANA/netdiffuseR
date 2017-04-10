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

  # Both the edgelist and the data must span the same time period otherwise you get an error.
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

  # Computing exposure
  fctcdiffnet2[["cohexp"]]  <- exposure(fctcdiffnet2)
  #fctcdiffnet2[["seexp"]]   <- exposure(fctcdiffnet2, alt.graph = "se", valued = TRUE)

  # Calculate degree scores
  fctcdiffnet2[["indegree"]]  <- dgr(fctcdiffnet2, cmode= "indegree")
  fctcdiffnet2[["outdegree"]]  <- dgr(fctcdiffnet2, cmode= "outdegree")

  test <- diffnet.attrs(fctcdiffnet2, as.df = TRUE)
  test$adopted <- as.integer(with(test, toa == per))
  #drop post adoption cases
  test <- subset(test, per <=  toa)
  # Estimate regression for each network separately
  mod_all <- as.formula(paste("adopted ~ factor(per) + cohexp + factor(continent) + population "))
  out_all <- glm(mod_all, data=test, family = binomial(link="logit"))

  print(summary(out_all))  # Print within the loop

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

  # Getting the exposure and degree score data
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

# Merging threshold
threshold6 <- fctc_diffnets[[6]]
threshold6[["thr6"]] <- threshold(threshold6)
threshold6 <- subset(diffnet.attrs(threshold6, as.df=TRUE),
                     select=c(id, per, thr6))
threshold6$year <- threshold6$per
fctc_att <- merge(fctc_att, threshold6, by=c("id", "year"))

# Merging threshold
infection6 <- fctc_diffnets[[6]]
infection6[["inf6"]] <- infection(infection6)
infection6 <- subset(diffnet.attrs(infection6, as.df=TRUE),
                     select=c(id, per, inf6))
infection6$year <- infection6$per
fctc_att <- merge(fctc_att, infection6, by=c("id", "year"))

# Running the model on adoption
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

# Running the model on thresholds
mod_all <- formula(
  thr6 ~ factor(year) + factor(continent) +
  population + democracy + GL_total + tobac_prod +
  degree_01 + degree_02 + degree_03 + degree_12 + degree_13 + degree_14
)
out_all <- glm(mod_all, data=fctc_att, family = gaussian())
summary(out_all)

# Running the model on infection
mod_all <- formula(
  inf6 ~ factor(year) + factor(continent) +
    population + democracy + GL_total + tobac_prod +
    degree_01 + degree_02 + degree_03 + degree_12 + degree_13 + degree_14
)
out_all <- glm(mod_all, data=fctc_att, family = gaussian())
summary(out_all)

# Adding time varying co-variates


# Implementation Regression
imp_mod_all <- as.formula(paste(
  "meanall ~ factor(who_reg2) + toa_year_fctc + perc_male_smoke + labor +
  population + democracy + GL_total + tobac_prod +
  degree_01 + degree_02 + degree_03 + degree_12 + degree_13 + degree_14"
)
)
imp_out_all <- glm(imp_mod_all, data=fctc_att, family = gaussian(link="identity"),
               subset = ((year ==  10) & (degree_14>1)) & (id!=183) & (id!=163))
summary(imp_out_all)

