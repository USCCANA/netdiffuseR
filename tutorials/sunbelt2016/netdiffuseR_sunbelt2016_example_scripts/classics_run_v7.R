# This R script Demonstrates how to USE NetdifuseR to analyze the
# 3 classic network diffusion datasets and estimate cohesion vs
# Structural equivalence influences on adoption

rm(list = ls())
library(foreign)
library(netdiffuseR)

# Read MI data note the "groupvar" option
# mi_att <- read.dta("c:/misc/diffnet/mi_v2.dta")
mi_att <- read.dta("playground/20160305_tom_classic_results/mi_v2.dta")
midiffnet <- survey_to_diffnet(mi_att, idvar="id", netvars=c("net11", "net12", "net13",
                                                             "net21", "net22", "net23"),
                               toavar="toa", groupvar = "city", warn.coercion = FALSE)
# Make sure it makes sense
summary(midiffnet)

# Calculate exposures
midiffnet[["cohexp"]] <- exposure(midiffnet)
midiffnet[["seexp2"]] <- exposure(midiffnet, alt.graph = "se", groupvar = "city")

midiffnet.df <- diffnet.attrs(midiffnet, as.df = TRUE)                  # Convert to dataframe
midiffnet.df$adopted <- as.integer(with(midiffnet.df, ado == per))          # Set adoption variable
midiffnet.df <- midiffnet.df[midiffnet.df$per <=  midiffnet.df$toa, ]        # Keep pre-adoption time only
mod_all <- as.formula(paste("adopted ~ factor(per) + proage + journ2 + science + detail + cohexp + seexp2  "))
out_all <- glm(mod_all, data=midiffnet.df, family = binomial(link="logit"))
summary(out_all)
# Draw a cumulate and new adopters graph
plot_adopters(midiffnet)
# Plot the the diffusion process
plot_diffnet(midiffnet)
plot_diffnet(midiffnet, slices=c(1, 6, 12, 18))



# BF data
bf_att <- read.dta("playground/20160305_tom_classic_results/brfarmers.dta")
bfdiffnet <- survey_to_diffnet(bf_att, idvar="id", netvars=c("net11", "net12", "net13",
                                                             "net21", "net22", "net23",
                                                             "net31", "net32", "net33"),
                               toavar="toa", groupvar = "village")

summary(bfdiffnet)

bfdiffnet[["cohexp"]] <- exposure(bfdiffnet)
bfdiffnet[["seexp"]] <-  exposure(bfdiffnet, alt.graph="se", groupvar="village",
                                  valued = TRUE)
# Store village variable with diffnet object
bfdiffnet[["village"]] <- bf_att$village

bfdiffnet.df <- diffnet.attrs(bfdiffnet, as.df = TRUE)
bfdiffnet.df$adopted <- as.integer(with(bfdiffnet.df, ado == per))
bfdiffnet.df <- bfdiffnet.df[bfdiffnet.df$per <=  bfdiffnet.df$toa, ]
mod_all <- as.formula(paste("adopted ~ factor(per) + visits + news1 + immexp + cohexp + seexp  "))
out_all <- glm(mod_all, data=bfdiffnet.df, family = binomial(link="logit"))
summary(out_all)
# Draw a cumulate and new adopters graph
plot_adopters(bfdiffnet)
# Plot the the diffusion process but just one time point
plot_diffnet(bfdiffnet, slices=10)
# Plot the the diffusion process but just one time point and one village
#plot_diffnet((with(bfdiffnet, village==10)), slices=10)

# Plot the the diffusion process for just one village
bfdiffnet10<-bfdiffnet[["village"]]==10
plot_diffnet(bfdiffnet[bfdiffnet10], slices=10)

# KFP Data
kfp_att<- read.dta("playground/20160305_tom_classic_results/kfp_v3.dta")
kfpdiffnet <- survey_to_diffnet(kfp_att, idvar="id", netvars=c("net11", "net12", "net13", "net14", "net15",
                                                 "net21", "net22", "net23", "net24", "net25",
                                                 "net31", "net32", "net33", "net34", "net35"),
                  toavar="toa", groupvar = "village")

summary(kfpdiffnet)

kfpdiffnet[["cohexp"]] <- exposure(kfpdiffnet)
kfpdiffnet[["seexp"]] <-  exposure(kfpdiffnet, alt.graph="se", groupvar="village",
                                   valued=TRUE)


kfpdiffnet.df <- diffnet.attrs(kfpdiffnet, as.df = TRUE)
kfpdiffnet.df$adopted <- as.integer(with(kfpdiffnet.df, ado == per))
kfpdiffnet.df <- kfpdiffnet.df[kfpdiffnet.df$per <=  kfpdiffnet.df$toa, ]
mod_all <- as.formula(paste("adopted ~ factor(per) + sons + pregs + cohexp + seexp  "))
out_all <- glm(mod_all, data=kfpdiffnet.df, family = binomial(link="logit"))
summary(out_all)

########################################################################
#                       The End                                        #
########################################################################

