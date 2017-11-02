rm(list = ls())

library(EpiModel)

set.seed(12345)

nw <- network::network.initialize(n = 1e3, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "risk", rep(0:1, each=500))

formation <- ~edges + nodefactor("risk") + nodematch("risk") + concurrent

target.stats <- c(250, 375, 225, 100)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 80)
coef.diss

# Step 1
est1 <- netest(nw, formation, target.stats, coef.diss)

# Step 2
dx <- netdx(est1, nsims=10, nsteps = 1000)
dx

# Step 3
init    <- init.net(i.num = 50)
param   <- param.net(inf.prob = 0.1, act.rate=5, rec.rate = .02)
control <- control.net(type = "SIS", nsteps=500, nsims = 10, epi.by = "risk")

sim1    <- netsim(est1, param, init, control)


plot(sim1)
summary(sim1, at = 500)

plot(sim1, y = c("i.num.risk0", "i.num.risk1"), legend=TRUE)

oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1, 2))
plot(sim1, type="network", at=1, col.status = TRUE)
plot(sim1, type="network", at=500, col.status = TRUE)
par(oldpar)
