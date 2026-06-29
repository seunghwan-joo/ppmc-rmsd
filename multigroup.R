setwd("~/Dropbox/Research/RMSD Item Fit/PPMC RMSD")
library(mirt)
library(ggplot2)
source("Functions/twopl.gen.R")
source("Functions/mg.twopl.mcmc.fast.R")
source("Functions/rmsd.ppmc.fast3.R")

# Data generation 
n.samples <- 500
n.items <- 10
n.groups <- 2
dif.item <- 10
dif.group <- 2
b.dif.size <- 1
a.dif.size <- 0
mn.diff <- 0

# Generating item parameters
a <- runif(n.items,0.5,2)
b <- runif(n.items,-2,2)
apars <- matrix(NA, nr=n.items, nc=n.groups)
bpars <- matrix(NA, nr=n.items, nc=n.groups)
for(g in 1:n.groups){
  apars[,g] <- a
  bpars[,g] <- b
}

# Generating DIF
apars[dif.item,dif.group] <- apars[dif.item,dif.group] + a.dif.size
bpars[dif.item,dif.group] <- bpars[dif.item,dif.group] + b.dif.size

# Generating group-specific theta
mn <- c(0,rep(mn.diff, n.groups-1))
theta <- matrix(NA, nr=n.samples, nc=n.groups)
for(g in 1:n.groups){
  theta[,g] <- rnorm(n.samples,mn[g],1)
}

# Generating response data and group variable
for(g in 1:n.groups){
  tmp <- twopl.gen(a=apars[,g], b=bpars[,g], theta=theta[,g])
  if(g==1) dat <- tmp
  if(g>1) dat <- rbind(dat, tmp)
}
groups <- rep(1:n.groups, each=n.samples)

# MCMC estimation
chains <- 2 # 2 chains
burn <- 500 # 500 burn ins
it <- 1000 # 1000 iterations
fit <- mg.twopl.mcmc(data=dat, group=groups, it=it, burn=burn, chains=chains)
mcmc.pars <- fit$item.par
mcmc.scores <- fit$theta.est
mcmc.mn <- fit$mn.est
mcmc.std <- fit$SD.est
fit$item.accept
fit$item.conv
fit$mn.est
fit$SD.est
fit$model.fit

# mirt estimation
itemnames <- colnames(dat)
grp <- as.character(groups)
mgfit <- multipleGroup(dat, 
                       1, 
                       group = grp, 
                       itemtype = "2PL",
                       invariance = c(itemnames, 'free_means', 'free_var'))
results <- coef(mgfit, simplify=TRUE, IRTpars=TRUE)
mirt.pars <- results$`1`$items[,1:2]
mirt.scores <- fscores(mgfit)
mirt.rmsd <- RMSD_DIF(mgfit)
mirt.rmsd
# Compare results
# Items
genpars <- cbind(a,b)
cbind(genpars,mcmc.pars,mirt.pars)
genpars - mcmc.pars
genpars - mirt.pars
# Thetas
mcmc.cors <- rep(NA,n.groups)
mirt.cors <- rep(NA,n.groups)
for(g in 1:n.groups){
  mcmc.cors[g] <- cor(mcmc.scores[groups==g],theta[,g])
  mirt.cors[g] <- cor(mirt.scores[groups==g],theta[,g])
}
mcmc.cors
mirt.cors

# RMSD with PPMC
# Read in MCMC iteration samples
a.sam <- matrix(NA,n.items,(it-burn))
b.sam <- matrix(NA,n.items,(it-burn))
th.sam <- matrix(NA,n.samples,(it-burn))
for(c in 1:chains){
  if(c==1){ 
    a.sam <- fit$a.sample[,,c]
    b.sam <- fit$b.sample[,,c]
    th.sam <- fit$th.sample[,,c]
  }
  if(c>1){
    a.sam <- cbind(a.sam, fit$a.sample[,,c])
    b.sam <- cbind(b.sam, fit$b.sample[,,c])
    th.sam <- cbind(th.sam, fit$th.sample[,,c])
  }
}

# Compute RMSD for each group
rmsd.ppp <- matrix(NA, nr=n.items, nc=n.groups)
rmsd.obs <- array(NA, dim=c(chains*(it-burn), n.items, n.groups))
rmsd.rep <- array(NA, dim=c(chains*(it-burn), n.items, n.groups))
for(g in 1:n.groups){
  ppmc <- rmsd.ppmc(data=dat[groups==g,], mn=mcmc.mn[g], std=mcmc.std[g],
                    a.sample=a.sam, 
                    b.sample=b.sam, 
                    th.sample=th.sam[groups==g,])
  rmsd.ppp[,g] <- ppmc$ppp
  rmsd.obs[,,g] <- ppmc$obs.rmsd
  rmsd.rep[,,g] <- ppmc$rep.rmsd
}
rmsd.ppp

# Figure 1
# Histogram for nonDIF
obs <- rmsd.obs[,1,dif.group]
rep <- rmsd.rep[,1,dif.group]
df <- data.frame(score=c(obs,rep),
                 RMSD=c(rep("observed",length(obs)),rep("replicated",length(rep))))
p1 <- ggplot(df, aes(x = score, fill = RMSD)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(x = "RMSD", y = "Frequency",title = "(a) Distributions of RMSDs (nonDIF)") +
  theme_minimal()

# Scatter plot for nonDIF
df <- data.frame(x=obs,y=rep)
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0, .15), ylim = c(0, .15)) +
  labs(
    title = "(b) Scatterplot of Observed and Replicated RMSDs",
    x = "Observed RMSD",
    y = "Replicated RMSD"
  ) +
  theme_minimal()

# Histogram for DIF
obs <- rmsd.obs[,dif.item,dif.group]
rep <- rmsd.rep[,dif.item,dif.group]
df <- data.frame(score=c(obs,rep),
                 RMSD=c(rep("observed",length(obs)),rep("replicated",length(rep))))
p3 <- ggplot(df, aes(x = score, fill = RMSD)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(x = "RMSD", y = "Frequency",title = "(c) Distributions of RMSDs (DIF)") +
  theme_minimal()

# Scatter plot for DIF
df <- data.frame(x=obs,y=rep)
p4 <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0, .15), ylim = c(0, .15)) +
  labs(
    title = "(d) Scatterplot of Observed and Replicated RMSDs",
    x = "Observed RMSD",
    y = "Replicated RMSD"
  ) +
  theme_minimal()

library(patchwork)
# Arrange in 2x2 grid
(p1 | p2) / (p3 | p4)