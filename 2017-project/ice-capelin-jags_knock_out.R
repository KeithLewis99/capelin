#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-03-19, R version 3.3.3 (2017-03-06)   

# The purpose of this file is to test the sensitivity of the results from the models in "ice-capelin-sequential-jag.R".  Specifically, the purpose is to test the sensitivity of any one year. 

## libraries----
library(R2jags)
library(rjags)
library(readr)
library(tidyr)
library(dplyr)
library(psych)
library(ggplot2)
library(lattice)
library(magrittr)
library(loo)

rm(list=ls())

## read in source code-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

df <- read_csv("Bayesian/biomass_cond_ag1/all_data_a1.csv")
df$year

## R/M3 knock-off----
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])"

m.RM3 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   
log_lik[i] <- logdensity.norm(N2[i], mu[i], sigma^-2)

# 3. Discrepancy measures
expY[i] <- mu[i]
varY[i] <- sigma^2
PRes[i] <- (N2[i] - expY[i]) / sigma
PResNew[i] <- (N2_new[i] - expY[i]) / sigma
#Squared residuals
D[i] <- pow(PRes[i], 2) #SSQ
DNew[i] <- pow(PResNew[i], 2)
#CD[i] <- 
}
#Sum of squared Pearson residuals:
Fit <- sum(D[1:N]) # look at overdispersion
FitNew <- sum(DNew[1:N])

# 2. Priors
alpha ~ dnorm(0, 20^-2) 
beta ~ dnorm(0, 10^-2) 
gamma ~ dgamma(5, 1/3) 
delta ~ dgamma(2.8, 1)
epsilon ~ dnorm(0, 20^-2) 
sigma ~ dunif(0, 100) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#delta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
df2 <- df
View(df2)
df3 <- subset(df2, year>2002)
#View(df3)


#df3: 2003-2017
x <- 12 #1:14
insert <- df3[x,]
df3 <- df3[-x, ]

#df2: 1999-2017
y <- x+5
insert_year <- 2015
insert_row  <- x

yaxis1 = "ln_biomass_med" 
ylab1 = "ln(capelin biomass(ktons))"

st <- df2[c("year", "surface_tows_lag2")][19,]
ti <- df2[c("year", "tice")][19,]
co <- df2[c("year", "meanCond_lag")][19,]
co[,2]

num_forecasts =  1# 1 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, st[
                        ,2]),
                   TI=c(df3$tice, ti[,2]), 
                   CO=c(df3$meanCond_lag, co[,2]),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM3))

#UPDATE WITH MORE BURN INS
run_RM3 <-update(run_RM3, n.iter = 300000, n.thin = 50, n.burnin = 100000)


##R/M3-plot 5----
# DIAGNOSTICS
out <- run_RM3$BUGSoutput 


# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon')
MyBUGSChains(out, vars)
#ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm_3/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
#ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm_3/auto_corr.pdf", width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/rm_3/fit_obs.pdf')
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf('Bayesian/rm_3/resid_covar.pdf')
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "ST", ylab = "Pearson resids")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM3$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median')
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))

#generate prediciton intevals using N2_new
y_new = run_RM3$BUGSoutput$sims.list$N2_new
p_med = apply(y_new,2,'median')
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm_3/pi.csv")


## RM_3-Results----
dic_RM3 <- dic.samples(run_RM3$model, n.iter=1000, type="pD")
dic_RM3sum <- sum(dic_RM3$deviance)

log_lik_rm3 = run_RM3$BUGSoutput$sims.list$log_lik
w_rm3 <- waic(log_lik_rm3)
w_rm3 <- loo(log_lik_rm3)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm_3/posteriors.pdf", width=10, height=8, units="in")

#plot credible and prediction intervals
plotCredInt2(df3, insert, 
             yaxis = yaxis1, 
             ylab = ylab1, 
             y_line = y_med, 
             ci = ci_df3, 
             dpp = p_med,
             dpi = pi_df3, 
             insert_year = insert_year,
             insert_row = insert_row)

ggsave(paste0("Bayesian/leave_out/credInt", insert_year, ".png"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/rm_3/params.csv")


# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)
delta <- posterior_fig(out$sims.list$delta)
epsilon <- posterior_fig(out$sims.list$epsilon)

"log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta))"


#alpha
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#beta
priorsd <- 10
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma
priormean <- 5
priorsd <- 1/3
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#delta
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#epsilon
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(epsilon$df)-0.3, max(epsilon$df) + 0.3)
x_label <- "epsilon"
bin_1 <- mean(epsilon$df)/100

p5 <- postPriors(df = epsilon$df, df2 = prior, df3 = epsilon$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, p5, labels = c("A", "B", "C", "D", "E"), ncol=2)

ggsave("Bayesian/rm_3/priorPost.png", width=10, height=8, units="in")



#df3[13, 2:10] <- NA
#df3[14, ] <- NA
df3 <- df3[-15, ]#2017
df3 <- df3[-14, ]#2016
df3 <- df3[-13, ]#2015
df3 <- df3[-12, ]
df3 <- df3[-11, ]
df3 <- df3[-10, ]

df3[8, "ln_biomass_med"] <- NA
df3[8, "ln_biomass_lci"] <- NA
df3[8, "ln_biomass_uci"] <- NA
