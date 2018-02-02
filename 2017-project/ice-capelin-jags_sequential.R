#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-01-15, R version 3.3.3 (2017-03-06)             #

# The ice-capelin-covariate file reveled some problems with the approach to date.  We had hoped to create "lumpy vaults" that would help us predict capelin abundance.  However the we merely created tilted vaults.  This does little to improve our understanding over Ale's original tice model.  Paul then suggested a sequential approach in a Bayesian framework.  Model recuitment, convert to spawning biomass, and then predict adult mortality.
# The purpose of this file is to:
# 1) build these sequential models in the Bayesian framework using JAGS
# 2) estimate the number of immature capelin in teh age 2 class
# However it is not clear if this will work conceptually because of the timing and limitations of the capelin surveys.

#So now, we are simply using JAGS as a convenient way to calculate confidence/credible intervals and prediction intervals

# The purpose of this file is to implement

## libraries----
library(R2jags)
library(rjags)
library(readr)
library(tidyr)
library(dplyr)
library(psych)
library(ggplot2)
library(lattice)

rm(list=ls())

## read in source code-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")

## load data----
## capelin (1985-2017)----
# source: Age disaggregate abundance for Keith Lewis - 2017 added_v1.xlsx
cap <- read_csv('data/capelin-2017.csv')
glimpse(cap)
#View(cap)
cap$ln_abun_med <- log(cap$abundance_med)
cap$ln_ab_lci <- log(cap$ab_lci)
cap$ln_ab_uci <- log(cap$ab_uci)
cap$ln_biomass_med <- log(cap$biomass_med)
cap$ln_bm_lci <- log(cap$bm_lci)
cap$ln_bm_uci <- log(cap$bm_uci)

## ice (1969-2017)----
#source: ice-chart-processing-data-v3.R
ice <- read_csv('output-processing/capelin-m1.csv')
glimpse(ice)

## condition (1995-2015)-----
# source: capelin-condition.R
#cond <- read_csv('data/condition_ag1_out.csv')
glimpse(cond)
#View(cond)

#add years
cond[22:23,4] <- NA
cond[c(22, 23), 1] <- matrix(c(2016, 2017), ncol = 1) 

#lag data
cond$meanCond_lag <- lag(cond$meanCond, 1)
cond$medCond_lag <- lag(cond$medCond, 1)

# fill in for 2017 - remove line when new data are available
cond$meanCond_lag[23] <- cond$meanCond_lag[22]

#Fran's original
condResids <- read_csv('data/archive/condition_resids.csv')
glimpse(condResids)
condResids[19:20,2] <- NA
condResids[c(19, 20), 1] <- matrix(c(2016, 2017), ncol = 1) 
condResids$resids_lag <- lag(condResids$resids, 1)
condResids$resids_lag[20] <- condResids$resids_lag[19]
View(condResids)

## larval data (1990-2017)----
# source "capelin_age_disaggregate_abundance.xlsx":Larval indices
larvae <- read_csv('data/capelin_larval_indices.csv')
glimpse(larvae)
larvae$surface_tows_lag2 <- lag(larvae$surface_tows, 2)

## Pseudocalanus data (1999-2016)----
# source "Copy of Copy of PSEUSP27_1999_2016.xlsx"
pscal <- read_csv("data/pseudocal_1999_2016.csv")
glimpse(pscal)
pscal$year

# filter out year and use just the "total" stage.  Calcuate mean and sd by year
ps_tot <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage == "c6") %>%
     group_by(year) %>%
     summarise(ps_meanTot = mean(density), ps_sdTot = sd(density))

#add years
ps_tot[19,] <- NA
ps_tot[19, 1] <- matrix(c(2017), ncol = 1) 

# lag mean and sd and /1000 to scale to tice
ps_tot$ps_meanTot_lag2 <- lag(ps_tot$ps_meanTot, 2)
ps_tot$ps_sdTot_lag2 <- lag(ps_tot$ps_sdTot, 2)
#View(ps_tot)

##join the above data sets----
capelin_join <- left_join(cap, ice, by = "year")
#capelin_join <- left_join(capelin_join, cond, by = "year")
capelin_join <- left_join(capelin_join, condResids, by = "year")
capelin_join <- left_join(capelin_join, larvae, by = "year")
capelin_join <- left_join(capelin_join, ps_tot, by = "year")

df <- capelin_join

glimpse(df)
#df[c("biomass_med", "ln_biomass_med", "tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")]
df[c("biomass_med", "ln_biomass_med", "tice", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")]

# normalize the data set except for year----
#df1 <- subset(df, year>1995)
df1 <- subset(df, year>1998)
glimpse(df1)
View(df1)
#cols <- c("tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")
cols <- c("tice", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")
df1 %<>%
     mutate_each_(funs(scale),cols)
glimpse(df1)

# relationships and correlations among RV and EV
pairs.panels(df1[c("ln_biomass_med", "tice", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses,
             cex.labels = 1.5,
             cex.cor = 2
)


## not sure of use----
df_age <- df[c("year", "logcapelin", "age2", "age3", "age4", "age5", "age6", "age2PerMat", "resids_adj")]

mean(df_age$age2PerMat, na.rm = T)
summary(df_age$age2PerMat)
plot(df_age$logcapelin, df_age$age2PerMat)
plot(df_age$year, df_age$age2PerMat)
summary(lm(age2PerMat ~ logcapelin, data = df_age))
m1 <- lm(age2PerMat ~ logcapelin, data = df_age)

layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(m1)
layout(matrix(1,1))

#test to see if lag occurs
plot(lag(df_age$logcapelin, 1), df_age$age2PerMat)
summary(lm(age2PerMat ~ lead(logcapelin, 1), data = df_age))
# no lag occurs

Age2PMat <- df_age$age2PerMat
df_age <- df_age %>%
     mutate(N3_6 = age3 + age4 + age5 + age6)
N3_6 <- df_age$N3_6

ST <- df$Ssurface_tows_lag2
TI <- df$tice
N <- nrow(df)

## test 1 first try----
jags_code = '
model {
# Likelihood
for (i in 1:N) {
     #recruitment
     mu_n2[i] <- beta + alpha*ST[i]
     N2[i] ~ dnorm(mu_n2[i], sigma^-2)
     # proportion N2
     Nsp[i] <- N2[i]*Age2Mat + N3_6
     # mortality
     mu_sp[i] <- gamma*TI[i]*(1-(TI[i]/delta))
     Nsp[i] ~ dnorm(mu_sp[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnamma(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
delta ~ dnorm(0, 100^-2) 
sigma_r ~ dunif(0, 100)
sigma_m ~ dunif(0, 100)
}'

num_forecasts = 2 # 10 extra years
model_data <- list(L = c(df$logcapelin, rep(NA, num_forecasts)), 
                   ST=c(df$Ssurface_tows_lag2, c(98, 41)), #from capelin_larval_indices 
                    RA=c(df$resids_adj, c(15, 15)), #made up - need new data
                   N = nrow(df) + num_forecasts)

jags_run <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma'),
                 model.file = textConnection(jags_code))

y_pred = jags_run$BUGSoutput$sims.list$mu
y_med = apply(y_pred,2,'median')
plot(c(test$year,2017:2018),y_med,type='l', ylim=c(0,9))
points(test$year, test$logcapelin)
points(test$year, y_med[1:14])
points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
lines(c(test$year,2017:2018), ci_test[1, ], lty = 3)
lines(c(test$year,2017:2018), ci_test[2, ], lty = 3)


str(test)
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] #these are the extra 2 years
ci_test <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_test[1, ]


## test 2 ----
# not sequential
#"Beta + Alpha*ST+ Gamma*RA",
jags_code = '
model {
# Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- beta + alpha*ST[i] + gamma*RA[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100)
}'


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df1$logcapelin, rep(NA, num_forecasts)), 
                   ST=c(df1$Ssurface_tows_lag2, c(98, 41)), #from capelin_larval_indices 
                   RA=c(df1$resids_adj, c(15, 15)), #made up - need new data
                   N = nrow(df1) + num_forecasts)

jags_run <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma', 'N2_new'),
                 model.file = textConnection(jags_code))

y_pred = jags_run$BUGSoutput$sims.list$mu
y_med = apply(y_pred,2,'median')

str(df1)
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] #these are the extra 2 years
ci_df1 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df1[1, ]

y_new = jags_run$BUGSoutput$sims.list$N2_new
pi_df1 <- apply(y_new,2,'quantile', c(0.05, 0.95))

plot(c(df1$year,2017:2018),y_med,type='l', ylim=c(0,9))
points(df1$year, df1$logcapelin)
#points(df1$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df1$year,2017:2018), ci_df1[1, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[1, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), ci_df1[2, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[2, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), pi_df1[1, ], lty = 2)
lines(c(2017:2018), pi_df1[1, 15:16], lty = 2)
#lines(c(df1$year,2017:2018), pi_df1[2, ], lty = 2)
lines(c(2017:2018), pi_df1[2, 15:16], lty = 2)

polygon(x = c(2017, 2017, 2018, 2018), y = c(pi_df1[1, 15], pi_df1[2,15], pi_df1[2, 16], pi_df1[1, 16]), col="grey75")
text(2007,8, "log(caplein) = Beta + Alpha*ST+ Gamma*RA")



## test 3 ----
# not sequential
#"Beta + Alpha*ST+ Gamma*PS",
jags_code = '
model {
# Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- beta + alpha*ST[i] + gamma*PS[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100)
}'


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df1$logcapelin, rep(NA, num_forecasts)), 
                   ST=c(df1$Ssurface_tows_lag2, c(98, 41)), #from capelin_larval_indices 
                   PS=c(df1$ps_meanTot_lag2, c(80, 80)), #made up - need new data
                   N = nrow(df1) + num_forecasts)

jags_run <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma', 'N2_new'),
                 model.file = textConnection(jags_code))

y_pred = jags_run$BUGSoutput$sims.list$mu
y_med = apply(y_pred,2,'median')

apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] #these are the extra 2 years
ci_df1 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df1[1, ]

y_new = jags_run$BUGSoutput$sims.list$N2_new
pi_df1 <- apply(y_new,2,'quantile', c(0.05, 0.95))

plot(c(df1$year,2017:2018),y_med,type='l', ylim=c(0,9))
points(df1$year, df1$logcapelin)
#points(df1$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df1$year,2017:2018), ci_df1[1, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[1, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), ci_df1[2, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[2, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), pi_df1[1, ], lty = 2)
lines(c(2017:2018), pi_df1[1, 15:16], lty = 2)
#lines(c(df1$year,2017:2018), pi_df1[2, ], lty = 2)
lines(c(2017:2018), pi_df1[2, 15:16], lty = 2)

polygon(x = c(2017, 2017, 2018, 2018), y = c(pi_df1[1, 15], pi_df1[2,15], pi_df1[2, 16], pi_df1[1, 16]), col="grey75")
text(2007,8, "log(caplein) = Beta + Alpha*ST+ Gamma*PS")


## test 4 ----
# not sequential
#"Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj"
jags_code = '
model {
# Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha*TI[i]*(1-(TI[i])*beta) + gamma*RA[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100)
}'


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df1$logcapelin, rep(NA, num_forecasts)), 
                   TI=c(df1$tice, c(70, 70)), #from capelin_larval_indices 
                   RA=c(df1$resids_adj, c(15, 15)), #made up - need new data
                   N = nrow(df1) + num_forecasts)

jags_run <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma', 'N2_new'),
                 model.file = textConnection(jags_code))

y_pred = jags_run$BUGSoutput$sims.list$mu
y_med = apply(y_pred,2,'median')

apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] #these are the extra 2 years
ci_df1 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df1[1, ]

y_new = jags_run$BUGSoutput$sims.list$N2_new
pi_df1 <- apply(y_new,2,'quantile', c(0.05, 0.95))

plot(c(df1$year,2017:2018),y_med,type='l', ylim=c(0,9))
points(df1$year, df1$logcapelin)
#points(df1$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df1$year,2017:2018), ci_df1[1, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[1, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), ci_df1[2, ], lty = 3)
lines(c(df1$year[1:14]), ci_df1[2, 1:14], lty = 3)
#lines(c(df1$year,2017:2018), pi_df1[1, ], lty = 2)
lines(c(2017:2018), pi_df1[1, 15:16], lty = 2)
#lines(c(df1$year,2017:2018), pi_df1[2, ], lty = 2)
lines(c(2017:2018), pi_df1[2, 15:16], lty = 2)

polygon(x = c(2017, 2017, 2018, 2018), y = c(pi_df1[1, 15], pi_df1[2,15], pi_df1[2, 16], pi_df1[1, 16]), col="grey75")

text(2007,8, "log(caplein) = Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj")

## Mortality----
# not sequential
#"delta + Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj"

m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- delta + alpha*TI[i]*(1-TI[i]/beta) + gamma*CO[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dgamma(2, 1/3) 
beta ~ dgamma(8.1, 1/11.11) 
gamma ~ dnorm(0, 100^-2) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
df2 <- df1
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(0, 0)), #from capelin_larval_indices 
                   CO=c(df2$resids_lag, c(0, 0)), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mortality <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                 model.file = textConnection(m.mortality))

#UPDATE WITH MORE BURN INS

## M-Plot 5----
# DIAGNOSTICS
print(run_mortality, intervals=c(0.025, 0.975), digits = 3)
plot(run_mortality)

out <- run_mortality$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_mortality$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma', 'delta')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mortality$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_mortality$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df2 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df2[1, ]

#generate prediciton intevals using N2_new
y_new = run_mortality$BUGSoutput$sims.list$N2_new
pi_df2 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df2$year,2018:2019),y_med,type='l', ylim = c(2,8))
points(df2$year, df2$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df2$year[1:22]), ci_df2[1, 1:22], lty = 3)
#lines(c(df2$year,2017:2018), ci_df2[2, ], lty = 3)
lines(c(df2$year[1:22]), ci_df2[2, 1:22], lty = 3)
#lines(c(df2$year,2017:2018), pi_df2[1, ], lty = 2)
lines(c(2018:2019), pi_df2[1, 23:24], lty = 2)
#lines(c(df2$year,2017:2018), pi_df2[2, ], lty = 2)
lines(c(2018:2019), pi_df2[2, 23:24], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df2[1, 23], pi_df2[2,23], pi_df2[2, 24], pi_df2[1, 24]), col="grey75")

lines(c(df2$year), ci_df2[1, 1:19], lty = 3)
lines(c(df2$year), ci_df2[2, 1:19], lty = 3)
lines(c(2018:2019), pi_df2[1, 20:21], lty = 2)
lines(c(2018:2019), pi_df2[2, 20:21], lty = 2)
polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df2[1, 20], pi_df2[2,20], pi_df2[2, 21], pi_df2[1, 21]), col="grey75")


#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(resids)")
#plot credible and prediction intervals
# better plots

p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df2$year, 2018:2019), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df2$year[1:22], 
                         ymax = ci_df2[2, 1:22], 
                         ymin = ci_df2[1, 1:22]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2018:2019), 
                         ymax = pi_df2[2, 23:24], 
                         ymin = pi_df2[1, 23:24]),
                     alpha = 0.3)        
p


# plot posterior against expected distribution
alpha_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dnorm(0, 100^-2))
str(alpha_post)
dist <- as.data.frame(rnorm(10000, 0, 10))
head(dist)
names(dist)[names(dist) == "rnorm(10000, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

alpha_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dgamma(2, 1/3))
str(alpha_post)
dist <- as.data.frame(rgamma(10000, 2, 1/3))
head(dist)
names(dist)[names(dist) == "rgamma(10000, 2, 1/3)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

beta_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$beta)
beta_prior <- as.data.frame(dgamma(8.1, 1/11.11))
str(beta_post)
dist <- as.data.frame(rgamma(10000, 8, 1/11.11))
head(dist)
names(dist)[names(dist) == "rgamma(10000, 8, 1/11.11)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = beta_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

gamma_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$gamma)
gamma_prior <- as.data.frame(dnorm(0, 100^-2))
dist <- as.data.frame(rnorm(10000, 0, 10))
head(dist)
names(dist)[names(dist) == "rnorm(10000, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = gamma_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

delta_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$delta)
delta_prior <- as.data.frame(dnorm(0, 100^-2))
dist <- as.data.frame(rnorm(10000, 0, 10))
head(dist)
names(dist)[names(dist) == "rnorm(10000, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = delta_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

## Recruitment model----
# not sequential
#"delta + Alpha*surface_tows_lag2) + Gamma*Pseudocal_lag2"

m.recuit = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*PS[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(0, 0)), #from capelin_larval_indices 
                   PS=c(df3$ps_meanTot_lag2, c(0, 0)), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_recruit <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                      model.file = textConnection(m.recruit))

#UPDATE WITH MORE BURN INS

##R-plot 5----
# DIAGNOSTICS
print(run_recruit, intervals=c(0.025, 0.975), digits = 3)
plot(run_recruit)

out <- run_recruit$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_recruit$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$PS, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_recruit$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_recruit$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df3[1, ]

#generate prediciton intevals using N2_new
y_new = run_recruit$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df3$year,2018:2019),y_med,type='l', ylim = c(2,9))
points(df3$year, df3$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[1, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), ci_df3[2, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[2, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), pi_df3[1, ], lty = 2)
lines(c(2018:2019), pi_df3[1, 16:17], lty = 2)
#lines(c(df3$year,2017:2018), pi_df3[2, ], lty = 2)
lines(c(2018:2019), pi_df3[2, 16:17], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df3[1, 16], pi_df3[2,16], pi_df3[2, 17], pi_df3[1, 17]), col="grey75")


text(2008,8, "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*pseudocal_lag2")

#plot credible and prediction intervals
# better plots

p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df2$year, 2018:2019), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df2$year[1:22], 
                         ymax = ci_df2[2, 1:22], 
                         ymin = ci_df2[1, 1:22]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2018:2019), 
                         ymax = pi_df2[2, 23:24], 
                         ymin = pi_df2[1, 23:24]),
                     alpha = 0.3)        
p


# plot posterior against expected distribution
alpha_post <- as.data.frame(run_recruit$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dnorm(0, 100^-2))
str(alpha_post)
dist <- as.data.frame(rnorm(10000000, 0, 10))
head(dist)
names(dist)[names(dist) == "rnorm(1e+07, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
pnorm(1, mean = 0, sd = 100)
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p


## R/M1 model----
# not sequential
#"delta + Alpha*surface_tows_lag2) + Gamma*Pseudocal_lag2"

m.RM1 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*CO[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
#df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(0, 0)), #from capelin_larval_indices 
                   CO=c(df3$resids_lag, c(0, 0)), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_RM1 <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                      model.file = textConnection(m.RM1))

#UPDATE WITH MORE BURN INS

##R/m-plot 5----
# DIAGNOSTICS
print(run_RM1, intervals=c(0.025, 0.975), digits = 3)
plot(run_RM1)

out <- run_RM1$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_RM1$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM1$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_RM1$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df3[1, ]

#generate prediciton intevals using N2_new
y_new = run_RM1$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df3$year,2018:2019),y_med,type='l', ylim = c(2,9))
points(df3$year, df3$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[1, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), ci_df3[2, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[2, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), pi_df3[1, ], lty = 2)
lines(c(2018:2019), pi_df3[1, 16:17], lty = 2)
#lines(c(df3$year,2017:2018), pi_df3[2, ], lty = 2)
lines(c(2018:2019), pi_df3[2, 16:17], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df3[1, 16], pi_df3[2,16], pi_df3[2, 17], pi_df3[1, 17]), col="grey75")

text(2008,8, "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*resids_adj")

#plot credible and prediction intervals
# better plots

p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df2$year, 2018:2019), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df2$year[1:22], 
                         ymax = ci_df2[2, 1:22], 
                         ymin = ci_df2[1, 1:22]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2018:2019), 
                         ymax = pi_df2[2, 23:24], 
                         ymin = pi_df2[1, 23:24]),
                     alpha = 0.3)        
p


# plot posterior against expected distribution
alpha_post <- as.data.frame(run_RM1$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dnorm(0, 100^-2))
str(alpha_post)
dist <- as.data.frame(rnorm(10000000, 0, 10))
names(dist)[names(dist) == "rnorm(1e+07, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p


## R/M2 model----
# not sequential
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta)"

m.RM2 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta)
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2) 
gamma ~ dgamma(2, 1/3) 
delta ~ dgamma(8.1, 1/11.11) 

sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
#df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(0, 0)), #from capelin_larval_indices 
                   TI=c(df3$tice, c(0, 0)), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_RM2 <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                      model.file = textConnection(m.RM2))

#UPDATE WITH MORE BURN INS

##R/M2-plot 5----
# DIAGNOSTICS
print(run_RM2, intervals=c(0.025, 0.975), digits = 3)
plot(run_RM2)

out <- run_RM2$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_RM2$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$TI, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM2$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_RM2$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df3[1, ]

#generate prediciton intevals using N2_new
y_new = run_RM2$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df3$year,2018:2019),y_med,type='l', ylim = c(2,9))
points(df3$year, df3$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[1, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), ci_df3[2, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[2, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), pi_df3[1, ], lty = 2)
lines(c(2018:2019), pi_df3[1, 16:17], lty = 2)
#lines(c(df3$year,2017:2018), pi_df3[2, ], lty = 2)
lines(c(2018:2019), pi_df3[2, 16:17], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df3[1, 16], pi_df3[2,16], pi_df3[2, 17], pi_df3[1, 17]), col="grey75")

text(2008,8, "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta))")

#plot credible and prediction intervals
# better plots

p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df2$year, 2018:2019), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df2$year[1:22], 
                         ymax = ci_df2[2, 1:22], 
                         ymin = ci_df2[1, 1:22]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2018:2019), 
                         ymax = pi_df2[2, 23:24], 
                         ymin = pi_df2[1, 23:24]),
                     alpha = 0.3)        
p


# plot posterior against expected distribution
alpha_post <- as.data.frame(run_RM2$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dnorm(0, 100^-2))
str(alpha_post)
dist <- as.data.frame(rnorm(10000000, 0, 10))
names(dist)[names(dist) == "rnorm(1e+07, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p


## M-Beta(uniform)----
# not sequential
#"delta + Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj"

m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- delta + alpha*TI[i]*(1-TI[i]/beta) + gamma*CO[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dgamma(2, 1/3) 
beta ~ dunif(0, 365) 
gamma ~ dnorm(0, 100^-2) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
df2 <- df1
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(0, 0)), #from capelin_larval_indices 
                   CO=c(df2$resids_lag, c(0, 0)), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mortality <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                      model.file = textConnection(m.mortality))

#UPDATE WITH MORE BURN INS

## M-Beta(uniform)-Plot 5----
# DIAGNOSTICS
print(run_mortality, intervals=c(0.025, 0.975), digits = 3)
plot(run_mortality)

out <- run_mortality$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_mortality$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma', 'delta')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mortality$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_mortality$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df2 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df2[1, ]

#generate prediciton intevals using N2_new
y_new = run_mortality$BUGSoutput$sims.list$N2_new
pi_df2 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df2$year,2018:2019),y_med,type='l', ylim = c(2,8))
points(df2$year, df2$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df2$year[1:22]), ci_df2[1, 1:22], lty = 3)
#lines(c(df2$year,2017:2018), ci_df2[2, ], lty = 3)
lines(c(df2$year[1:22]), ci_df2[2, 1:22], lty = 3)
#lines(c(df2$year,2017:2018), pi_df2[1, ], lty = 2)
lines(c(2018:2019), pi_df2[1, 23:24], lty = 2)
#lines(c(df2$year,2017:2018), pi_df2[2, ], lty = 2)
lines(c(2018:2019), pi_df2[2, 23:24], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df2[1, 23], pi_df2[2,23], pi_df2[2, 24], pi_df2[1, 24]), col="grey75")


beta_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$beta)
beta_prior <- as.data.frame(dnorm(0, 100^-2))
beta_prior <- as.data.frame(dunif(0, 365))
dist <- as.data.frame(runif(1000000, 0, 365))
head(dist)
names(dist)[names(dist) == "runif(1e+06, 0, 365)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = beta_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

## R/M3 model----
# not sequential
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])"

m.RM3 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   

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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2) 
gamma ~ dgamma(2, 1/3) 
delta ~ dgamma(8.1, 1/11.11) 
epsilon ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
#df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(0, 0)), #from capelin_larval_indices 
                   TI=c(df3$tice, c(0, 0)), #made up - need new data
                   CO=c(df3$resids_lag, c(0, 0)),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D'),
                model.file = textConnection(m.RM3))

#UPDATE WITH MORE BURN INS

##R/M3-plot 5----
# DIAGNOSTICS
print(run_RM3, intervals=c(0.025, 0.975), digits = 3)
plot(run_RM3)

out <- run_RM3$BUGSoutput 
str(out)
out$mean

# Asess mixing of chains
post <- run_RM3$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')
plot(post[, 'epsilon'], type = 'l')

# or
vars <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon')

MyBUGSChains(out, vars)
#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)

#overdispersion
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333 and values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes
str(out$mean)
F1 <- out$mean$expY
N2 <- out$mean$N2
D <- out$mean$D

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))


# Residuals v covariates
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "ST", ylab = "Pearson resids")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM3$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
str(run_RM3$BUGSoutput$sims.list)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df3[1, ]

#generate prediciton intevals using N2_new
y_new = run_RM3$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))


#PLOT credible and prediction intervals
plot(c(df3$year,2018:2019),y_med,type='l', ylim = c(2,9))
points(df3$year, df3$ln_biomass_med)
#points(df2$year, y_med[1:14], col = "blue")
#points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
#lines(c(df2$year,2017:2018), ci_df2[1, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[1, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), ci_df3[2, ], lty = 3)
lines(c(df3$year[1:15]), ci_df3[2, 1:15], lty = 3)
#lines(c(df3$year,2017:2018), pi_df3[1, ], lty = 2)
lines(c(2018:2019), pi_df3[1, 16:17], lty = 2)
#lines(c(df3$year,2017:2018), pi_df3[2, ], lty = 2)
lines(c(2018:2019), pi_df3[2, 16:17], lty = 2)

polygon(x = c(2018, 2018, 2019, 2019), y = c(pi_df3[1, 16], pi_df3[2,16], pi_df3[2, 17], pi_df3[1, 17]), col="grey75")

text(2008,8, "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta + epsilon*resids_adj))")

#plot credible and prediction intervals
# better plots

p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df2$year, 2018:2019), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df2$year[1:22], 
                         ymax = ci_df2[2, 1:22], 
                         ymin = ci_df2[1, 1:22]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2018:2019), 
                         ymax = pi_df2[2, 23:24], 
                         ymin = pi_df2[1, 23:24]),
                     alpha = 0.3)        
p


# plot posterior against expected distribution
alpha_post <- as.data.frame(run_RM3$BUGSoutput$sims.list$alpha)
alpha_prior <- as.data.frame(dnorm(0, 100^-2))
str(alpha_post)
dist <- as.data.frame(rnorm(10000000, 0, 10))
names(dist)[names(dist) == "rnorm(1e+07, 0, 10)"] <- "v2" 

p <- ggplot()
p <- p + geom_histogram(data = alpha_post, aes(x=V1, y = ..ncount..), colour="black", fill="white") 
p <- p + geom_density(data = dist, aes(v2), colour = "red")
p

