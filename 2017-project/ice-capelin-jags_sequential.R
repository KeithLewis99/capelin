#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-01-15, R version 3.3.3 (2017-03-06)             #

# The ice-capelin-covariate file reveled some problems with the approach to date.  We had hoped to create "lumpy vaults" that would help us predict capelin abundance.  However the we merely created tilted vaults.  This does little to improve our understanding over Ale's original tice model.  Paul then suggested a sequential approach in a Bayesian framework.  Model recuitment, convert to spawning biomass, and then predict adult mortality.
# The purpose of this file is to:
# 1) build these sequential models in the Bayesian framework using JAGS
# 2) estimate the number of immature capelin in teh age 2 class
# However it is not clear if this will work conceptually because of the timing and limitations of the capelin surveys.

#So now, we are simply using JAGS as a convenient way to calculate confidence/credible intervals and prediction intervals

# The purpose of this file is to implement
library(R2jags)
library(rjags)
library(readr)
library(dplyr)

rm(list=ls())
## load data----
# this from the ice-capelin-covariates file
df <- read.csv('figs/covariates/capelin_covariates_2001.csv',header=T)
str(df)

df1 <- df[c("year", "logcapelin", "tice", "resids_adj", "Ssurface_tows_lag2", "ps_meanTot_lag2")]
str(df1)

df1$resids_adj[14] <- 35



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

## test 5----
# not sequential
#"delta + Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj"
m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- delta + alpha*TI[i]*(1-(TI[i])*beta) + gamma*RA[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2)
}
# 2. Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
delta ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100)
}'


m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- delta + alpha*TI[i]*(1-(TI[i])*beta) + gamma*RA[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2)

# 3. Discrepancy measures
YNew[i] ~dnorm(mu[i], 100^-2)
expY[i] <- mu[i]
varY[i] <- mu[i]
PRes[i] <- (N2[i] - expY[i]) / sqrt(varY[i])
PResNew[i] <- (YNew[i] - expY[i]) / sqrt(varY[i])
D[i] <- pow(PRes[i], 2)
DNew[i] <- pow(PResNew[i], 2)
}

Fit <- sum(D[1:N])
FitNew <- sum(DNew[1:N])

# 2. Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2)
delta ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100)
}'

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df1$logcapelin, rep(NA, num_forecasts)), 
                   TI=c(df1$tice, c(70, 70)), #from capelin_larval_indices 
                   RA=c(df1$resids_adj, c(15, 15)), #made up - need new data
                   N = nrow(df1) + num_forecasts)

run_mortality <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma', 'N2_new', 'alpha', 'beta', 'gamma', 'delta'),
                 model.file = textConnection(m.mortality))

#UPDATE WITH MORE BURN INS

# DIAGNOSTICS
print(run_mortality, intervals=c(0.025, 0.975), digits = 3)
plot(run_mortality)

# Asess mixing of chains
post <- run_mortality$BUGSoutput$sims.matrix
head(post)
plot(post[, 'alpha'], type = 'l')
plot(post[, 'beta'], type = 'l')
plot(post[, 'gamma'], type = 'l')
plot(post[, 'delta'], type = 'l')

#overdispersion


# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mortality$BUGSoutput$sims.list$mu
str(y_pred)
head(y_pred)
#look at output
y_med = apply(y_pred,2,'median')
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] 
ci_df1 <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_df1[1, ]

#generate prediciton intevals using N2_new
y_new = run_mortality$BUGSoutput$sims.list$N2_new
pi_df1 <- apply(y_new,2,'quantile', c(0.05, 0.95))

#PLOT credible and prediction intervals
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

text(2009,8, "log(caplein) = delta + Alpha*tice*(1-(tice/Beta)) + Gamma*resids_adj")

#plot credible and prediction intervals
# better plots
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = df1, 
                    aes(y = logcapelin, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("Year") + ylab("ln(capelin)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(aes(x = c(df1$year, 2017:2018), 
                       y = y_med))
p <- p + geom_ribbon(aes(x = df1$year, 
                         ymax = ci_df1[2, 1:14], 
                         ymin = ci_df1[1, 1:14]),
                     alpha = 0.5)
p <- p + geom_ribbon(aes(x = c(2017:2018), 
                         ymax = pi_df1[2, 14:16], 
                         ymin = pi_df1[1, 14:16]),
                     alpha = 0.3)        
p
