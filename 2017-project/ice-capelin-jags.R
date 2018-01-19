
library(R2jags)
library(rjags)
library(readr)

rm(list=ls())
## load data----
test <- read.csv('figs/covariates/capelin_covariates_2001.csv',header=T)
str(test)

sheep = read.csv('D:/Keith/R/time-series/ecots/data/sheep.csv')
str(sheep)
# set values----
y = test$logcapelin
x = test$tice
N <- nrow(test)



## test 1 first try----
jags_code = '
model {
     # Likelihood
     for (i in 1:N) {
     y[i] ~ dnorm(alpha*x[i], sigma^-2)
     }
     # Priors
  alpha ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 100)
}'

model_data <- list(y = y, x = x, N = N)
model_parameters <- c("alpha", "sigma")

jags_run <- jags(data=model_data,
                 parameters.to.save = model_parameters,
                 model.file = textConnection(jags_code))

print(jags_run)
plot(jags_run)

post <- jags_run$BUGSoutput$sims.matrix
head(post)

plot(post[,'alpha'], type = 'l')

alpha_mean <- mean(post[,'alpha'])
plot(test$tice, test$capelin)
lines(test$tice, alpha_mean*test$tice, col = 'red')

apply(post, 2, quantile, probs = c(0.025, 0.975))
hist(post[, 'alpha'], breaks = 30)


## test 2 simple model----
jags_code = '
model {
     # Likelihood
     for (i in 1:N) {
     y[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)
     }
     # Priors
  alpha ~ dnorm(0, 100^-2)
  beta ~ dnorm(200, 10)     
  sigma ~ dunif(0, 100)
}'

model_data <- list(y = y, x = x, N = N)
model_parameters <- c("alpha", "beta", "sigma")

jags_run <- jags(data=model_data,
                 parameters.to.save = model_parameters,
                 model.file = textConnection(jags_code))

print(jags_run)
plot(jags_run)

post <- jags_run$BUGSoutput$sims.matrix
head(post)

plot(post[,'alpha'], type = 'l')
plot(post[,'beta'], type = 'l')

head(post[,'alpha'])
alpha_mean <- mean(post[,'alpha'])
beta_mean <- mean(post[,'beta'])
pred <- alpha_mean*test$tice*(1-(test$tice/beta_mean))

cbind(test$year, pred)
plot(test$year, test$logcapelin, ylim = c(0,10))
lines(test$year, pred, col = 'red')

ci_test <- apply(post, 2, quantile, probs = c(0.025, 0.975))
hist(post[, 'alpha'], breaks = 30)

alpha_2.5 <- ci_test[1,1]
beta_2.5 <- ci_test[1,2]
alpha_97.5 <- ci_test[2,1]
beta_97.5 <- ci_test[2,2]

ci_2.5 <- alpha_2.5*test$tice*(1-(test$tice/beta_2.5))

ci_97.5 <- alpha_97.5*test$tice*(1-(test$tice/beta_97.5))

lines(test$year, ci_2.5, col = 'blue')
lines(test$year, ci_97.5, col = 'blue')

## test 3 y_pred----
jags_code = '
model {
     # Likelihood
     for (i in 1:N) {
     y[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)
     y_pred[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)

     }
     # Priors
  alpha ~ dnorm(0, 100^-2)
  beta ~ dnorm(200, 100^-2)     
  sigma ~ dunif(0, 100)
}'


model_data <- list(y = y, x = x, N = N)
model_parameters <- c(c('y_pred'))

jags_run <- jags(data=model_data,
                 parameters.to.save = model_parameters,
                 model.file = textConnection(jags_code))

print(jags_run)
plot(jags_run)

#Posterior predictive outputs
str(jags_run)
pars = jags_run$BUGSoutput$sims.list$y_pred
head(pars)
plot(test$year, apply(pars,2,'mean'), ylim = c(4, 6.5))
lines(test$year, pred, col = 'red')


## test 4 2 step----
jags_code = '
model {
# Likelihood
for (i in 1:N) {
y[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)
y_one_ahead[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)
}
for (i in 2:N) {
y_two_ahead[i] ~ dnorm(alpha*x[i]*(1-(x[i]/beta)), sigma^-2)

}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dnorm(200, 100^-2)     
sigma ~ dunif(0, 100)
}'


model_data <- list(y = y, x = x, N = N)
model_parameters <- c(c('y_one_ahead', 'y_two_ahead'))

jags_run <- jags(data=model_data,
                 parameters.to.save = model_parameters,
                 model.file = textConnection(jags_code))

one_ahead = jags_run$BUGSoutput$sims.list$y_one_ahead
two_ahead = jags_run$BUGSoutput$sims.list$y_two_ahead
plot(test$year, apply(one_ahead,2,'mean'), ylim = c(4, 6.5))
points(test$year[2:nrow(test)], apply(two_ahead,2,'mean'), col = 'blue', pch = 19)
lines(test$year, pred, col = 'red')

## test 5 NA trick----
jags_code = '
model {
     # Likelihood
     for (i in 1:N) {
          mu[i] <- alpha*x[i]*(1-(x[i]/beta))
          y[i] ~ dnorm(mu[i], sigma^-2)
     }
     # Priors
     alpha ~ dnorm(0, 100^-2)
     beta ~ dnorm(200, 10)     
     sigma ~ dunif(0, 100)
}'

num_forecasts = 2 # 10 extra years
model_data <- list(y = c(test$logcapelin, rep(NA, num_forecasts)), 
                   x=c(test$tice, c(65, 77)),
                   N = nrow(test) + num_forecasts)

jags_run <- jags(data=model_data,
                 parameters.to.save = c('mu', 'sigma'),
                 model.file = textConnection(jags_code))

y_pred = jags_run$BUGSoutput$sims.list$mu
y_med = apply(y_pred,2,'median')
plot(c(test$year,2017:2018),y_med,type='l', ylim=c(3,7))
points(test$year, y_med[1:14])
points(c(2017:2018), y_med[15:16], col = 'red', pch=19)
# this makes the credible interval 90%
lines(c(test$year,2017:2018), ci_test[1, ], lty = 3)
lines(c(test$year,2017:2018), ci_test[2, ], lty = 3)


str(test)
apply(y_pred,2,'quantile', c(0.05, 0.95))[, 15:16] #these are the extra 2 years
ci_test <- apply(y_pred,2,'quantile', c(0.05, 0.95))
ci_test[1, ]


## test gamma on simple model----

jags_code = '
model {
# Likelihood
for (i in 1:N) {
mu[i] <- alpha*x[i]*(1-(x[i]/beta))
y[i] ~ dnorm(mu[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dgamma(5, 2)     
sigma ~ dunif(0, 100)
}'

num_forecasts = 2 # 10 extra years
model_data <- list(y = c(test$logcapelin, rep(NA, num_forecasts)), 
                   x=c(test$tice, c(65, 77)),
                   N = nrow(test) + num_forecasts)

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


## sequential model----

jags_code = '
model {
# Likelihood
for (i in 1:N) {
mu[i] <- alpha*x[i]*(1-(x[i]/beta))
y[i] ~ dnorm(mu[i], sigma^-2)
}
# Priors
alpha ~ dnorm(0, 100^-2)
beta ~ dgamma(5, 2)     
sigma ~ dunif(0, 100)
}'

num_forecasts = 2 # 10 extra years
model_data <- list(y = c(test$logcapelin, rep(NA, num_forecasts)), 
                   x=c(test$tice, c(65, 77)),
                   N = nrow(test) + num_forecasts)

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
