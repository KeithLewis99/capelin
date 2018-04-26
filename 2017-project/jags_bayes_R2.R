#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-04-25, R version 3.3.3 (2017-03-06)             #

# this file is to confirm that the function "rsq_bayes" in JAGS produces results that are equivalent to Gelmen's "bayes_R2" written in STAN: they are now!!!

library(R2jags)
library(rjags)

source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

m.jags = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*X[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
#N2[i] ~ dnorm(mu[i])

# 3. Discrepancy measures
expY[i] <- mu[i]
Res[i] <- (N2[i] - expY[i])

#Squared residuals
}
#Sum of squared Pearson residuals:

# 2. Priors
alpha ~ dnorm(0, 25) 
beta ~ dnorm(1, 25)
sigma ~ dunif(0, 100) 
}'

# in stan, the priors are specified in terms of a "location" and "scale" parameter.  According to wikipedia "https://en.wikipedia.org/wiki/Scale_parameter" and "https://en.wikipedia.org/wiki/Scale_parameter", the location parameter is .... which is equivalent to the mean and sd (or sometimes var).  In stan, this is the sd - see ?stan_glm
# JAGS uses "precision" which is 1/var
# Gelmen uses scale = 0.2.  So scale = 1/sqrt(precision) therefore precision = 1/scale^2 and therefore 1/0.2^2 = 25: 
# Gelmen also uses a flat uniform prior which is essentially what sigma is in the JAGS code.
# see also "https://en.wikipedia.org/wiki/Location%E2%80%93scale_family"
# the results are virtually identical

model_data <- list(N2 = xy$y, 
                   X= xy$x,
                   N = nrow(xy))

run_jags <- jags(data=model_data,
                    parameters.to.save = c('mu', 'sigma', 'N2', 'alpha', 'beta', 'Res'),
                    model.file = textConnection(m.jags))

rsq_bayes <- function(ypred = ypred, out=out){
  #browser()
  # variance of predicted values: use 1:15 bc 15 observations of capelin
  y_var = apply(y_pred, 1, 'var')
  y_Var <- sum(y_var)/nrow(y_pred) # this is the sum of variances from each posterior draw / the number of posterior draws
  
  # variance of residuals
  res_pred = out$BUGSoutput$sims.list$Res # as above but for resids
  #res_med = apply(res_pred,2,'median')
  res_var = apply(res_pred,1,'var')
  
  #res_Var <- sum(res_var)/nrow(y_pred)
  
  # Bayesian R-squared
  #bay_rsq <- y_Var/(sum(res_var + y_var)/nrow(y_pred))
  bay_rsq <- y_var/(res_var + y_var)
  
  return(bay_rsq)
}

y_pred = run_jags$BUGSoutput$sims.list$mu
R1_r2 <- rsq_bayes(ypred = y_pred, out = run_jags)
print(median(R1_r2))
print(sd(R1_r2))
