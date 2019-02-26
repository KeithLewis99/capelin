#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-03-19, R version 3.3.3 (2017-03-06)   

# The purpose of this file is to validate the models in "ice-capelin-sequential-jag.R".  Specifically, the purpose is to test the influence of any one year and how much the absence of that year influences the ability of the model to predict the year.  This is a jack-knife type of cross validation. See Olden et al. 2002.


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

## read in source code (Zuur and Lewis functions)-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")


## data----
filepath_gen <- "Bayesian/biomass_cond_ag1_2_DIC_new/" #"ag2_cond_ag1_DIC" "Bayesian/biomass_cond_ag1_DIC/"
folder <- "sensitivity"

dataset <- "biomass" #"age"
if(dataset == "biomass"){
     df <- read_csv(paste0(filepath_gen, "all_data_a1.csv"))
     df2 <- df
     #View(df2)
     df3 <- subset(df2, year>2002)
     yaxis1 = "ln_biomass_med" 
     ylab1 = "ln(capelin biomass(ktons))"
     var_x <- df3$ln_biomass_med
     
} else if(dataset == "age") {
     df <- read_csv("Bayesian/age2_cond_ag1/all_data_a1.csv")
     df2 <- df
     #View(df2)
     df3 <- subset(df2, year>2002)
     yaxis1 = "age2_log10" 
     ylab1 = "log10 capelin - age2"
     var_x <- df3$age2_log10
}

##lists----
# Make list of dataframe with JUST the knockout year (this may not be needed)
#make the list of dataframes
insert_ <- rep(list(list()), 15)
#fill the list
for (i in 1:15){
     df3_temp <- df3[i, ]
     insert_[[i]] <- df3_temp
     names(insert_) <- 2003:2017
     list(insert_ = insert_)
}

# Make list of dataframe without the knockout year (this may not be needed)
#make the list of dataframes
df3_ <- rep(list(list()), 15)

# fill the list
for(i in 1:15){
     df3_temp <- df3[-i, ]
     df3_[[i]] <- df3_temp
     names(df3_) <- 2003:2017
     list(df3_=df3_)
}

# make a vector of years
insert_year <- c(2003:2017)     

# make list of values for the "prediction" of each year
# make vectors
st_ <- rep(NA, 15)
#ti_ <- rep(NA, 15)
#co_ <- rep(NA, 15)
zo_ <- rep(NA, 15)

#fill vectors in a list
for(i in 1:15){
     st <- df3[c("year", "surface_tows_lag2")][i,2]
 #    ti <- df3[c("year", "tice")][i,2]
  #   co <- df3[c("year", "meanCond_lag")][i,2]
     zo <- df3[c("year", "ps_meanTot_lag2")][i,2]
     st_[i] <- st
   #  ti_[i] <- ti
    # co_[i] <- co
     zo_[i] <- zo
     names(st_) <- 2003:2017
     #names(ti_) <- 2003:2017
     #names(co_) <- 2003:2017
     names(zo_) <- 2003:2017
     #list(st_=st_, ti_=ti_, co_=co_)
     list(st_=st_, zo_=zo_)
}

## Recruitment----
#"alpha + beta*ST[i] + gamma*PS[i]"

m.recruit = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*PS[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   
log_lik[i] <- logdensity.norm(N2[i], mu[i], sigma^-2)

# 3. Discrepancy measures
expY[i] <- mu[i]
varY[i] <- sigma^2
Res[i] <- (N2[i] - expY[i])
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
alpha ~ dnorm(0, 100^-2) # Int
beta ~ dnorm(0, 100^-2)  # larval abundance
gamma ~ dnorm(0, 100^-2) #Psuedocalanus
sigma ~ dunif(0, 100) 
}'


num_forecasts =  1# 1 extra years

##lists for models----
# create a bunch of empty lists
run_RM3_ <- rep(list(list()), 15)
y_pred_ <- rep(list(list()), 15)
y_med_ <- rep(list(list()), 15)
ci_df3_ <- rep(list(list()), 15)
p_med_ <- rep(list(list()), 15)
pi_df3_ <- rep(list(list()), 15)
err_df3_ <- rep(list(list()), 15)
err_new_df3_ <- rep(list(list()), 15)
err_mean_new_df3_ <- rep(list(list()), 15)

##loop for knockout----
for(i in 1:15){
     model_data <- list(N2 = c(df3_[[i]]$ln_biomass_med, rep(NA, num_forecasts)), #var_x spacified in data section - ln_biomass_med
                        ST=c(df3_[[i]]$surface_tows_lag2, st_[[i]]),
                        #TI=c(df3_[[i]]$tice, ti_[[i]]), 
                        #CO=c(df3_[[i]]$meanCond_lag, co_[[i]]),
                        PS=c(df3_[[i]]$ps_meanTot_lag2, zo_[[i]]),
                        N = nrow(df3_[[i]]) + num_forecasts)
     run_RM3 <- jags(data=model_data,
                         parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', 'log_lik'),
                         model.file = textConnection(m.recruit))
     #run_RM3_[[i]] <- run_RM3   
      #UPDATE WITH MORE BURN INS
     #run_RM3 <-update(run_RM3, n.iter = 300000, n.thin = 50, n.burnin = 100000)
     run_RM3_[[i]] <- run_RM3
     y_pred_ = run_RM3$BUGSoutput$sims.list$mu
     y_pred_[, i]
     y_med_[[i]] = apply(y_pred_,2,'median')
     ci_df3_[[i]] <- apply(y_pred_,2,'quantile', c(0.025, 0.975))
     y_new_ = run_RM3$BUGSoutput$sims.list$N2_new
     p_med_[[i]] = apply(y_new_,2,'median')
     pi_df3_[[i]] <- apply(y_new_,2,'quantile', c(0.025, 0.1, 0.9, 0.975))
     err_df3_[[i]] <- insert_[[i]]$ln_biomass_med - y_pred_[, 15]
     err_new_df3_[[i]] <- insert_[[i]]$ln_biomass_med - y_new_[, 15]
     err_mean_new_df3_[[i]] <- insert_[[i]]$ln_biomass_med - mean(y_new_[, 15])
     list(run_RM3_[[i]], y_med_[[i]], ci_df3_[[i]], p_med_[[i]], pi_df3_[[i]], err_df3_[[i]], err_new_df3_[[i]], err_mean_new_df3_[[i]])     
}     

# preditive score for RM3 - to get all scores, run code but with different models - do it brute force for now to save time
# deviation for the estiamte from the missing year to the value of the missing year
all_err <- unlist(err_df3_)
mse <- mean(all_err ^ 2, na.rm = T)
mse

# CI for mse_new
mse_ci <- quantile(all_err^ 2, c(0.025, 0.975), na.rm = T)

#Paul thinks that this is the right value
all_err <- unlist(err_mean_new_df3_)
mse <- mean(all_err ^ 2, na.rm = T)
mse
mse_ci <- quantile(all_err^ 2, c(0.1, 0.8), na.rm = T)
mse_ci


for(i in 1:15){
     plotCredInt3(df3_[[i]], insert_[[i]], 
                  yaxis = yaxis1, 
                  ylab = ylab1, 
                  y_line = y_med_[[i]], 
                  ci = ci_df3_[[i]], 
                  dpp = p_med_[[i]],
                  dpi = pi_df3_[[i]], 
                  insert_year = insert_year[i], type = "CI")
     ggsave(paste0(filepath_gen, folder, "/CS1/credInt", insert_year[i], "v2.png"), width=10, height=8, units="in")
}
