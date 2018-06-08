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
library(magrittr)
library(loo)
library(forcats)

rm(list=ls())

## read in source code-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

# make new folders
folder_path1 <- "Bayesian/biomass_cond_ag1_2_DIC/" #"Bayesian/biomass_cond_ag1_DIC/""Bayesian/ag2_cond_ag1_DIC/"

make_direct1(folder_names, folder_path1)

## Set parameters for analysis----
capelin_data_set <- "biomass" # type of capelin index - alt value ==  "biomass" "age"
#cond_data_a1_2 <- "cond" # type of condition index - alt value == "resids" or "cond_dat_a1_2"
dic_run <- "yes" # is this a DIC run? "yes" or "no"

# sets the values for the predition interval
# these values are normalized.  The normalization was performed after the data were all brought together.  The calculations are performed on lines 159-164 unless stated otherwise
STpred <- c(-1.1798, -0.5525)
PSpred <- c(-0.024866638, 0) 
TIpred <- c(0.57, 0.788)  #first value estimated from this year (line 71-72), second, mean from centered data (line 168)

x <- "condition" #"resids"
y <- "a1" # "a1_2"
COpred <- if(x=="condition" & y == "a1"){
     COpred <- c(1.67, 0) # from lines 159-164
} else if(x=="condition" & y == "a1_2"){
     COpred <- c(1.90, 0) # from lines 159-164
} else if(x =="resids"){
          COpred <- c(0, 0)
}


# sets the subfolder
filepath_gen <- "biomass_cond_ag1_2_DIC"  # for filepath for datasets and pairs Plots - "biomass_cond_ag1_DIC" "ag2_cond_ag1_DIC"

# changes of the axis and axis labels for credInt
ylab1 = "ln(capelin biomass(ktons))" #alt: "log10 capelin - age2"  "ln(capelin biomass(ktons))"
yaxis1 = "ln_biomass_med" # alt: "age2_log10" "ln_biomass_med"

# change "type = [something]" or = NA  
# change "CO=c(df3$meanCond_lag alt: resids_lag


# tice Date
# to inform TIpred values
tice_day <- as.Date('2018-02-26') # this is the lowest ice extent as of 2018-04-17
lubridate::yday(tice_day)

## load data----
## capelin (1985-2017)----
# source(biomass): "Age disaggregate abundance for Keith Lewis - 2017 added_v1.xlsx"
# source(age2): "capelin_age_disaggregate_abundance.xlsx" worksheet:Age disagg acoustic index
cap <- capelin_data(capelin_data_set)

## ice (1969-2017)----
#source: ice-chart-processing-data-v3.R
ice <- read_csv('output-processing/capelin-m1.csv')
glimpse(ice)
ice1 <- subset(ice, year > 2002)
mean(ice1$tice/100)
var(ice1$tice/100)

sd(ice1$tice/100)
sh = var(ice1$tice/100)/mean(ice1$tice/100)
scale = mean(ice1$tice/100)^2/var(ice1$tice/100)
ra = 1/scale

## condition (1995-2017)-----
# source: for "cond" see capelin-condition.R for original and derived datasets
# source: for "resids" see Fran's original "capelin_condition_maturation.xlsx"
#cond <- condition_data(cond, 'data/condition_ag1_out.csv')
if(capelin_data_set == "age"){
     cond_dat <- read_csv('data/condition_ag1_MF_out.csv')
} else if (capelin_data_set == "biomass"){
     cond_dat <- read_csv('data/condition_ag1_2_MF_out.csv')     
}
cond_dat$meanCond_lag <- lag(cond_dat$meanCond, 1)
#cond_dat$medCond_lag <- lag(cond_dat$medCond, 1)
#condResids <- condition_data("resids", 'data/condition_ag1_out.csv')


## larval data (2003-2017)----
# source "capelin_age_disaggregate_abundance.xlsx":Larval indices
larvae <- read_csv('data/capelin_larval_indices.csv')
#glimpse(larvae)
larvae$surface_tows_lag2 <- lag(larvae$surface_tows, 2)
#View(larvae)

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
capelin_join <- left_join(capelin_join, cond_dat, by = "year")
#capelin_join <- left_join(capelin_join, condResids, by = "year")
capelin_join <- left_join(capelin_join, larvae, by = "year")
capelin_join <- left_join(capelin_join, ps_tot, by = "year")

df <- capelin_join
glimpse(df)
if(capelin_data_set == "biomass"){
     df[c("biomass_med", "ln_biomass_med", "tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")]
} else if (capelin_data_set == "age"){
     df[c("age2", "age2_log10", "tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")]    
}


# normalize the data set----
#df1 <- subset(df, year>1995)
df1 <- subset(df, year>1998)
glimpse(df1)

if(x == "condition"){
     cols <- c("meanCond_lag", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")
} else if (x == "resids"){
     cols <- c("meanCond_lag", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")
}

df1 %<>%
     mutate_each_(funs(scale),cols)
glimpse(df1)
df1$meanCond
df1$meanCond_lag
co_temp <- df1[c("year", "meanCond")]
co_temp$meanCond_n <- scale(co_temp$meanCond)
ps_temp <- df1[c("year", "ps_meanTot")]
ps_temp$ps_meanTot_n <- scale(ps_temp$ps_meanTot)
st_temp <- df1[c("year", "surface_tows")]
st_temp$surface_tows_n <- scale(st_temp$surface_tows)

# divide by 100
df1$tice <- df1$tice/100
mean(df1$tice)

write_csv(df1, paste0("Bayesian/", filepath_gen, "/all_data_a1.csv"))
# shrink df - not normalized but small and easy to interpret parameters


# pairs-plot: relationships and correlations among RV and EV----
df_name <- name_pairPlot(df1, "age")
z <- "ln_biomass_med" # alt "log10_age2 index" ln capelin biomass

pdf(paste0("Bayesian/", filepath_gen, "/pairs_resids.pdf"))
pairs.panels(df_name[c(z, "tice", "condition l1", "larval abundance l2", "zooplankton abun l2")], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses,
             cex.labels = 1,
             cex.cor = 1
)
dev.off()

#######Set data sets----
# subset the data depending on whether the analysis if for "biomass" or "age" and if it is for a DIC run ("yes" or "no")
if(dic_run == "yes" & capelin_data_set == "biomass"){
     df3 <- subset(df1, year>2002)
     df2 <- subset(df1, year>2002)
     
     pred3 <- df3$ln_biomass_med 
     pred2 <- pred3
     
} else if (dic_run == "no" & capelin_data_set == "biomass"){
     df3 <- subset(df1, year>2002)
     df2 <- df1
     
     pred3 <- df3$ln_biomass_med 
     pred2 <- df2$ln_biomass_med 
     
} else if (dic_run == "yes" & capelin_data_set == "age"){
     df3 <- subset(df1, year>2002)
     df2 <- subset(df1, year>2002)
     
     pred3 <- df3$age2_log10
     pred2 <- pred3
     
} else if (dic_run == "no" & capelin_data_set == "age"){
     df3 <- subset(df1, year>2002)
     df2 <- df1
     
     pred3 <- df3$age2_log10
     pred2 <- df2$age2_log10
}

#########Bayesian models#############################----
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   PS=c(df3$ps_meanTot_lag2, PSpred), #see df_norm
                   N = nrow(df3) + num_forecasts)

run_recruit <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', 'log_lik'),
                      model.file = textConnection(m.recruit))

#UPDATE WITH MORE BURN INS

##R-Diagnostics----
print(run_recruit, intervals=c(0.025, 0.975), digits = 3)
out <- run_recruit$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma')
filepath <- paste0(filepath_gen, "/recruitment_1")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data") # should follow the line
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
df_diag <- as.data.frame(model_data)
df_diag <- cbind(df_diag, E1)
plot(df_diag$ST, df_diag$EI, xlab = "Surface Tows", ylab = "Pearson resids")
plot(df_diag$PS, df_diag$EI, xlab = "Pseudo-calanus", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_recruit$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_recruit$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))


## R-Results----
# Get the log likelihood
log_lik_r1 = run_recruit$BUGSoutput$sims.list$log_lik
w_r1 <- waic(log_lik_r1)
w_r1 <- loo(log_lik_r1)
plot(w_r1)
x <- print(w_r1$paretok_k)
cbind(w_r1$pointwise, x)
str(w_r1)
pareto_k_table(w_r1)
pareto_k_ids(w_r1)
print(w_r1$pareto_k)
compare(w_r1, w_m1)
dic_R1 <- dic.samples(run_recruit$model, n.iter=1000, type="pD")
str(x)
dic_R1sum <- sum(dic_R1$deviance)

R1_r2m <- rsq_bayes(ypred = y_pred, out = run_recruit)
R1_r2 <- c(median(R1_r2m), sd(R1_r2m))

#PLOT credible and prediction intervals
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*pseudocal_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))


# plot posterior against expected distribution
# gamma dist
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)


#alpha
priormean <- 0
priorsd <- 10
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#beta
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## Mortality----
#"alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]"
m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#mortality
mu[i] <- alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]
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
alpha ~ dnorm(0, 20^-2) 
beta ~ dgamma(5, 1/3) 
gamma ~ dgamma(11.5, 5.7)
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 10) 
}'

# alpha and delta are based on normalized means and large uniformed priors
#beta based on Bolker pg 132 - Fig4.13 - trying for an uniformative prior 
#gamma: see gammaDist_mode_shape_rate_calculator.R - 
# mode = 1.86 - based on Ale's value of 93 for beta*2 -> 186: but this is the mode - Ale's value is the maximum
# sd = 1
#shape(a) is mean^2/var; scale(s) equal Var/mean - but we use rate so mean/Var
# sigma: uninformative for condition


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred2, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, TIpred), #from capelin_larval_indices 
                   CO=c(df2$meanCond_lag, COpred), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mortality <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', 'log_lik'),
                      model.file = textConnection(m.mortality))

run_mortality <-update(run_mortality, n.iter = 300000, n.thin = 50, n.burnin = 100000)

## M-DIAGNOSTICS----
print(run_mortality, intervals=c(0.025, 0.975), digits = 3)
# saved as table_1
out <- run_mortality$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta')
filepath <- paste0(filepath_gen, "/mortality_1")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# Need to implement the code for Cook's D.  The code is in ONeNote (Capelin- "To do list")
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
# this may not be quite right as N2 is not quite the observed values because NA has been estiamted.
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()


# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
df_diag <- as.data.frame(model_data)
df_diag <- cbind(df_diag, E1)
#Myxyplot(df_diag, MyVar, "E1")
plot(df_diag$TI, df_diag$EI, xlab = "Tice", ylab = "Pearson resids")
plot(df_diag$CO, df_diag$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mortality$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df2 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_mortality$BUGSoutput$sims.list$N2_new
pi_df2 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], paste0("Bayesian/", filepath, "/pi.csv"))

## M-Results----
dic_M1 <- dic.samples(run_mortality$model, n.iter=1000, type="pD")
str(x)
dic_M1sum <- sum(dic_M1$deviance)
M1_r2m <- rsq_bayes(ypred = y_pred, out = run_mortality)
M1_r2 <- c(median(M1_r2m), sd(M1_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. Also, is this matched up properly - i haven't used PanelNames
MyBUGSHist1(out, vars, transform = "yes")
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")


# plot the credible and prediction intervals
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
#txt <- "log(capelin) = alpha +  beta*tice*(1-(tice/gamma)) \n + delta*Condition(resids)"

plotCredInt(df2, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2008, y = 7, type = "CI")

ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig1(out$sims.list$beta, transform = "yes", parm = "slope")
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "width")
delta <- posterior_fig(out$sims.list$delta)


#alpha
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#beta
priormean <- 5
priorsd <- 1/3
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)/10
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#delta
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## M-(uniform)----
#"alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]"
#"alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]"
m.mort_unif = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dunif(0, 100) 
gamma ~ dunif(0, 10) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#alpha based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred2, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, TIpred), #from capelin_larval_indices 
                   CO=c(df2$meanCond_lag, COpred), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mort_unif <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                      model.file = textConnection(m.mort_unif))

run_mort_unif <-update(run_mort_unif, n.iter = 300000, n.thin = 50, n.burnin = 100000)


## M(uniform)-diagnostics----
print(run_mort_unif, intervals=c(0.025, 0.975), digits = 3)
out <- run_mort_unif$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta')
filepath <- paste0(filepath_gen, "/mortality_2")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# Need to implement the code for Cook's D.  The code is in ONeNote (Capelin- "To do list")
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
# this may not be quite right as N2 is not quite the observed values because NA has been estiamted.
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
df_diag <- as.data.frame(model_data)
df_diag <- cbind(df_diag, E1)
#Myxyplot(df_diag, MyVar, "E1")
plot(df_diag$TI, df_diag$EI, xlab = "Tice", ylab = "Pearson resids")
plot(df_diag$CO, df_diag$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mort_unif$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median')
ci_df2 <- apply(y_pred,2,'quantile', c(0.05, 0.95))

#generate prediciton intevals using N2_new
y_new = run_mort_unif$BUGSoutput$sims.list$N2_new
pi_df2 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], paste0("Bayesian/", filepath, "/pi.csv"))

## M(uniform)-Results----
# Get the log likelihood
dic_M_unif <- dic.samples(run_mort_unif$model, n.iter=1000, type="pD")
dic_M_unifsum <- sum(dic_M_unif$deviance)

M_unif_r2m <- rsq_bayes(ypred = y_pred, out = run_mort_unif)
M_unif_r2 <- c(median(M_unif_r2m), sd(M_unif_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist1(out, vars, transform = "yes")
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

txt <- "log(caplein) = alpha +  beta*tice*(1-(tice/gamma(unif))) \n + delta*Condition(resids)"

plotCredInt(df2, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2006, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

# R-squared - see Gelmen paper

# plot posterior against expected distribution - BUT DO THIS ONLY FOR INFORMATIVE PRIORS
# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig1(out$sims.list$beta, transform = "yes", parm = "slope")
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "width")
delta <- posterior_fig(out$sims.list$delta)

#alpha
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#beta
prior <- runif(n = 10000, 0, 100)/10
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma
prior <- runif(n = 10000, 0, 10)*100
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)
mean(gamma$df)/2

alpha ~ dnorm(0, 100^-2) 
beta ~ dunif(0, 100) 
gamma ~ dunif(0, 10) 
delta ~ dnorm(0, 100^-2) 

#delta
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## RM_1----
#"mu[i] <- alpha + beta*ST[i] + gamma*CO[i]"

summary(lm())
m.RM1 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*CO[i]
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2)
gamma ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

# priors are uninformative for condition

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   CO=c(df3$meanCond_lag, COpred), #see df_norm
                   N = nrow(df3) + num_forecasts)

run_RM1 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM1))


##RM_1-Diagnostics----
print(run_RM1, intervals=c(0.025, 0.975), digits = 3)
out <- run_RM1$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma')
filepath <- paste0(filepath_gen, "/rm_1")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
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
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
test <- as.data.frame(model_data)
test <- cbind(test, E1)
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM1$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_RM1$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))

## RM_1-Results----
dic_RM1 <- dic.samples(run_RM1$model, n.iter=1000, type="pD")
dic_RM1sum <- sum(dic_RM1$deviance)

RM1_r2m <- rsq_bayes(ypred = y_pred, out = run_RM1)
RM1_r2 <- c(median(RM1_r2m), sd(RM1_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*resids_adj"
#PLOT credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))


# plot posterior against expected distribution
# gamma dist
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)

#alpha
priormean <- 0
priorsd <- 10
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#beta
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## RM_2----
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta)"

m.RM2 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta)
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
alpha ~ dnorm(0, 100^-2) 
#alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2) 
gamma ~ dgamma(5, 1/3) 
delta ~ dgamma(11.5, 5.7)
#gamma ~ dunif(0, 100) 
#delta ~ dunif(0, 10) 
sigma ~ dunif(0, 10) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#gamma: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean

num_forecasts = 2 # 2 extra years

model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   TI=c(df3$tice, TIpred), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_RM2 <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                      model.file = textConnection(m.RM2))

#UPDATE WITH MORE BURN INS
run_RM2 <-update(run_RM2, n.iter = 300000, n.thin = 50, n.burnin = 100000)


##RM_2-Diagnostics----
print(run_RM2, intervals=c(0.025, 0.975), digits = 3)
out <- run_RM2$BUGSoutput 


#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta')
filepath <- paste0(filepath_gen, "/rm_2")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# Need to implement the code for Cook's D.  The code is in ONeNote (Capelin- "To do list")
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
# this may not be quite right as N2 is not quite the observed values because NA has been estiamted.
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()


# Residuals v covariates
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Surface tow", ylab = "Pearson resids")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM2$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median')
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95))

#generate prediciton intevals using N2_new
y_new = run_RM2$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))


## RM_2-Results----
dic_RM2 <- dic.samples(run_RM2$model, n.iter=1000, type="pD")
dic_RM2sum <- sum(dic_RM2$deviance)

RM2_r2m <- rsq_bayes(ypred = y_pred, out = run_RM2)
RM2_r2 <- c(median(RM2_r2m), sd(RM2_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist1(out, vars, transform = "yes")
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta))"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2010, y = 7.75, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))


# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "slope")
delta <- posterior_fig1(out$sims.list$delta, transform = "yes", parm = "width")

#alpha
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = -bin_1)

#beta
priorsd <- 10
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma
priormean <- 5
priorsd <- 1/3
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)/10
str(prior)
plot(density(prior))
plot(density(prior), xlim = c(0,5))

limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#delta
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
plot(density(prior))
plot(density(prior), xlim = c(0,2))

limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## R/M3 model----
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2) 
gamma ~ dgamma(5, 1/3) 
delta ~ dgamma(11.5, 5.7)
epsilon ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#delta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
mean(df3$tice)^2/var(df3$tice)
var(df3$tice)/mean(df3$tice)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   TI=c(df3$tice, TIpred), #made up - need new data
                   CO=c(df3$meanCond_lag, COpred),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM3))

#UPDATE WITH MORE BURN INS
run_RM3 <-update(run_RM3, n.iter = 300000, n.thin = 50, n.burnin = 100000)


##R/M3-plot 5----
# DIAGNOSTICS
print(run_RM3, intervals=c(0.025, 0.975), digits = 3)
out <- run_RM3$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon')
filepath <- paste0(filepath_gen, "/rm_3")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
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
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
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
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))


## RM_3-Results----
dic_RM3 <- dic.samples(run_RM3$model, n.iter=1000, type="pD")
dic_RM3sum <- sum(dic_RM3$deviance)

RM3_r2m <- rsq_bayes(ypred = y_pred, out = run_RM3)
RM3_r2 <- c(median(RM3_r2m), sd(RM3_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 

MyBUGSHist1(out, vars, transform = "yes")
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta)) + epsilon*CO"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2010, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.png"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))


# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "slope")
delta <- posterior_fig1(out$sims.list$delta, transform = "yes", parm = "width")
epsilon <- posterior_fig(out$sims.list$epsilon)

#alpha
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- abs(mean(alpha$df)/100)

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#beta
priorsd <- 10
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma

priormean <- 5/10
priorsd <- 1/3/10
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)/10
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#delta
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
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

ggsave(paste0("Bayesian/", filepath, "/priorPost.png"), width=10, height=8, units="in")


## R1- informative prior----
#"alpha + beta*ST[i] + gamma*PS[i]"

m.recruit.prior = '
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0.35, 7.25) 
gamma ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'

# For beta: I ran the same regression as Murphy et al. (2018) in ice-capelin-covariates but for normalized values.  I got somewhat different values. ST slope = 0.35, SE = 0.138 (inverse = 7.25)
# Where is this analysis????

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices 
                   PS=c(df3$ps_meanTot_lag2, PSpred), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_recruit_prior <- jags(data=model_data,
                    parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                    model.file = textConnection(m.recruit.prior))

#UPDATE WITH MORE BURN INS

##R1-infor_prior-Diagnostics----
print(run_recruit_prior, intervals=c(0.025, 0.975), digits = 3)
out <- run_recruit_prior$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma')
filepath <- paste0(filepath_gen, "/recruitment_2")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$PS, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_recruit_prior$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_recruit_prior$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))


## R1-inform_prior-Results----
# Get the log likelihood
dic_R2 <- dic.samples(run_recruit_prior$model, n.iter=1000, type="pD")
str(x)
dic_R2sum <- sum(dic_R2$deviance)

R2_r2m <- rsq_bayes(ypred = y_pred, out = run_recruit_prior)
R2_r2 <- c(median(R2_r2m), sd(R2_r2m))

#PLOT credible and prediction intervals
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

txt <- "log(caplein) = alpha +  beta(prior)*surface_tows_lag2 \n + gamma*pseudocal_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

# plot posterior against expected distribution - BUT DO THIS ONLY FOR INFORMATIVE PRIORS
# plot posterior against expected distribution
# gamma dist
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)


#alpha
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#beta
priormean <- 0.35
priorsd <- 0.138
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
plot(density(prior))
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## M_NULL----
# use the Paul/Ale measure of condition
#"beta*TI[i]*(1-TI[i]/gamma)"

m.mortality_null = '
model {
# 1. Likelihood
for (i in 1:N) {
#mortaliyt
mu[i] <- beta*TI[i]*(1-TI[i]/gamma)
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
beta ~ dgamma(5, 1/3) 
gamma ~ dgamma(2.8, 1)
sigma ~ dunif(0, 10) 
}'


#gamma based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 180 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
#df2 <- subset(df1, year>2002)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred2, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, TIpred), #from capelin_larval_indices 
                   N = nrow(df2) + num_forecasts)

run_mortality_null <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'beta', 'gamma', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', 'log_lik'),
                      model.file = textConnection(m.mortality_null))

run_mortality_null <- update(run_mortality_null, n.iter = 50000, n.thin = 50)


#Do we need to remove the burn-in?  How much?  UPDATE the chain? Figure out and apply to all

## M_NULL-DIAGNOSTICS----
print(run_mortality_null, intervals=c(0.025, 0.975), digits = 3)
# saved as table_1
out <- run_mortality_null$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('beta', 'gamma')
filepath <- paste0(filepath_gen, "/mortality_0")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# Need to implement the code for Cook's D.  The code is in ONeNote (Capelin- "To do list")
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
# this may not be quite right as N2 is not quite the observed values because NA has been estiamted.
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()


# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(1,1), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_mortality_null$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df2 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_mortality_null$BUGSoutput$sims.list$N2_new
pi_df2 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], paste0("Bayesian/", filepath, "/pi.csv"))

## M-Null_Results----
dic_Mo <- dic.samples(run_mortality_null$model, n.iter=1000, type="pD")
dic_Mosum <- sum(dic_Mo$deviance)

Mo_r2m <- rsq_bayes(ypred = y_pred, out = run_mortality_null)
Mo_r2 <- c(median(Mo_r2m), sd(Mo_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. Also, is this matched up properly - i haven't used PanelNames
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

# plot the credible and prediction intervals
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
#txt <- "log(caplein) = beta*tice*(1-(tice/gamma))"

plotCredInt(df2, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2006, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=11, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

# plot posterior against expected distribution
# gamma dist
beta <- posterior_fig1(out$sims.list$beta, transform = "yes", parm = "slope")
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "width")

#beta
priormean <- 5
priorsd <- 1/3
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)/10
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100
mean(gamma$df)/2

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p2, p3, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## R_Null----
#"alpha + beta*ST[i] + gamma*PS[i]"

m.recruit.null = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i]
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dnorm(0, 100^-2)
sigma ~ dunif(0, 100) 
}'


num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   N = nrow(df3) + num_forecasts)

run_recruit_null <- jags(data=model_data,
                    parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', 'log_lik'),
                    model.file = textConnection(m.recruit.null))

#UPDATE WITH MORE BURN INS

##R_NULL-Diagnostics----
print(run_recruit_null, intervals=c(0.025, 0.975), digits = 3)
out <- run_recruit_null$BUGSoutput 
out$mean

#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta')
filepath <- paste0(filepath_gen, "/recruitment_0")

MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(1,1), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "Tice", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_recruit_null$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median') # median values of y_pred
ci_df3 <- apply(y_pred,2,'quantile', c(0.05, 0.95)) #credible interval a subscript at end of code can pull out specific values [, 15:16]

#generate prediciton intevals using N2_new
y_new = run_recruit_null$BUGSoutput$sims.list$N2_new
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], paste0("Bayesian/", filepath, "/pi.csv"))

## R_NUll-Results----
dic_Ro <- dic.samples(run_recruit_null$model, n.iter=1000, type="pD")
dic_Rosum <- sum(dic_Ro$deviance)

Ro_r2m <- rsq_bayes(ypred = y_pred, out = run_recruit_null)
Ro_r2 <- c(median(Ro_r2m), sd(Ro_r2m))

#PLOT credible and prediction intervals
MyBUGSHist1(out, vars, transform = "yes")
ggsave(MyBUGSHist(out, vars), filename = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = yaxis1, 
            ylab = ylab1, 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

# plot posterior against expected distribution
# gamma dist
alpha <- posterior_fig1(out$sims.list$alpha, transform = "yes", parm = "slope")
beta <- posterior_fig1(out$sims.list$beta, transform = "yes", parm = "width")

#alpha
priormean <- 0
priorsd <- 10
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df)-0.3, max(alpha$df) + 0.3)
x_label <- "alpha"
bin_1 <- mean(alpha$df)/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#beta
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorpost.pdf"), width=10, height=8, units="in")


## DIC table ----
index <- c(1:9)
tab1 <- as.data.frame(rbind(
     dic_R1sum,
      dic_M1sum,
     dic_M_unifsum,
      dic_RM1sum,
      dic_RM2sum,
      dic_RM3sum,
     dic_R2sum, # informed recruitment 
     dic_Mosum, 
     dic_Rosum
      ))
str(tab1)
tab1 <- tibble::rownames_to_column(tab1)
tab1 <- cbind(index, tab1)
names(tab1)[names(tab1)=="V1"] <- "DIC"
names(tab1)[names(tab1)=="rowname"] <- "model"
tab1 <- tab1[order(tab1$DIC), , drop = F]
tab1$dDIC <- tab1$DIC-tab1$DIC[1]
tab1
write.csv(tab1, paste0("Bayesian/", filepath_gen, "/DIC.csv"))

########################Sensitivity Analysis------
# sets the values for the predition interval

filepath <- paste0(filepath_gen, "/sensitivity")

temp <- df1$tice
str(temp)

# this value is from line 71-72
temp[20] <- 0.57
mean(temp)
sd(temp)
tice_m <- round(mean(temp), 2)
tice_h <- round(mean(temp) + sd(temp), 2)
tice_l <- round(mean(temp) - 1*sd(temp), 2)
tice_vl <- round(mean(temp) - 2*sd(temp), 2)


# where did this value come from??  Normalized condition data
temp1 <- df1$meanCond_lag
temp1[20] <- COpred[1]
mean(temp1)
sd(temp1)
cond_m <- round(mean(temp1), 2)
cond_h <- round(mean(temp1) + sd(temp1), 2)
cond_l <- round(mean(temp1) - sd(temp1), 2)

#create an array of the possible value combinations
row.names <- c("high", "med", "low", "veryLow")
col.names <- "val"
y.names <- c("2018", "2019")
mat.names <- c("tice", "cond")

T <- c(tice_h, tice_m, tice_l, tice_vl)
C <- c(cond_h, cond_m, cond_l, NA)
#A <- array(c(T, C), dim = c(3, 1, 2, 2), dimnames = list(row.names, col.names, y.names, mat.names))
A <- array(c(T, C), dim = c(4, 1, 2), dimnames = list(row.names, col.names, mat.names))

# the values TIpred and COpred are from lines 46-48
comb_ls <- matrix(c("MO", TIpred[1], COpred[1],
                    "HH", A[1,1,1], A[1,1,2],
                    "HM", A[1,1,1], A[2,1,2],
                    "HL", A[1,1,1], A[3,1,2],
                    "MH", A[2,1,1], A[1,1,2],
                    "MM", A[2,1,1], A[2,1,2],
                    "ML", A[2,1,1], A[3,1,2],
                    "LH", A[3,1,1], A[1,1,2],
                    "LM", A[3,1,1], A[2,1,2],
                    "LL", A[3,1,1], A[3,1,2],
                    "vLH", A[4,1,1], A[1,1,2],
                    "vLM", A[4,1,1], A[2,1,2],
                    "vLL", A[4,1,1], A[3,1,2]), 
             nrow=13,
             ncol=3, byrow=T)

# model for sensitivity
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])"

RM3_sens = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
#mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])
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

num_forecasts = 1
# 2 extra years

# create a bunch of empty lists
run_RM3_ <- rep(list(list()), 13)
y_pred_ <- rep(list(list()), 13)
y_med_ <- rep(list(list()), 13)
ci_df3_ <- rep(list(list()), 13)
p_med_ <- rep(list(list()), 13)
pi_df3_ <- rep(list(list()), 13)

# a loop for automating the process of assessing sensitivity
# the values TIpred and COpred are from lines 46-48
for(i in 1:13){
     model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)), #10^pred3 puts it on arithmetic scale
                        ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                        TI=c(df3$tice, comb_ls[i,2]), #made up - need new data
                        CO=c(df3$meanCond_lag, comb_ls[i,3]),
                        N = nrow(df3) + num_forecasts)
#     model_data_[[i]] <- model_data
     run_RM3 <- jags(data=model_data,
                     parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
                     model.file = textConnection(RM3_sens))
     #list(model_data_)
     #run_RM3 <-update(run_RM3, n.iter = 300000, n.thin = 50, n.burnin = 100000)
     run_RM3_[[i]] <- run_RM3
     #list(run_RM3_)
#     y_pred_[[i]] = run_RM3_[[i]]$BUGSoutput$sims.list$mu
     y_pred_ = run_RM3$BUGSoutput$sims.list$mu
     y_med_[[i]] = apply(y_pred_,2,'median')
     ci_df3_[[i]] <- apply(y_pred_,2,'quantile', c(0.05, 0.95))
     y_new_ = run_RM3$BUGSoutput$sims.list$N2_new
     p_med_[[i]] = apply(y_new_,2,'median')
     pi_df3_[[i]] <- apply(y_new_,2,'quantile', c(0.05, 0.95))
     list(run_RM3_[[i]], y_med_[[i]], ci_df3_[[i]], p_med_[[i]], pi_df3_[[i]])
}

# create a table of prediction intervals
# create the blank matrix
model <- c("MO", "HH", "HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL", "vLH", "vLM", "vLL")
model <- c("MOD", "LG", "LM", "LP", "MG", "MM", "MP", "EG", "EM", "EP", "vEG", "vEM", "vEP")

#model <- c("MO", "MH", "MM", "ML")
pi_pt <- rep(NA, 13)
per_2_5 <- rep(NA, 13)
per_97_5 <- rep(NA, 13)
pi_tabl <- cbind(model, pi_pt, per_2_5, per_97_5)

# loop to collect all of the prection intervals
#for(i in seq(model)){
 #    pi_tabl[i, 2] <- round(p_med_[[i]][16], 2)
  #   pi_tabl[i, 3:4] <- round(pi_df3_[[i]][,16], 2)
#}

# loop to collect all of the prection intervals
for(i in seq(model)){
     pi_tabl[i, 2] <- p_med_[[i]][16]
     pi_tabl[i, 3:4] <- pi_df3_[[i]][,16]
}

pi_tabl
pi_tabl <- data.frame(pi_tabl)
str(pi_tabl)
pi_tabl$pi_pt <- as.numeric(levels(pi_tabl$pi_pt)[pi_tabl$pi_pt])
pi_tabl$pi_pt <- 10^pi_tabl$pi_pt
pi_tabl$per_2_5 <- as.numeric(levels(pi_tabl$per_2_5)[pi_tabl$per_2_5])
pi_tabl$per_2_5 <- 10^pi_tabl$per_2_5
pi_tabl$per_97_5 <- as.numeric(levels(pi_tabl$per_97_5)[pi_tabl$per_97_5])
pi_tabl$per_97_5 <- 10^pi_tabl$per_97_5

write.csv(pi_tabl, paste0("Bayesian/", filepath, "/pi_tabl.csv"))

# graph summarizing the influence of the knockouts
p <- ggplot(data = data.frame(pi_tabl), aes(x = fct_inorder(model))) 
p <- p + geom_point(aes(y=as.numeric(pi_pt)), size = 3)
p <- p + geom_errorbar(aes(ymax = as.numeric(per_97_5), ymin = as.numeric(per_2_5)), width = 0.5)
p <- p + xlab("Scenario") + ylab("Prediction point estimate \n & interval")
p <- p + theme_bw(base_size = 20) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
p <- p + geom_hline(aes(yintercept = pi_tabl$pi_pt[1]), colour = "red", size = 2)
#p <- p + scale_x_discrete(breaks = c("HH", "HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL"), labels = c("HH", "HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL"))
p
ggsave(paste0("Bayesian/", filepath, "/sensitivity_arith.pdf"), width=10, height=8, units="in")

# This shows the reason why the middle values for tice are slightly higher than the high values.  Also shows the influence of high and low values on model
p3 <- ggplot()
p3 <- p3 + geom_line(aes(x = c(temp[5:19], 0.57), y = y_med_[[1]]))
p3 <- p3 + geom_line(aes(x = c(temp[5:19], 0.57), y = y_med_[[5]]), colour = "red")
p3 <- p3 + geom_line(aes(x = c(temp[5:19], 0.57), y = y_med_[[9]]), colour = "blue")
p3

# R-squared----
tab2 <- as.data.frame(rbind(
     R1_r2,
     M1_r2,
     M_unif_r2,
     RM1_r2,
     RM2_r2,
     RM3_r2,
     R2_r2,
     Mo_r2,
     Ro_r2
))

tab2 <- tibble::rownames_to_column(tab2)
tab2 <- cbind(index, tab2)
str(tab2)
names(tab2)[names(tab2)=="rowname"] <- "model"
names(tab2)[names(tab2)=="V1"] <- "R_sq"
names(tab2)[names(tab2)=="V2"] <- "SD"
tab2 <- tab2[order(tab2$R_sq), , drop = F]
tab2 <- tab2[order(tab2$R_sq, decreasing = T), , drop = F]
tab2
write.csv(tab2, paste0("Bayesian/", filepath_gen, "/Rsquared.csv"))

temp <- left_join(tab1, tab2, by = "index")
temp <- temp[c("model.x", "model.y","DIC", "dDIC", "R_sq", "SD")]

write.csv(temp, paste0("Bayesian/", filepath_gen, "/DIC_Rsq.csv"))


##############
##PI - proof that the PI for 2019 is wider than 2019 based on the arithmetic scale!
# from RM#
m <- matrix(c(2.707571, 3.720766,
              6.047425, 6.399664), nrow=2, ncol =2)

m[1,2]-m[1,1]
m[2,2]-m[2,1]

m1 <- exp(m)
m1[1,2]-m1[1,1]
m1[2,2]-m1[2,1]
# conclude that the PI is wider on arithmetic scale but not on the ln scale.
