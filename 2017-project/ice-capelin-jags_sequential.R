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

rm(list=ls())

## read in source code-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")


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

#cap[26,11:13]
#temp <- subset(cap, year>1998)
#x <- (temp$ln_biomass_med[4] - temp$ln_biomass_med[12])/2


#cap$ln_biomass_med[26] <- cap$ln_biomass_med[26] + x
#cap$ln_bm_lci[26] <- log(cap$bm_lci)[26] + x
#cap$ln_bm_uci[26] <- log(cap$bm_uci)[26] + x

## ice (1969-2017)----
#source: ice-chart-processing-data-v3.R
ice <- read_csv('output-processing/capelin-m1.csv')
glimpse(ice)

## condition (1995-2015)-----
# source: capelin-condition.R
cond <- read_csv('data/condition_ag1_out.csv')
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
#View(condResids)

## larval data (2003-2017)----
# source "capelin_age_disaggregate_abundance.xlsx":Larval indices
larvae <- read_csv('data/capelin_larval_indices.csv')
glimpse(larvae)
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
capelin_join <- left_join(capelin_join, cond, by = "year")
capelin_join <- left_join(capelin_join, condResids, by = "year")
capelin_join <- left_join(capelin_join, larvae, by = "year")
capelin_join <- left_join(capelin_join, ps_tot, by = "year")

df <- capelin_join
glimpse(df)
#df[c("biomass_med", "ln_biomass_med", "tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")]
df[c("biomass_med", "ln_biomass_med", "tice", "meanCond_lag", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")]



# normalize the data set----
#df1 <- subset(df, year>1995)
df1 <- subset(df, year>1998)
glimpse(df1)

View(df1)
#cols <- c("tice", "meanCond_lag", "surface_tows_lag2", "ps_meanTot_lag2")

cols <- c("tice", "meanCond_lag", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")
df1 %<>%
     mutate_each_(funs(scale),cols)
write_csv(df1, "Bayesian/all_data.csv")

cols1 <- c("surface_tows", "ps_meanTot")
df_norm <- df1 %<>%
     mutate_each_(funs(scale),cols1)
write_csv(df_norm, "Bayesian/all_data_norm.csv")

# no tice
cols <- c("meanCond_lag", "resids_lag", "surface_tows_lag2", "ps_meanTot_lag2")
df1 %<>%
     mutate_each_(funs(scale),cols)
glimpse(df1)
df1$tice <- df1$tice/100
write_csv(df1, "Bayesian/all_data.csv")
# shrink df - not normal but small and easy to interpret parameters


# names for pairs-plot
df_name <- df1
names(df_name)[names(df_name) == "ln_biomass_med"] <- "ln capelin biomass"
names(df_name)[names(df_name) == "resids_lag"] <- "condition l1"
names(df_name)[names(df_name) == "surface_tows_lag2"] <- "larval abundance l2"
names(df_name)[names(df_name) == "ps_meanTot_lag2"] <- "zooplankton abun l2"

str(df_name)
df_name <- subset(df_name, year >= 2003)
# pairs-plot: relationships and correlations among RV and EV----
pdf("Bayesian/pairs_resids.pdf")
pairs.panels(df_name[c("ln capelin biomass", "tice", "condition l1", "larval abundance l2", "zooplankton abun l2")], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses,
             cex.labels = 1,
             cex.cor = 2
)
dev.off()

# calculator for normalization----
# Buren et al showed that the max of tice was ~ 94.  Need to plug this value into the prior but normalize it and do the same for variance.

test <- subset(df, year>1998)
gmean <- (94-mean(test$tice))/sd(test$tice, na.rm = T)


# variance - want variance that is large but was limited to 365 days - a larger variance for tice makes no sense.
var(test$tice)
gvar <- (241-mean(test$tice))/sd(test$tice, na.rm = T)

gmean^2/gvar # shape 

gvar/gmean # scale

#prior for beta
gmean <- (70-mean(test$tice))/sd(test$tice, na.rm = T)

(43-mean(test$tice, na.rm = T))/sd(test$tice, na.rm = T)

#########Bayesian models#############################
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

#subset data
df3 <- subset(df1, year>2002)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   PS=c(df3$ps_meanTot_lag2, c(-0.024866638, 0
)), #see df_norm
                   N = nrow(df3) + num_forecasts)

run_recruit <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'PRes', 'expY', 'D', 'log_lik'),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/recruitment_1/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/recruitment_1/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/recruitment_1/fit_obs.pdf')
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf('Bayesian/recruitment_1/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/recruitment_1/pi.csv")


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

#PLOT credible and prediction intervals
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/recruitment_1/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*pseudocal_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8)
ggsave("Bayesian/recruitment_1/credInt.png", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/recruitment_1/params.csv")


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

ggsave("Bayesian/recruitment_1/priorpost.pdf", width=10, height=8, units="in")


## Mortality----
#"alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]"
m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#mortaliyt
mu[i] <- alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]
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
beta ~ dgamma(5, 1/3) 
gamma ~ dgamma(2.8, 1)
delta ~ dnorm(0, 10^-2) 
sigma ~ dunif(0, 10) 
}'

# alpha and delta are based on normalized means and large uniformed priors
#beta based on Bolker pg 132 - Fig4.13 - trying for an uniformative prior 
#delta: see gammaDist_mode_shape_rate_calculator.R - 
# mode = 1.86 - based on Ale's value of 93 for beta*2 -> 186
# sd = 1
#shape(a) is mean^2/var; scale(s) equal Var/mean - but we use rate so mean/Var
# sigma: uninformative for condition
#df2 <- subset(df1, year>2002)
df2 <- df1
x <- as.Date('2018-02-26') # this is a minimum - I looked at the ice maps on 2018-03-02 and the ice is comming!
lubridate::yday(x)

x <- as.Date('2018-03-17')
lubridate::yday(x)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(0.5, 0.788)), #from capelin_larval_indices 
                   CO=c(df2$resids_lag, c(0, 0)), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mortality <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', 'log_lik'),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/mortality_1/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/mortality_1/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/mortality_1/fit_obs.pdf')
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
pdf('Bayesian/mortality_1/resid_covar.pdf')
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
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], "Bayesian/mortality_1/pi.csv")

## M-Results----

dic_M1 <- dic.samples(run_mortality$model, n.iter=1000, type="pD")
str(x)
dic_M1sum <- sum(dic_M1$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. Also, is this matched up properly - i haven't used PanelNames
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/mortality_1/posteriors.pdf", width=10, height=8, units="in")


# plot the credible and prediction intervals
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
#txt <- "log(capelin) = alpha +  beta*tice*(1-(tice/gamma)) \n + delta*Condition(resids)"

plotCredInt(df2, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2008, y = 7)

ggsave("Bayesian/mortality_1/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/mortality_1/params.csv")

# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)
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
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)
mean(gamma$df)/2

#delta
priormean <- 0
priorsd <- 20
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2)

ggsave("Bayesian/mortality_1/priorpost.pdf", width=10, height=8, units="in")




## M-(uniform)----
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
#df2 <- subset(df1, year>2002)
df2 <- df1
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(0.5, 0.788)), #from capelin_larval_indices 
                   CO=c(df2$resids_lag, c(0, 0)), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mort_unif <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/mortality_2/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/mortality_2/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/mortality_2/fit_obs.pdf')
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
pdf('Bayesian/mortality_2/resid_covar.pdf')
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
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], "Bayesian/mortality_2/pi.csv")

## M(uniform)-Results----
# Get the log likelihood
dic_M_unif <- dic.samples(run_mort_unif$model, n.iter=1000, type="pD")
dic_M_unifsum <- sum(dic_M_unif$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/mortality_2/posteriors.pdf", width=10, height=8, units="in")

txt <- "log(caplein) = alpha +  beta*tice*(1-(tice/gamma(unif))) \n + delta*Condition(resids)"

plotCredInt(df2, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2006, y = 8)
ggsave("Bayesian/mortality_2/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/mortality_2/params.csv")

# R-squared - see Gelmen paper

# plot posterior against expected distribution - BUT DO THIS ONLY FOR INFORMATIVE PRIORS
# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)
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
prior <- runif(n = 10000, 0, 100)
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma
prior <- runif(n = 10000, 0, 10)
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

ggsave("Bayesian/mortality_2/priorpost.pdf", width=10, height=8, units="in")



## RM_1----
#"mu[i] <- alpha + beta*ST[i] + gamma*CO[i]"

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
#df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   CO=c(df3$resids_lag, c(0, 0
)), #see df_norm
                   N = nrow(df3) + num_forecasts)

run_RM1 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm_1/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm_1/auto_corr.pdf", width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/rm_1/fit_obs.pdf')
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
pdf('Bayesian/rm_1/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm_1/pi.csv")

## RM_1-Results----
dic_RM1 <- dic.samples(run_RM1$model, n.iter=1000, type="pD")
dic_RM1sum <- sum(dic_RM1$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm_1/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2 \n + gamma*resids_adj"
#PLOT credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8)
ggsave("Bayesian/rm_1/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/rm_1/params.csv")


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

ggsave("Bayesian/rm_1/priorpost.pdf", width=10, height=8, units="in")




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
#delta ~ dunif(0, 3) 
sigma ~ dunif(0, 10) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#gamma: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
df3 <- subset(df1, year>2002)
num_forecasts = 2 # 2 extra years

model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   TI=c(df3$tice, c(0.5, 0.788)), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_RM2 <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm_2/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm_2/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/rm_2/fit_obs.pdf')
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
pdf('Bayesian/rm_2/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm_2/pi.csv")


## RM_2-Results----
dic_RM2 <- dic.samples(run_RM2$model, n.iter=1000, type="pD")
dic_RM2sum <- sum(dic_RM2$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm_2/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta))"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2010, y = 7.75)
ggsave("Bayesian/rm_2/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/rm_2/params.csv")


# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)
delta <- posterior_fig(out$sims.list$delta)


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
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
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
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
plot(density(prior))
plot(density(prior), xlim = c(0,2))

limits <- c(min(delta$df)-0.3, max(delta$df) + 0.3)
x_label <- "delta"
bin_1 <- mean(delta$df)/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2)

ggsave("Bayesian/rm_2/priorpost.pdf", width=10, height=8, units="in")


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
df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   TI=c(df3$tice, c(0.5, 0.788)), #made up - need new data
                   CO=c(df3$resids_lag, c(0, 0)),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm_3/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm_3/auto_corr.pdf", width=10, height=8, units="in")


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
pi_df3 <- apply(y_new,2,'quantile', c(0.05, 0.95))
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm_3/pi.csv")


## RM_3-Results----
dic_RM3 <- dic.samples(run_RM3$model, n.iter=1000, type="pD")
dic_RM3sum <- sum(dic_RM3$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm_3/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta)) + epsilon*CO"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2010, y = 8)
ggsave("Bayesian/rm_3/credInt.png", width=10, height=8, units="in")

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
df3 <- subset(df2, year>2002)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices 
                   PS=c(df3$ps_meanTot_lag2, c(-0.024866638, 0)), #made up - need new data
                   N = nrow(df3) + num_forecasts)

run_recruit_prior <- jags(data=model_data,
                    parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/recruitment_2/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/recruitment_2/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/recruitment_2/fit_obs.pdf')
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf('Bayesian/recruitment_2/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/recruitment_2/pi.csv")


## R1-inform_prior-Results----
# Get the log likelihood
dic_R2 <- dic.samples(run_recruit_prior$model, n.iter=1000, type="pD")
str(x)
dic_R2sum <- sum(dic_R2$deviance)


#PLOT credible and prediction intervals
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/recruitment_2/posteriors.pdf", width=10, height=8, units="in")

txt <- "log(caplein) = alpha +  beta(prior)*surface_tows_lag2 \n + gamma*pseudocal_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8)
ggsave("Bayesian/recruitment_2/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/recruitment_2/params.csv")

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

ggsave("Bayesian/recruitment_2/priorpost.pdf", width=10, height=8, units="in")


## Mortality_COalt----
# use the Paul/Ale measure of condition
#"alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]"

m.mortality = '
model {
# 1. Likelihood
for (i in 1:N) {
#mortaliyt
mu[i] <- alpha + beta*TI[i]*(1-TI[i]/gamma) + delta*CO[i]
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
alpha ~ dnorm(0, 100^-2) 
beta ~ dgamma(2, 1/3) 
gamma ~ dgamma(8.1, 1/11.11) 
delta ~ dnorm(0, 100^-2) 
sigma ~ dunif(0, 100) 
}'


#gamma based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#beta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
df2 <- df1
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(-2.190133, 0)), #from capelin_larval_indices 
                   CO=c(df2$meanCond_lag, c(0, 0)), #made up - need new data
                   N = nrow(df2) + num_forecasts)

run_mortality <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', 'log_lik'),
                      model.file = textConnection(m.mortality))

#Do we need to remove the burn-in?  How much?  UPDATE the chain? Figure out and apply to all

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
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/mortality_3/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/mortality_3/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/mortality_3/fit_obs.pdf')
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
pdf('Bayesian/mortality_3/resid_covar.pdf')
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
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
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], "Bayesian/mortality_3/pi.csv")

## M-Results----
# Get the log likelihood
log_lik_m1 = run_mortality$BUGSoutput$sims.list$log_lik
w_m1 <- waic(log_lik_m1)
w_m1 <- loo(log_lik_m1)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. Also, is this matched up properly - i haven't used PanelNames
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/mortality_3/posteriors.pdf", width=10, height=8, units="in")


# plot the credible and prediction intervals
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
txt <- "log(caplein) = alpha +  beta*tice*(1-(tice/gamma)) \n + delta*Condition(index)"

plotCredInt(df2, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2006, y = 8)
ggsave("Bayesian/mortality_3/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/mortality_3/params.csv")

# R-squared - see Gelmen paper

# plot posterior against expected distribution - BUT DO THIS ONLY FOR INFORMATIVE PRIORS
beta_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$beta)
dist <- as.data.frame(rgamma(1e+05, 2, 1/3))
names(dist)[names(dist) == "rgamma(1e+05, 2, 1/3)"] <- "v2" 
x <- range(beta_post)
p1 <- priorPosterior(beta_post, dist, x)

gamma_post <- as.data.frame(run_mortality$BUGSoutput$sims.list$gamma)
dist <- as.data.frame(rgamma(1000000, 8.1, 1/11.11))
names(dist)[names(dist) == "rgamma(1e+06, 8.1, 1/11.11)"] <- "v2" 
x <- range(gamma_post)
p2 <- priorPosterior(gamma_post, dist, x)

cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol=2)
ggsave("Bayesian/mortality_3/priorPost.pdf", width=10, height=8, units="in")


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
df2 <- df1
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df2$ln_biomass_med, rep(NA, num_forecasts)), 
                   TI=c(df2$tice, c(0.5, 0.789)), #from capelin_larval_indices 
                   N = nrow(df2) + num_forecasts)

run_mortality_null <- jags(data=model_data,
                      parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'beta', 'gamma', 'Fit', 'FitNew', 'PRes', 'expY', 'D', 'log_lik'),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/mortality_0/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/mortality_0/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/mortality_0/fit_obs.pdf')
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
pdf('Bayesian/mortality_0/resid_covar.pdf')
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
write.csv(pi_df2[, (ncol(pi_df2)-1):ncol(pi_df2)], "Bayesian/mortality_0/pi.csv")

## M-Null_Results----
dic_Mo <- dic.samples(run_mortality_null$model, n.iter=1000, type="pD")
dic_Mosum <- sum(dic_Mo$deviance)

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. Also, is this matched up properly - i haven't used PanelNames
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/mortality_0/posteriors.pdf", width=10, height=8, units="in")

# plot the credible and prediction intervals
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag2)")
#text(2004,7, "log(caplein) = delta +  Alpha*tice*(1-(tice/Beta)) \n + Gamma*Condition(ag1)")
#txt <- "log(caplein) = beta*tice*(1-(tice/gamma))"

plotCredInt(df2, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df2, pi_df2, 
            model = txt, x = 2006, y = 8)
ggsave("Bayesian/mortality_0/credInt.pdf", width=10, height=8, units="in")
ggsave("Bayesian/mortality_0/credInt.png", width=11, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/mortality_0/params.csv")

# plot posterior against expected distribution
# gamma dist
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig(out$sims.list$gamma)

#beta
priormean <- 5
priorsd <- 1/3
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(beta$df)-0.3, max(beta$df) + 0.3)
x_label <- "beta"
bin_1 <- mean(beta$df)/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)


#gamma
priormean <- 2.8
priorsd <- 1
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)
limits <- c(min(gamma$df)-0.3, max(gamma$df) + 0.3)
x_label <- "gamma"
bin_1 <- mean(gamma$df)/100
mean(gamma$df)/2

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p2, p3, ncol=2)

ggsave("Bayesian/mortality_0/priorpost.pdf", width=10, height=8, units="in")



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

#subset data
#df2 <- df1
df3 <- subset(df2, year>2002)

num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   N = nrow(df3) + num_forecasts)

run_recruit_null <- jags(data=model_data,
                    parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'Fit', 'FitNew', 'PRes', 'expY', 'D', 'log_lik'),
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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/recruitment_0/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this could be a problem!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/recruitment_0/auto_corr.pdf", width=10, height=8, units="in")

# Model Validation (see Zuuer et al. 2013 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/recruitment_0/fit_obs.pdf')
par(mfrow = c(2,1), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf('Bayesian/recruitment_0/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/recruitment_0/pi.csv")



## R_NUll-Results----
dic_Ro <- dic.samples(run_recruit_null$model, n.iter=1000, type="pD")
dic_Rosum <- sum(dic_Ro$deviance)

#PLOT credible and prediction intervals
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/recruitment_0/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha +  beta*surface_tows_lag2"

# plot the credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2008, y = 8)
ggsave("Bayesian/recruitment_0/credInt.pdf", width=10, height=8, units="in")
ggsave("Bayesian/recruitment_0/credInt.png", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/recruitment_0/params.csv")

# plot posterior against expected distribution
# gamma dist
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)

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

ggsave("Bayesian/recruitment_0/priorpost.pdf", width=10, height=8, units="in")


## RM3_projection----
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])"

RM3_1p_med = '
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
df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
TI=c(df3$tice, c(0.788, 0.788)), #made up - need new data
CO=c(df3$resids_lag, c(0, 0)),
N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
model.file = textConnection(RM3_1p_med))

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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm3_1p_med/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm3_1p_med/auto_corr.pdf", width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/rm3_1p_med/fit_obs.pdf')
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
pdf('Bayesian/rm3_1p_med/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm3_1p_med/pi.csv")


## RM_3-_proj_Results----
# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm3_1p_med/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta)) + epsilon*CO"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
ylab = "ln(capelin)", 
y_line = y_med, ci_df3, pi_df3, 
model = txt, x = 2010, y = 8)
ggsave("Bayesian/rm3_1p_med/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/rm3_1p_med/params.csv")


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

ggsave("Bayesian/rm3_1p_med/priorPost.pdf", width=10, height=8, units="in")




## RM3_projection----
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta + epsilon*CO[i])"

RM3_1p_high = '
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
df3 <- subset(df2, year>2002)
num_forecasts = 2 # 2 extra years
model_data <- list(N2 = c(df3$ln_biomass_med, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, c(-1.1798, -0.5525)), #from capelin_larval_indices - see df_norm
                   TI=c(df3$tice, c(0.98, 0.98)), #made up - need new data
                   CO=c(df3$resids_lag, c(0, 0)),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(RM3_1p_high))

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
MyBUGSChains(out, vars)
ggsave(MyBUGSChains(out, vars), filename = "Bayesian/rm3_1p_high/chains.pdf", width=10, height=8, units="in")

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
ggsave(MyBUGSACF(out, vars), filename = "Bayesian/rm3_1p_high/auto_corr.pdf", width=10, height=8, units="in")


# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf('Bayesian/rm3_1p_high/fit_obs.pdf')
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
pdf('Bayesian/rm3_1p_high/resid_covar.pdf')
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
write.csv(pi_df3[, (ncol(pi_df3)-1):ncol(pi_df3)], "Bayesian/rm3_1p_high/pi.csv")


## RM_3-_proj_Results----
# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 
MyBUGSHist(out, vars)
ggsave(MyBUGSHist(out, vars), filename = "Bayesian/rm3_1p_high/posteriors.pdf", width=10, height=8, units="in")

#txt <- "log(caplein) = alpha + beta*Surface_tows_lag2 \n + gamma*tice*(1-(tice/delta)) + epsilon*CO"

#plot credible and prediction intervals
plotCredInt(df3, yaxis = "ln_biomass_med", 
            ylab = "ln(capelin)", 
            y_line = y_med, ci_df3, pi_df3, 
            model = txt, x = 2010, y = 8)
ggsave("Bayesian/rm3_1p_high/credInt.pdf", width=10, height=8, units="in")

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), "Bayesian/rm3_1p_high/params.csv")


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

ggsave("Bayesian/rm3_1p_high/priorPost.pdf", width=10, height=8, units="in")


## DIC table ----
dic_R1sum
dic_M1sum
dic_RM1sum
dic_RM1sum
dic_RM2sum
dic_RM2sum
dic_RM3sum
dic_RM3sum
dic_Rosum
dic_Mosum
tab1 <- as.data.frame(rbind(
     dic_R1sum,
      dic_M1sum,
      dic_RM1sum,
      dic_RM2sum,
      dic_RM3sum,
      dic_Rosum,
      dic_Mosum,
      dic_M_unifsum,
     dic_R2sum
     
      ))
str(tab1)
names(tab1)[names(tab1)=="V1"] <- "DIC"
tab1 <- tab1[order(tab1$DIC), , drop = F]
tab1$dDIC <- tab1$DIC-tab1$DIC[1]
tab1
