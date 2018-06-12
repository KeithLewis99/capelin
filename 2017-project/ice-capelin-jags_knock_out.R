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


## data knock out----
filepath_gen <- "Bayesian/biomass_cond_ag1_2_DIC/" #"ag2_cond_ag1_DIC" "Bayesian/biomass_cond_ag1_DIC/"

dataset <- "biomass" #"age"
if(dataset == "biomass"){
     df <- read_csv("Bayesian/biomass_cond_ag1/all_data_a1.csv")
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

##ITERATIONS START HERE----
#df3: 2003-2017
df3$year
x <- 4 #1:15 1 = 2003,  8 = 2010, 15 = 2017: DO THIS FOR ALL OF THE YEARS INCLUDING 2006 AND 2016!!!
insert <- df3[x,] # this is the value of the year to insert into the graph
df3 <- df3[-x, ]

#df2: 1999-2017
y <- x+4
insert_year <- 2006 # value of year to align with credible interval
insert_row  <- x
folder <- "sensitivity"

st <- df2[c("year", "surface_tows_lag2")][y,]
ti <- df2[c("year", "tice")][y,]
co <- df2[c("year", "meanCond_lag")][y,]
st[,2]
ti[,2]
co[,2]

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
beta ~ dnorm(0, 10^-2) 
gamma ~ dgamma(5, 1/3) 
delta ~ dgamma(2.8, 1)
epsilon ~ dnorm(0, 20^-2) 
sigma ~ dunif(0, 100) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#delta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition


num_forecasts =  1# 1 extra years
# df3$ln_biomass_med
# df3$age2_log10
model_data <- list(N2 = c(var_x, rep(NA, num_forecasts)), 
                   ST=c(df3$surface_tows_lag2, st[,2]),
                   TI=c(df3$tice, ti[,2]), 
                   CO=c(df3$meanCond_lag, co[,2]),
                   N = nrow(df3) + num_forecasts)

run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM3))

#UPDATE WITH MORE BURN INS
run_RM3 <-update(run_RM3, n.iter = 300000, n.thin = 50, n.burnin = 100000)


##R/M3-plot 5----
# DIAGNOSTICS
out <- run_RM3$BUGSoutput 

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


## RM_3-Results----

#plot credible and prediction intervals
plotCredInt2(df3, insert, 
             yaxis = yaxis1, 
             ylab = ylab1, 
             y_line = y_med, 
             ci = ci_df3, 
             dpp = p_med,
             dpi = pi_df3, 
             insert_year = insert_year, type = NA)

ggsave(paste0(filepath_gen, folder, "/credInt", insert_year, ".png"), width=10, height=8, units="in")

# this creates a file that accumulates the predicted value of the knockout values
per_diff_file <- data.frame(year = integer(), 
                            per_diff = numeric(), # percent difference between predicted knockout value
                            p_med = numeric(), # median value of PI
                            pi025 = numeric(), # 2.5% value of PI
                            pi975 = numeric()) # 97.5% value of PI


if(dataset == "biomass"){
     per_diff <- ((df2$ln_biomass_med[y] - p_med[15])/df2$ln_biomass_med[y])*100
     #per_diff_file <- as.data.frame(cbind(insert_year, per_diff))
     per_diff_file[1,2] <- per_diff
} else if (dataset == "age"){
     per_diff <- ((df2$age2_log10[y] - p_med[15])/df2$age2_log10[y])*100
     #per_diff_file <- as.data.frame(cbind(insert_year, per_diff))
     per_diff_file[1,2] <- per_diff
     }

# this enters the above values into the dataframe "per_diff_fileo
# re-run this for every year to create the file
per_diff_file[1,1] <- insert_year
per_diff_file[1,3] <- p_med[15]
per_diff_file[1,4] <- pi_df3[1,15]
per_diff_file[1,5] <- pi_df3[2,15]

# save file as in next line, then hash tag and do the next three lines for all subsequent analyses
#write_csv(per_diff_file, paste0(filepath_gen, folder, "/perdiff_all.csv"))

temp <- read_csv(paste0(filepath_gen, folder, "/perdiff_all.csv"))

temp <- rbind(per_diff_file, temp)

write_csv(temp, paste0(filepath_gen, folder, "/perdiff_all.csv"))


## graph of values v. knockout----
# use final values from above - change "per_diff_file" in ggplot code to temp
# need to generalize this for biomass and age!

type <- "CI"

temp$data <- "PI"
df3$data <- "CI"

df4 <- df3[,c("year", "ln_biomass_med", "ln_bm_lci", "ln_bm_uci", "data")]
temp1 <- temp[, -2]
temp1[4, 2:4] <- NA
temp1[14, 2:4] <- NA
names(df4) <- names(temp1)
temp2 <- rbind(temp1, df4)
View(temp2)
dg <- 0.5
p <- ggplot(data = temp2, aes(x = year, y = p_med, group = data, color = data)) 
p <- p + geom_point(position = position_dodge(width = dg), size = 1.5)
#p <- p + geom_linerange(aes(ymin = (pi025), ymax = (pi975)), position = position_dodge(width = dg), size = 0.7)
p <- p + geom_errorbar(aes(ymin = (pi025), ymax = (pi975)), position = position_dodge(width = dg), size = 0.7)

p <- p + labs(x = 'Year', y = "ln(capelin biomass(ktons))")
p <- p + scale_color_manual(values = c("red", "black"))
p <- p + theme_bw(base_size = 25)
p <- p + theme(legend.title = element_blank())#, axis.text.x  = element_text(angle=90, vjust=0.5))
ggsave(paste0(filepath_gen, "/sensitivity/loo.pdf"), width=10, height=8, units="in")


###########
# can probably delete all of this below here
p <- p + geom_point(data = temp2, 
                    aes(y=p_med, x = year), size = 3)
p <- p + geom_linerange(data = temp2,
                       aes(ymax = pi975, ymin = pi025, x = year), width = 0.5)
p <- p + geom_point(data = df3, 
                    aes(y=ln_biomass_med, x = year), #age2_log10
                    colour = "red", size = 3,
                    position = position_nudge(x = -0.25))
if(!is.na(type)){
     p <- p + geom_linerange(data = df3, 
                            aes(x = year, ymin=ln_bm_lci, ymax=ln_bm_uci), 
                            width = 0.3, colour = "red", 
                            position = position_nudge(x = -0.25, y = 0))
} else if (is.na(type)) {
     
     p
}
p <- p + xlab("Year") + ylab("Estimate and prediction")
p <- p + theme_bw(base_size = 20) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
#p <- p + scale_x_discrete(breaks = c("HH", "HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL"), labels = c("HH", "HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL"))
p
ggsave(paste0("Bayesian/", filepath_gen, "/sensitivity/loo.pdf"), width=10, height=8, units="in")
