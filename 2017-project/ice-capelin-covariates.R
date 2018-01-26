#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-08, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to:
# 1) Use function optim to test additional covariates to improve the Tice model fit
# 2) focus on capelin recruitment and pseudocalanus abundance. As well as capelin condition

rm(list=ls())

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(plotly)

## read in source code-----
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")

## READ/MANIPULATE DATA----
## read original data (Ale)----
capelin <- read.csv('capelin-ice-2014.csv',header=T)
# get rid of extra columns in csv file
capelin <- capelin[c("year", "maxarea", "minlat", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
capelin_ale <- capelin
# from "acoustic_estimates_2017.xlsx" in capelin folder
# add capelin, years, and lncapelin
capelin[c(47, 48, 49), 5] <- matrix(c(662, NA, 158), ncol = 1) 
capelin[c(47, 48, 49), 1] <- matrix(c(2015, 2016, 2017), ncol = 1) 
capelin[c(47, 48, 49), 8] <- matrix(c(log(662), NA, log(158)), ncol = 1)

# data set for joining to m1-m6 datasets
capelin_join <- capelin[c("year", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]

## read in Age2 capelin----
#from "capelin_age_disaggregate_abundance.xlsx" worksheet:Age disagg acoustic index
capelinAbun <- read.csv('data/capelin_age_disaggregate_abundance.csv',header=T)
str(capelinAbun)
capelinAbun$age2_log <- log(capelinAbun$age2)
capelinAbun$age2_log10 <- log10(capelinAbun$age2)
capelinAbun$capelin_log10 <- log10(capelinAbun$capelin)
capelinAbun$Ln_adult_abun <- log(sum(capelinAbun$age3, capelinAbun$age4, capelinAbun$age5, capelinAbun$age6, na.rm=T) + capelinAbun$age2*capelinAbun$age2PerMat)
capelinAbun$adult_abun <- sum(capelinAbun$age3, capelinAbun$age4, capelinAbun$age5, capelinAbun$age6, na.rm=T) + capelinAbun$age2*capelinAbun$age2PerMat

#View(capelinAbun)

# read in capelin condition----
# from "capelin_condition_maturation.xlsx"
capelinCond <- read.csv('data/capelin_condition_maturation.csv',header=T)
str(capelinCond)
#View(capelinCond)
capelinCond$resids_adj <- lag(capelinCond$resids*1000, 1)
capelinCond[19:21,3] <- NA

## read in larval data----
# from "capelin_age_disaggregate_abundance.xlsx":Larval indices
larvae <- read.csv('data/capelin_larval_indices.csv',header=T)
str(larvae)
larvae$surface_tows_lag2 <- lag(larvae$surface_tows, 2)
larvae$surface_tows_lag1 <- lag(larvae$surface_tows, 1)

#normalize the surface_tows 
larvae$Nsurface_tows_lag2 <- (larvae$surface_tows_lag2 - mean(larvae$surface_tows_lag2, na.rm=T)) / sd(larvae$surface_tows_lag2, na.rm = T)

## read in Pseudocalanus data----
# from "Copy of Copy of PSEUSP27_1999_2016.xlsx"
pscal <- read_csv("data/pseudocal_1999_2016.csv")
str(pscal)

# filter out year and use just the "total" stage.  Calcuate mean and sd by year
ps_tot <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage == "c6") %>%
     group_by(year) %>%
     summarise(ps_meanTot = mean(density), ps_sdTot = sd(density))

# lag mean and sd and /1000 to scale to tice
ps_tot$ps_meanTot_lag1 <- lag(ps_tot$ps_meanTot, 1)/100
ps_tot$ps_sdTot_lag1 <- lag(ps_tot$ps_sdTot, 1)/100
ps_tot$ps_meanTot_lag2 <- lag(ps_tot$ps_meanTot, 2)/100
ps_tot$ps_sdTot_lag2 <- lag(ps_tot$ps_sdTot, 2)/100
range(ps_tot$ps_meanTot_lag2, na.rm = TRUE)
ps_tot$ps_meanLog_lag2 <- lag(log10(ps_tot$ps_meanTot), 2)

##join the above data sets----
capelin_join <- left_join(capelin_join, larvae, by = "year")
capelin_join <- left_join(capelin_join, ps_tot, by = "year")
capelin_join <- left_join(capelin_join, capelinAbun, by = "year")
capelin_join <- left_join(capelin_join, capelinCond, by = "year")
str(capelin_join)

### load files generated in area-ice----
## max area
# look at files in "output-processing" and create a file pattern
data_path <- "output-processing"
files <- list.files(data_path)
pattern <- grep("capelin", files, value = T)
# read the datasets and join them to capelin_join
cape <- loadSubsetDatasets1(df = capelin_join, name = "capelin_m", pat = pattern, 6, var1 = "area", var2 = "minlats", nvar1 = "max_area", nvar2 = "minlats")

# Divide area by 1000 and normalize Ntice
for(i in 1:length(cape)){
     cape[[i]]$max_area1000 <- cape[[i]]$max_area/1000
     cape[[i]]$Ntice <- ((cape[[i]]$tice - mean(cape[[i]]$tice))/sd(cape[[i]]$tice)) + 5
     cape[[i]]$logsurface_tows_lag2 <- log(cape[[i]]$surface_tows_lag2)
     cape[[i]]$log10surface_tows_lag2 <- log10(cape[[i]]$surface_tows_lag2)
     cape[[i]]$logtice <- log(cape[[i]]$tice)
     cape[[i]]$Htice <- cape[[i]]$tice*100
     cape[[i]]$Ssurface_tows_lag2 <- cape[[i]]$surface_tows_lag2/10
}
range(cape$capelin_m1$Ssurface_tows_lag2, na.rm = T)

# seperate cape into time components
cape_1991 <- map(cape, ~filter(.x, year <= 1991))
cape_2017 <- map(cape, ~filter(.x, year > 1991))
cape_2001 <- map(cape, ~filter(.x, year > 2002 & year < 2017))

write.csv(cape_2001$capelin_m1, file = "figs/covariates/capelin_covariates_2001.csv", row.names=F, na="")

## EXPLORATORY ANALYSIS----
# simplify dataset
sdf <- cape_2001$capelin_m1[c("year", "tice", "logcapelin", "age2", "age2_log", "age2_log10", "surface_tows", "surface_tows_lag2", "log10surface_tows_lag2", "ps_meanTot", "ps_meanTot_lag1", "ps_meanTot_lag2", "ps_meanLog_lag2", "resids_adj")]
#View(sdf)
str(cape_2001$capelin_m1)

# plot explanatory variables against response variable and each other - check for trends and colinearity

#patterns
windows()
plot(sdf$surface_tows_lag2, sdf$logcapelin)
plot(sdf$ps_meanLog_lag2, sdf$logcapelin)
plot(sdf$ps_meanLog_lag2, sdf$log10capelin)
plot(sdf$ps_meanTot_lag1, sdf$logcapelin)
plot(sdf$ps_meanTot_lag2, sdf$logcapelin)
plot(sdf$ps_meanTot_lag2, sdf$logcapelin)
plot(sdf$resids_adj, sdf$logcapelin)

# colinearity
plot(sdf$tice, sdf$surface_tows_lag2)
plot(sdf$tice, sdf$surface_tows_lag2)
plot(sdf$tice, sdf$ps_meanTot_lag2)
plot(sdf$surface_tows_lag2, sdf$ps_meanTot_lag2)
plot(sdf$resids_adj, sdf$tice)
plot(sdf$resids_adj, sdf$surface_tows_lag2)
plot(sdf$resids_adj, sdf$ps_meanTot_lag2)


# compare logcapelin (biomass) against age2 (abundance) - these are highly correlated
plot(sdf$logcapelin, sdf$age2_log)
summary(lm(logcapelin ~ age2_log, data=sdf))

ggplot(data=sdf) + geom_point(aes(x = age2_log, y = logcapelin)) +
     geom_text(aes(x = age2_log, y = logcapelin + 0.05, label=year, size=3)) 
     
ggplot(data=sdf) + geom_point(aes(x = age2, y = logcapelin)) +
     geom_text(aes(x = age2, y = logcapelin + 0.05, label=year, size=3)) 

# plot variables in more detail - main graphic
# capelin v STow
ggplot(data=cape_2001$capelin_m1) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=tice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3, colour = "red") +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=round(ps_meanTot_lag2, 0)), nudge_x = -25, size=3, colour = "blue")

# capelin v Pcalanus(lag1)
ggplot(data=cape_2001$capelin_m1) + geom_point(aes(x = ps_meanTot_lag1, y = logcapelin)) +
     geom_text(aes(x = ps_meanTot_lag1, y = logcapelin + 0.1, label=tice), size=3) + 
     geom_text(aes(x = ps_meanTot_lag1, y = logcapelin - 0.1, label=year), size=3, colour = "red")

# capelin v Pcalanus(lag2) - match with Murphy et al 2017
# this mimics Murphy et al 2017 perfectly (2001-2014)
plot(sdf$ps_meanLog_lag2, sdf$age2_log10)
summary(lm(age2_log10 ~ ps_meanLog_lag2, data = sdf ))

ggplot(data=cape_2001$capelin_m1) + geom_point(aes(x = ps_meanLog_lag2, y = age2_log10)) +
#     geom_text(aes(x = ps_meanTot_lag2, y = logcapelin + 0.1, label=tice), size=3) + 
     geom_text(aes(x = ps_meanLog_lag2, y = age2_log10 - 0.1, label=year-2), size=3)

# this is very close to Murphy (2000-2015)
plot(sdf$log10surface_tows_lag2, sdf$age2_log10)
ggplot(data=cape_2001$capelin_m1) + 
     geom_point(aes(x = log10surface_tows_lag2, y = age2_log10)) +
     geom_text(aes(x = log10surface_tows_lag2, y = age2_log10 - 0.1, label=year-2), size=3)



## OPTIMIZE
## new values for optim functions----
yearInt <- seq(2002, 2017, by=3)
yearLim <- c(2002, 2017)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)


titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")

##MaxTice1 - raw data, 2 parm, 1 var, dome----
MaxTice1a <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "tice", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                          x1_range = seq(0,150, 10),
                          x2_range = NULL,
                          x3_range = NULL, 
                          method = c("L-BFGS-B"))

rawTice <- optimSummary(MaxTice1a, titlenames = titlenames)
var(cape_2001$capelin_m1$logcapelin, na.rm = T)
x <- cape_2001$capelin_m1
SS <- sum((x$logcapelin - mean(x$logcapelin, na.rm=T))^2, na.rm=T)

rawTice[1,4]/SS
optimGraphs2_all(MaxTice1a, "tice", var2 = NULL, var2val = NULL, file_name = "rawTice", "no")

# Surface tows - really, a dome makes no sense here
MaxTice1b <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 300), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                          x1_range = seq(30, 500, 50),
                          x2_range = NULL,
                          x3_range = NULL)

scalST <- optimSummary(MaxTice1b, titlenames = titlenames)
optimGraphs2_all(MaxTice1b, "Ssurface_tows_lag2", var2 = NULL, var2val = NULL, "scalST", "yes")

# Pseudo calanus - dome makes less sense here
MaxTice1c <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "ps_meanTot_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                          x1_range = seq(10, 280, 20),
                          x2_range = NULL,
                          x3_range = NULL)

scalPS2 <- optimSummary(MaxTice1c, titlenames = titlenames)
optimGraphs2_all(MaxTice1c, "ps_meanTot_lag2", var2 = NULL, var2val = NULL, "scalPS2", "yes")

# Surface tows_linear model
MaxTice1d <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1",
                          x1_range = seq(30, 500, 50),
                          x2_range = NULL,
                          x3_range = NULL)

scalST_lin <- optimSummary(MaxTice1d, titlenames = titlenames)
optimGraphs2_all(MaxTice1d, "Ssurface_tows_lag2", var2 = NULL, var2val = NULL, "scalST_lin", "yes")

# Pseudo calanus_linear model
MaxTice1e <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "ps_meanTot_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1",
                          x1_range = seq(10, 280, 20),
                          x2_range = NULL,
                          x3_range = NULL)
scalPS2_lin <- optimSummary(MaxTice1e, titlenames = titlenames)
optimGraphs2_all(MaxTice1e, "ps_meanTot_lag2", var2 = NULL, var2val = NULL, "scalPS2_lin", "yes")

# Pseudo calanus_Holling II
MaxTice1f <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "ps_meanTot_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1/(1 + (tmp1*Beta*Alpha))",
                          x1_range = seq(10, 280, 20),
                          x2_range = NULL,
                          x3_range = NULL,
                          method = c("L-BFGS-B"))

scalPS2_H2 <- optimSummary(MaxTice1f, titlenames = titlenames)
optimGraphs2_all(MaxTice1f, "ps_meanTot_lag2", var2 = NULL, var2val = NULL, "scalPS2_H2", "yes")

# Capelin condition - 
MaxTice1g <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 5), 
                          var1 = "resids_adj", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1",
                          x1_range = seq(-60, 50, 20),
                          x2_range = NULL,
                          x3_range = NULL)
scalRA_lin <- optimSummary(MaxTice1g, titlenames = titlenames)
optimGraphs2_all(MaxTice1g, "resids_adj", var2 = NULL, var2val = NULL, "scalRA_lin", "no")


##MaxTice1h - age2_log: raw data, 2 parm, 1 var, dome
MaxTice1h <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200), 
                          var1 = "tice", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "age2_log",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                          x1_range = seq(0,150, 10),
                          x2_range = NULL,
                          x3_range = NULL, 
                          method = c("L-BFGS-B"))

rawTice_abun <- optimSummary(MaxTice1h, titlenames = titlenames)
optimGraphs2_all(MaxTice1h, "tice", var2 = NULL, var2val = NULL, file_name = "rawTice_abun", "yes")

# Pseudo calanus_log
MaxTice1j <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 4), 
                          var1 = "ps_meanLog_lag2", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "age2_log",
                          form1 = "Beta + Alpha*tmp1",
                          x1_range = seq(3.5, 5, 0.1),
                          x2_range = NULL,
                          x3_range = NULL)

logPS2 <- optimSummary(MaxTice1j, titlenames = titlenames)
optimGraphs2_all(MaxTice1j, "ps_meanLog_lag2", var2 = NULL, var2val = NULL, "logPS2", "yes")

# rawTice adult abun
MaxTice1k <- calcFit_all3(cape_2001, 
                         titlenames, 
                         par = c(1, 200), 
                         var1 = "tice", 
                         var2 = NULL,
                         var3 = NULL, 
                         rv = "Ln_adult_abun",
                         form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                         x1_range = seq(0, 150, 10),
                         x2_range = NULL,
                         x3_range = NULL)

rawTice_LNadult_abun <- optimSummary(MaxTice1k, titlenames = titlenames)
optimGraphs3_all(MaxTice1k, 
                 var = "tice", 
                 var2 = NULL, 
                 var2val = NULL, 
                 file_name = "rawTice_LNadult_abun",
                 saveGraph = "no", 
                 rv1 = "Ln_adult_abun",
                 rv2 = "adult_abun",
                 y_axis_lim_p1 = c(6,8), 
                 ylabel_p1 = "ln(Adult_abun(ktons))",
                 ylabel_p2 = "Adult_abun(ktons))", 
                 ci="no")


##MaxTice2a - normalized data, 2 parm, 1 var, dome----
#flattens out the line

MaxTice2a <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 7), 
                          var1 = "Ntice", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                          x1_range = seq(0, 10, 0.067),
                          x2_range = NULL,
                          x3_range = NULL)

normTice <- optimSummary(MaxTice2a, titlenames = titlenames)
optimGraphs2_all(MaxTice2a, "Ntice", var2 = NULL, var2val = 0, "normTice", "yes")


### MaxTice4 log data, 2 parm, 1 var, dome----
MaxTice4a <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 1), 
                          var1 = "logtice", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                          x1_range = seq(3, 5, 0.1),
                          x2_range = NULL,
                          x3_range = NULL,
                          method = c("L-BFGS-B"),
                          lowerLim = "yes")

logTice <- optimSummary(MaxTice4a, titlenames = titlenames)
optimGraphs2_all(MaxTice4a, "logtice", var2 = NULL, var2val = 4.0, file_name = "logTice", "no")


### MaxTice6 - scaled data, 3 parm, 2 var, dome + line----
MaxTice6 <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = NULL)

scalTiceSTow <- optimSummary(MaxTice6, titlenames = titlenames)
var2val(MaxTice6$optim_ls$`MaxTice-m1`$regime2$Ssurface_tows_lag2)

optimGraphs3_all(MaxTice6, 
                 var = "tice", 
                 var2 = "Ssurface_tows_lag2", 
                 var2val = 230, 
                 file_name = "scalTiceSTow",
                 saveGraph = "no", 
                 rv1 = "logcapelin",#"Ln_adult_abun",
                 rv2 = "capelin",#"adult_abun",
                 y_axis_lim_p1 = c(4, 9), 
                 ylabel_p1 = "ln(Adult_abun(ktons))",
                 ylabel_p2 = "Adult_abun(ktons))", 
                 ci="no")

library(plotly)

t <- MaxTice6$optim_ls$`MaxTice-m1`$regime2
head(t)
p <- plot_ly(t, x = ~tice, y = ~Ssurface_tows_lag2, z = ~ExpectedLogBiomass, type = 'scatter3d', mode = 'lines',
             opacity = 1, line = list(width = 6, reverscale = FALSE))


head(MaxTice6$optim_ls$`MaxTice-m1`$regime2)
### MaxTice6a
# psuedocalanus
#scaled data, 3 parm, 2 var, dome + line
MaxTice6a <- calcFit_all3(cape_2001, 
                         titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", 
                         var2 = "ps_meanTot_lag2",
                         var3 = NULL, 
                         rv = "logcapelin",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(10, 280, 20),
                         x3_range = NULL)

scalTicePS2 <- optimSummary(MaxTice6a, titlenames = titlenames)
var2val(MaxTice6a$optim_ls$`MaxTice-m1`$regime2$ps_meanTot_lag1)
optimGraphs2_all(MaxTice6a, "tice", var2 = "ps_meanTot_lag1", var2val = 90, file_name = "scalTicePS2", saveGraph = "yes")

### MaxTice6b scaled data, 3 parm, 2 var, dome + line
# as above, this really makes no sense biologically
MaxTice6b <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = "ps_meanTot_lag2",
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1 + Gamma*tmp2",
                          x1_range = seq(30, 500, 50),
                          x2_range = seq(10, 280, 20),
                          x3_range = NULL)
scalSTowPS2_lin <- optimSummary(MaxTice6b, titlenames = titlenames)
var2val(MaxTice6a$optim_ls$`MaxTice-m1`$regime2$ps_meanTot_lag1)
optimGraphs2_all(MaxTice6b, "Ssurface_tows_lag2", var2 = "ps_meanTot_lag2", var2val = 90, file_name = "scalSTowPS2_lin", saveGraph = "yes")

### MaxTice6c - Holling IV: scaled data, 3 parm, 1 var 
# this is clunky code - there is only one variable being used but I need to put two in bc of the clunky code
MaxTice6c <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1^2/(Beta + Gamma*tmp1+tmp1^2)",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = NULL)

scalTice_H4 <- optimSummary(MaxTice6c, titlenames = titlenames)
var2val(MaxTice6c$optim_ls$`MaxTice-m1`$regime2$Ssurface_tows_lag2)
optimGraphs2_all(MaxTice6c, "tice", var2 = NULL, var2val = 280, file_name = "scalTice_H4", saveGraph = "no")


### MaxTice6d - Holling IV: scaled data, 4 parm, 2 var
# as above - clunky code
MaxTice6d <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = "ps_meanTot_lag2", 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1^2/(Beta + Gamma*tmp1+tmp1^2)+ Delta*tmp2",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = seq(10, 280, 20))

scalTiceST_H4 <- optimSummary(MaxTice6d, titlenames = titlenames)
var2val(MaxTice6c$optim_ls$`MaxTice-m1`$regime2$Ssurface_tows_lag2)
optimGraphs2_all(MaxTice6d, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceST_H4", saveGraph = "yes")

### MaxTice6e - dome + Holling II: scaled data, 4 parm, 2 var
# as above - clunky code
MaxTice6e <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", 
                          var2 = "ps_meanTot_lag2",
                          var3 = "Ssurface_tows_lag2", 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))+ Gamma*tmp2/(1 + (tmp2*Delta*Gamma))",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(10, 200, 20),
                          x3_range = seq(30, 500, 50))

scalTicePS2_domeH2 <- optimSummary(MaxTice6e, titlenames = titlenames)
var2val(MaxTice6e$optim_ls$`MaxTice-m1`$regime2$ps_meanTot_lag2)
optimGraphs2_all(MaxTice6e, "tice", var2 = "ps_meanTot_lag2", var2val = 90, file_name = "scalTicePS2_domeH2", saveGraph = "yes")

MaxTice6e$optim_ls$`MaxTice-m1`$regime2


MaxTice6f <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 1, 1), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = "ps_meanLog_lag2",
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1+ Gamma*tmp2",
                          x1_range = seq(30, 500, 50),
                          x2_range = seq(10, 280, 20),
                          x3_range = NULL)
scalSTowPS2_lin <- optimSummary(MaxTice6f, titlenames = titlenames)
var2val(MaxTice6a$optim_ls$`MaxTice-m1`$regime2$ps_meanTot_lag1)
optimGraphs2_all(MaxTice6f, "Ssurface_tows_lag2", var2 = "ps_meanTot_lag2", var2val = 90, file_name = "scalSTowPS2_lin", saveGraph = "no")

optimGraphs3_all(MaxTice6f, 
                 var = "Ssurface_tows_lag2", 
                 var2 = "ps_meanLog_lag2", 
                 var2val = 90, 
                 file_name = "scalSTowPS2_lin",
                 saveGraph = "no", 
                 rv1 = "logcapelin",#"Ln_adult_abun",
                 rv2 = "capelin",#"adult_abun",
                 y_axis_lim_p1 = c(4, 9), 
                 ylabel_p1 = "ln(capelin(ktons))",
                 ylabel_p2 = "capelin(ktons))", 
                 ci="yes")

t1 <- MaxTice6f$optim_ls$`MaxTice-m1`$regime2
head(t1)
p <- plot_ly(t1, x = ~tice, y = ~Ssurface_tows_lag2, z = ~ExpectedLogBiomass, type = 'scatter3d', mode = 'lines',
             opacity = 1, line = list(width = 6, reverscale = FALSE))
t2 <- MaxTice6f$optim_ls$`MaxTice-m1`$df
q <- plot_ly(t2, x = ~ps_meanTot_lag2, y = ~Ssurface_tows_lag2, z = ~logcapelin, type = 'scatter3d') %>% add_markers()

p <- p + plot_ly(t2, x = ~ps_meanTot_lag2, y = ~Ssurface_tows_lag2, z = ~logcapelin, type = 'scatter3d') %>% add_markers()


MaxTice6g <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 1, 1), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = "resids_adj",
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1+ Gamma*tmp2",
                          x1_range = seq(30, 500, 50),
                          x2_range = seq(-60, 50, 20),
                          x3_range = NULL)
scalSTowRA_lin <- optimSummary(MaxTice6g, titlenames = titlenames)

optimGraphs3_all(MaxTice6g, 
                 var = "Ssurface_tows_lag2", 
                 var2 = "resids_adj", 
                 var2val = 90, 
                 file_name = "scalSTowRA_lin",
                 saveGraph = "no", 
                 rv1 = "logcapelin",#"Ln_adult_abun",
                 rv2 = "capelin",#"adult_abun",
                 y_axis_lim_p1 = c(4, 9), 
                 ylabel_p1 = "ln(capelin(ktons))",
                 ylabel_p2 = "capelin(ktons))", 
                 ci="yes")

MaxTice6h <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 1, 1, 1), 
                          var1 = "Ssurface_tows_lag2", 
                          var2 = "resids_adj",
                          var3 = "ps_meanLog_lag2", 
                          rv = "logcapelin",
                          form1 = "Beta + Alpha*tmp1+ Gamma*tmp2  + Delta*tmp3",
                          x1_range = seq(30, 500, 50),
                          x2_range = seq(-60, 50, 20),
                          x3_range = seq(10, 280, 20))
scalSTowRAPS_lin <- optimSummary(MaxTice6h, titlenames = titlenames)

optimGraphs3_all(MaxTice6h, 
                 var = "Ssurface_tows_lag2", 
                 var2 = "resids_adj", 
                 var2val = 90, 
                 file_name = "scalSTowRA_lin",
                 saveGraph = "no", 
                 rv1 = "logcapelin",#"Ln_adult_abun",
                 rv2 = "capelin",#"adult_abun",
                 y_axis_lim_p1 = c(3, 9), 
                 ylabel_p1 = "ln(capelin(ktons))",
                 ylabel_p2 = "capelin(ktons))", 
                 ci="yes")

MaxTice6i <- calcFit_all3(cape_2001, 
                         titlenames, 
                         par = c(1, 200, 1, 1), 
                         var1 = "tice", 
                         var2 = "Ssurface_tows_lag2",
                         var3 = NULL, 
                         rv = "logcapelin",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50),
                         x3_range = NULL)

scalTiceST_int <- optimSummary(MaxTice6i, titlenames = titlenames)

t1 <- MaxTice6i$optim_ls$`MaxTice-m1`$regime2
head(t1)
p <- plot_ly(t1, x = ~tice, y = ~Ssurface_tows_lag2, z = ~ExpectedLogBiomass, type = 'scatter3d', mode = 'lines',
             opacity = 1, line = list(width = 6, reverscale = FALSE))

### MaxTice8 scaled data, 4 parm, 3 var, dome + 2line----
MaxTice8 <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = "ps_meanTot_lag2", 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = seq(10, 280, 20))

scalTiceSTowPS2 <- optimSummary(MaxTice8, titlenames = titlenames)
optimGraphs2_all(MaxTice8, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowPS2", saveGraph = "yes")



#arrange(MaxTice8$optim_ls$`MaxTice-m1`$regime2, Ssurface_tows_lag2)
#arrange(filter(MaxTice8$optim_ls$`MaxTice-m1`$regime2, Ssurface_tows_lag2 == 230), tice)

# As above but with ps_lag2
MaxTice8a <- calcFit_all3(cape_2001, 
                         titlenames, 
                         par = c(1, 200, 1, 1), 
                         var1 = "tice", 
                         var2 = "Ssurface_tows_lag2",
                         var3 = "ps_meanTot_lag2", 
                         rv = "logcapelin",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50),
                         x3_range = seq(10, 280, 20))

scalTiceSTowPS_lag2 <- optimSummary(MaxTice8a, titlenames = titlenames)
optimGraphs2_all(MaxTice8a, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowPS_lag2", saveGraph = "yes")

# As above but with ps_lag2 and Holling II
# invalid model - not enough parameters
#MaxTice8b <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = "ps_meanTot_lag2", 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Alpha*tmp3/(1 + (tmp3*Delta*Beta))",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = seq(10, 200, 20))

#optimSummary(MaxTice8b, titlenames = titlenames)
#optimGraphs2_all(MaxTice8b, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowPS_l2", saveGraph = "no")

### MaxTice8d scaled data, 4 parm, 3 var, dome + 2line
# as above but with resids adjusted
range(cape_2001$capelin_m1$resids_adj, na.rm = T)
MaxTice8d <- calcFit_all3(cape_2001, 
                          titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", 
                          var2 = "Ssurface_tows_lag2",
                          var3 = "resids_adj", 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = seq(-60, 50, 20))

scalTiceSTowRA <- optimSummary(MaxTice8d, titlenames = titlenames)
optimGraphs2_all(MaxTice8d, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowRA", saveGraph = "no")


####
## Table of optimSummary results----
optimSummary_ls <- list(rawTice=rawTice, 
                        scalST = scalST, 
                        scalPS2 = scalPS2,
                        scalST_lin = scalST_lin,
                        scalPS2_lin = scalPS2_lin,
                        scalPS2_H2 = scalPS2_H2,
                        scalRA_lin = scalRA_lin,
                        rawTice_abun = rawTice_abun,
                        logPS2 = logPS2,
                        rawTice_adult_abun = rawTice_adult_abun,
                        #one var transform
                        normTice = normTice, 
                        logTice = logTice, 
                        #two vars
                        scalTiceSTow = scalTiceSTow,
                        scalTicePS2 = scalTicePS2,
                        scalSTowPS2_lin = scalSTowPS2_lin,
                        scalTice_H4 = scalTice_H4,
                        scalTiceST_H4 = scalTiceST_H4,
                        # three vars
                        scalTicePS2_domeH2 = scalTicePS2_domeH2, 
                        scalTiceSTowPS2 = scalTiceSTowPS2, 
                        scalTiceSTowPS_lag2 = scalTiceSTowPS_lag2, 
                        scalTiceSTowRA = scalTiceSTowRA)

#save as a list - pass to console and cut/paste into excel; then use space delimited
optimSummary_ls

## med area -----
# tells you almost nothing bc the ice runs from 1 Jan to 1 Jun and ice ends usually May.  Therefore, median values are always in 70s.
# read the datasets and join them to capelin_join
pattern <- grep("sub2017-m+[^rsq].csv", files, value = T)
med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(med_cape_2017)){
     med_cape_2017[[i]]$med_area1000 <- med_cape_2017[[i]]$med_area/1000
}

for(i in 1:length(med_cape_2017)){
     med_cape_2017[[i]]$Ntice <- ((med_cape_2017[[i]]$tice - mean(med_cape_2017[[i]]$tice))/sd(med_cape_2017[[i]]$tice)) + 5
     med_cape_2017[[i]]$logsurface_tows_lag2 <- log(med_cape_2017[[i]]$surface_tows_lag2)
     med_cape_2017[[i]]$logtice <- log(med_cape_2017[[i]]$tice)
     med_cape_2017[[i]]$Htice <- med_cape_2017[[i]]$tice*100
     med_cape_2017[[i]]$Ssurface_tows_lag2 <- med_cape_2017[[i]]$surface_tows_lag2/10
}

# seperate cape into time components
med_cape_2001 <- map(med_cape_2017, ~filter(.x, year > 2001))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- med_cape_2001$med_m1[c("year", "dtice", "logcapelin", "surface_tows_lag2")]

windows()
ggplot(data=med_cape_2001$med_m1) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=dtice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)

med_titlenames <- c("MedArea-m1", "MedArea-m2", "MedArea-m3", "MedArea-m4", "MedArea-m5", "MedArea-m6")

#optimize
var2val(med_cape_2001$med_m1$dtice)
MaxTice10 <- calcFit_all3(med_cape_2001, 
                         med_titlenames, 
                         par = c(1, 75, 1, 1), 
                         var1 = "dtice", 
                         var2 = "Ssurface_tows_lag2",
                         var3 = "ps_meanTot_lag1", 
                         rv = "logcapelin",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(70, 80, 0.5),
                         x3_range = seq(10, 280, 20),
                         lowerLim = "yes")

med_scalTiceSTowPS <- optimSummary(MaxTice10, titlenames = titlenames)
optimGraphs2_all(MaxTice10, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "med_scalTiceSTowPS", saveGraph = "yes")

med_cape_2001$med_m1

###
## D5 med area----
pattern <- grep("iceMedD5p2017-m.a", files, value = T)

d5med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(d5med_cape_2017)){
     d5med_cape_2017[[i]]$d5med_area1000 <- d5med_cape_2017[[i]]$d5med_area/1000
}


for(i in 1:length(d5med_cape_2017)){
     d5med_cape_2017[[i]]$Ntice <- ((d5med_cape_2017[[i]]$tice - mean(d5med_cape_2017[[i]]$tice))/sd(d5med_cape_2017[[i]]$tice)) + 5
     d5med_cape_2017[[i]]$logsurface_tows_lag2 <- log(d5med_cape_2017[[i]]$surface_tows_lag2)
     d5med_cape_2017[[i]]$logtice <- log(d5med_cape_2017[[i]]$tice)
     d5med_cape_2017[[i]]$Htice <- d5med_cape_2017[[i]]$tice*100
     d5med_cape_2017[[i]]$Ssurface_tows_lag2 <- d5med_cape_2017[[i]]$surface_tows_lag2/10
}

# seperate cape into time components
d5med_cape_2001 <- map(d5med_cape_2017, ~filter(.x, year > 2001))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- d5med_cape_2001$d5med_m1[c("year", "d5tice", "logcapelin", "surface_tows_lag2")]

windows()
ggplot(data=d5med_cape_2001$d5med_m1) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=d5tice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)



d5titlenames <- c("D5MedArea-m1", "D5MedArea-m2", "D5MedArea-m3", "D5MedArea-m4", "D5MedArea-m5", "D5MedArea-m6")

var2val(d5med_cape_2001$d5med_m1$d5tice)
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")

MaxTicetest <- calcFit_all3(d5med_cape_2001, 
                            d5titlenames, 
                          par = c(1, 200), 
                          var1 = "d5tice", 
                          var2 = NULL,
                          var3 = NULL, 
                          rv = "logcapelin",
                          form1 = "Alpha*tmp1*(1 - (tmp1/Beta))",
                          x1_range = seq(0,150, 10),
                          x2_range = NULL,
                          x3_range = NULL)

d5scalTiceSTowPS <- optimSummary(MaxTicetest, titlenames = d5titlenames)

# only works for one variable and convergence is crap - abandon!