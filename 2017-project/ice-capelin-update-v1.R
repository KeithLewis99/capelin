################################################################
#  Script written by Alejandro???? (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# The purpose of this file is to:
#1) Plot the optimization curves in comparing capelin abundance with timing of sea ice retreat.

# Notes:
# This does an outstanding job except for 2013 and 2014 - why?  What covariates might explain the difference?

rm(list=ls())
#setwd("D:/Keith/capelin/2017-project")

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)

## read in source code-----
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
## load Ale original data----
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

### load files generated in area-ice----
## max area
# look at files in "output-processing" and create a file pattern
data_path <- "output-processing"
files <- list.files(data_path)
pattern <- grep("capelin", files, value = T)
# read the datasets and join them to capelin_join
cape <- loadSubsetDatasets1(df = capelin_join, name = "capelin_m", pat = pattern, 6, var1 = "area", var2 = "minlats", nvar1 = "max_area", nvar2 = "minlats")

# seperate cape into time components
cape_1991 <- map(cape, ~filter(.x, year <= 1991))
cape_2017 <- map(cape, ~filter(.x, year > 1991))

## med area -----
pattern <- grep("sub1991-m+[^rsq].csv", files, value = T)
# read the datasets and join them to capelin_join

med_cape_1991 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

pattern <- grep("sub2017-m+[^rsq].csv", files, value = T)
med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

## D5 med area----
pattern <- grep("iceMedD5p2017-m.a", files, value = T)

d5med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

pattern <- grep("iceMedD5p1992-m.a", files, value = T)

d5med_cape_1992 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")


# merge these lists by data-----
med_cape_all <- map2(med_cape_1991, med_cape_2017, bind_rows)

d5med_cape_all <- map2(d5med_cape_1992, d5med_cape_2017, bind_rows)

########################################################################
## EXPLORATORY PLOTS------
## plots - capelin v. area
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")

titlenames <- c("m1", "m2", "m3", "m4", "m5", "m6")
for(i in 1:length(titlenames)){
     mm <- capelinAreaPlot(ls1=cape, ls2=med_cape_all, ls3=d5med_cape_all, i, titlenames=titlenames)
     ggsave(mm, filename = paste0("figs/capelinArea/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

## plots - lncapelin v. maxarea
titlenames <- c("m1-max", "m2-max", "m3-max", "m4-max", "m5-max", "m6-max")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(cape, cape_1991, cape_2017, "max_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

## plots - lncapelin v. medarea
titlenames <- c("m1-med", "m2-med", "m3-med", "m4-med", "m5-med", "m6-med")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(med_cape_all, med_cape_1991, med_cape_2017, "med_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

## plots - lncapelin v. d5med-area
titlenames <- c("m1-d5med", "m2-d5med", "m3-d5med", "m4-d5med", "m5-d5med", "m6-d5med")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(d5med_cape_all, d5med_cape_1992, d5med_cape_2017, "d5med_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}


########################################################################
## OPTIMIZATION PLOTS
## common features to all plots------
# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

yearInt <- seq(1982, 2014, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
#################################
## Ale's original-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
capelin_ale$myear <- capelin_ale$year
capelin_ale$myear[45:46] <- c(1960,1961)

optim_Ale <- calcFit(capelin_ale)

optimGraphs(optim_Ale$df, optim_Ale$regime1, optim_Ale$regime2, yearInt, lnbiomassInt, "max-values:Ale")
ggsave("figs/optimization/opt-maxAle.pdf", width=10, height=8, units="in")

## Optimization for max area----
titlenames <- c("opt_m1max", "opt_m2max", "opt_m3max", "opt_m4max", "opt_m5max", "opt_m6max")
optim_max <- calcFit_all(cape, titlenames)
str(optim_max, max.level = 3)
optim_max$optim_ls$opt_m1max$cdf

#create and save graphs
for(i in 1:length(optim_max$optim_ls)){
     df1 <- as.data.frame(optim_max$optim_ls[[i]]$df)
     df2 <- as.data.frame(optim_max$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(optim_max$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i])
     ggsave(mm, filename = paste0("figs/optimization/max", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

# compare Sums of Squares
# $value reports the minimized residual sums of squares or minimized Likelihood
# http://www.magesblog.com/2013/03/how-to-use-optim-in-r.html
# http://plantecology.syr.edu/fridley/bio793/likelihood.html - shows how to do AIC
for(i in 1:length(optim_max$optim_ls)){
     print(optim_max$optim_ls[[i]]$cdf$value)
}

for(i in 1:length(optim_max$optim_ls)){
     print(optim_max$optim_ls[[i]]$cdf$convergence)
}
#####################################################################################################
#graphics.off()
#####################################################################################################
par(mfrow=c(1,1))

capelin$logdiff <- NA
for(i in 2:nrow(capelin)){
  capelin$logdiff[i] <- (capelin$logcapelin[i]-capelin$logcapelin[i-1])*100
}

# new version of figure below
#png(filename = "ice-capelin-update2.png",width = 1000, height = 1000, units = "px", pointsize = 20, bg = "white", res = NA, family = "", restoreConsole = TRUE, type = c("cairo-png")) 

ggplot(capelin, aes(x = year, y = logdiff)) + 
  geom_point(pch = 16, size = 3) + 
  geom_hline(yintercept = 100, lty = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -100, lty = 3) +
  xlab("Year") +
  ylab("L%") + 
  #ylim(0,9) +
  theme_bw()

#plot of logdiff in capelin from previous year    
   with(capelin, plot(year,logdiff,xlab='year',ylab='L%',pch=16))
   abline(h=100,lty=3)
   abline(h=0)
   abline(h=-100,lty=3)
   
#dev.off()   
   
#####################################################################################################   
   
   a <- 3
   b <- 5
   g <- 2
   
   test <- a*capelin$tice
   (1-capelin$tice/b)

   ## Ale's original changed from tice to area-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present
   str(cape$capelin_m1)
   CapelinDomeFit_area <- optim(par = c(2,400000,0.6),
                           dataf = cape$capelin_m1[which(cape$capelin_m1$logcapelin!='NA'),
                                           c('year','max_area','logcapelin')],
                           fn = SSQCapelinDome, method=c("BFGS"))
   
   # before 2010
   CapelinDomeFitOld_area <- optim(par=c(2,400000,0.6),
                              dataf = cape$capelin_m1[which(cape$capelin_m1$year < 2011 & 
                                                                 cape$capelin_m1$logcapelin!='NA'),
                                              c('year','max_area','logcapelin')],
                              fn = SSQCapelinDome, 
                              method=c("L-BFGS-B"), control=list(maxit=50000))
   
   ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
   cape$capelin_m1$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit_area$par), dataf = cape$capelin_m1[,c('year','max_area')])
   
   cape$capelin_m1$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld_area$par), dataf = cape$capelin_m1[,c('year','max_area')])
   
   
   
   # attach the optimization curves of capelin abundance to ice data
   xarea <- expand.grid(year = c(1990,2000), area = seq(from = 50000, to = 400000, by = 500))
   xarea <- xarea[order(xarea$area),]
   xarea$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld_area$par),dataf = xarea)
   xarea$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit_area$par),dataf = xarea)
   #not sure what these are for but used in plots below but creates a data set where all values of year are the same????
   
   # make optimization graphs by year and in comparison to ice
   # plot of capelin biomass v. year with ice models  
   # plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for
   
   # set values
   regime1 <- xarea[which(xarea$year == 1990),]
   regime2 <- xarea[which(xarea$year == 2000),]
   
   yearInt <- seq(1982, 2018, by=4)
   lnbiomassInt <- seq(0, 10, by=2)
   biomassInt <- seq(0, 8500)
   labtice <- expression(paste(italic(t[area]), '(day of year)')) # label for figures
   
   windows()
ggplot() +
     geom_line(data = regime1, aes(x = area, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
     geom_line(data = regime2, aes(x = area, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
     geom_line(data = regime1, aes(x = area, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
     geom_line(data = regime2, aes(x = area, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
     geom_point(data = subset(cape$capelin_m1, year < 1991), aes(x = max_area, y = logcapelin), shape=2, size=3) +
     geom_point(data = subset(cape$capelin_m1, year > 1991), aes(x = max_area, y = logcapelin), shape=15, size=3) + 
     geom_errorbar(data = subset(cape$capelin_m1, year < 1991), aes(x = max_area, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
     geom_errorbar(data = subset(cape$capelin_m1, year > 1991), aes(x = max_area, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
     xlab("maxarea") +
     ylab("ln (Capelin biomass (ktons))") + 
     ylim(0,9) +
     theme_bw()   
   
   
   # AIC for models
   2* CapelinDomeFit_area$value+2*length(CapelinDomeFit_area$par)
   
   p1 <- ggplot(data = capelin, aes(x = maxarea, y = logcapelin)) + 
     geom_point() + 
     geom_smooth(method=lm)
   
   summary(lm(logcapelin ~ maxarea, data=capelin))
   summary(lm(logcapelin ~ maxarea, data=cape$capelin_m1))
   summary(lm(logcapelin ~ maxarea, data=cape$capelin_m2))
   summary(lm(logcapelin ~ maxarea, data=cape$capelin_m3))
   summary(lm(logcapelin ~ maxarea, data=cape$capelin_m4))
################################################################################