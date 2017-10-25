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
setwd("D:/Keith/capelin/2017-project")

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)

## read in data-----
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
## Ale original
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


# look at files in "output-processing" and create a file pattern
data_path <- "output-processing"
files <- list.files(data_path)
pattern <- grep("capelin", files, value = T)
# read the datasets and join them to capelin_join
cape <- loadSubsetDatasets1(df = capelin_join, name = "capelin_m", pat = pattern, 6, var1 = "area", var2 = "minlats", nvar1 = "max_area", nvar2 = "minlats")

###############################################################################
## common features to all plots------
# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

yearInt <- seq(1982, 2014, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
###############################################################################
## Ale's original-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
capelin_ale$myear <- capelin_ale$year
capelin_ale$myear[45:46] <- c(1960,1961)

optim_Ale <- calcFit(capelin_ale)

optimGraphs(optim_Ale$df, optim_Ale$regime1, optim_Ale$regime2, yearInt, lnbiomassInt, "max-values:Ale")
ggsave("figs/u1-maxAle.pdf")

## m1--------------
yearInt <- seq(1982, 2018, by=4)

optim_m1 <- calcFit(cape$capelin_m1)

#windows()
#note that multiple windows can be opened but cowplot does not seem to recognize the news windoes and plots to the old one.  if you want to see plots in multiple widows, use multiplot but cannot label
optimGraphs(optim_m1$df, optim_m1$regime1, optim_m1$regime2, yearInt, lnbiomassInt,  "max-values:m1")

ggsave("figs/u2-maxm1.pdf")

## m2--------------
#source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
optim_m2 <- calcFit(cape$capelin_m2)

optimGraphs(optim_m2$df, optim_m2$regime1, optim_m2$regime2, yearInt, lnbiomassInt,  "max-values:m2")

ggsave("figs/u3-maxm2.pdf")

## m3--------------
optim_m3 <- calcFit(cape$capelin_m3)

optimGraphs(optim_m3$df, optim_m3$regime1, optim_m3$regime2, yearInt, lnbiomassInt,  "max-values:m3")

ggsave("figs/u4-maxm3.pdf")

## m4--------------
optim_m4 <- calcFit(cape$capelin_m4)

optimGraphs(optim_m4$df, optim_m4$regime1, optim_m4$regime2, yearInt, lnbiomassInt,  "max-values:m4")

ggsave("figs/u5-maxm4.pdf")

## m5--------------
optim_m5 <- calcFit(cape$capelin_m5)

optimGraphs(optim_m5$df, optim_m5$regime1, optim_m5$regime2, yearInt, lnbiomassInt,  "max-values:m5")

ggsave("figs/u6-maxm5.pdf")

## m6--------------
optim_m6 <- calcFit(cape$capelin_m6)
windows()
optimGraphs(optim_m6$df, optim_m6$regime1, optim_m6$regime2, yearInt, lnbiomassInt,  "max-values:m6")

ggsave("figs/u7-maxm6.pdf")

# compare Sums of Squares
2*optim_Ale$cdf$value+2*length(optim_Ale$cdf$par)


optim_m1$cdf$value
optim_m2$cdf$value
optim_m3$cdf$value
optim_m4$cdf$value
optim_m5$cdf$value
optim_m6$cdf$value


########################################################################

# look at files in "output-processing" and create a file pattern
data_path <- "output-processing"
files <- list.files(data_path)
pattern <- grep("sub1991-m+[^rsq].csv", files, value = T)
# read the datasets and join them to capelin_join

med_cape_1991 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

pattern <- grep("sub2017-m+[^rsq].csv", files, value = T)
med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")


pattern <- grep("iceMedD5p2017-m.a", files, value = T)

d5med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

pattern <- grep("iceMedD5p1992-m.a", files, value = T)

d5med_cape_1992 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

View(d5med_cape_all$d5med_m1)
# need to merge these
med_cape_all <- map2(med_cape_1991, med_cape_2017, bind_rows)

d5med_cape_all <- map2(d5med_cape_1992, d5med_cape_2017, bind_rows)



titlenames <- c("m1", "m2", "m3", "m4", "m5", "m6")
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")

for(i in 1:length(titlenames)){
     mm <- capelinAreaPlot(ls1=cape, ls2=med_cape_all, ls3=d5med_cape_all, 
                              i, titlenames=titlenames)
     ggsave(mm, filename = paste0("figs/capelinArea/", titlenames[i], ".pdf"))
}

# seperate cape into time components
cape_1991 <- map(cape, ~filter(.x, year <= 1991))
cape_2017 <- map(cape, ~filter(.x, year > 1991))

titlenames <- c("m1-maxarea", "m2-maxarea", "m3-maxarea", "m4-maxarea", "m5-maxarea", "m6-maxarea")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(cape, cape_1991, cape_2017, "max_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"))
}

titlenames <- c("m1-medarea", "m2-medarea", "m3-medarea", "m4-medarea", "m5-medarea", "m6-medarea")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(med_cape_all, med_cape_1991, med_cape_2017, "med_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"))
}

titlenames <- c("m1-d5medarea", "m2-d5medarea", "m3-d5medarea", "m4-d5medarea", "m5-d5medarea", "m6-d5medarea")
for(i in 1:length(titlenames)){
     mm <- lnCapelinArea(d5med_cape_all, d5med_cape_1992, d5med_cape_2017, "d5med_area", "logcapelin", 1, titlenames)
     ggsave(mm, filename = paste0("figs/lncapelinArea/", titlenames[i], ".pdf"))
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
   CapelinDomeFit_area <- optim(par = c(2,400000,0.6),
                           dataf = cape$capelin_m1[which(cape$capelin_m1$logcapelin!='NA'),
                                           c('year','maxarea','logcapelin')],
                           fn = SSQCapelinDome, method=c("BFGS"))
   
   # before 2010
   CapelinDomeFitOld_area <- optim(par=c(2,400000,0.6),
                              dataf = cape$capelin_m1[which(cape$capelin_m1$year < 2011 & 
                                                                 cape$capelin_m1$logcapelin!='NA'),
                                              c('year','maxarea','logcapelin')],
                              fn = SSQCapelinDome, 
                              method=c("L-BFGS-B"), control=list(maxit=50000))
   
   ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
   cape$capelin_m1$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit_area$par), dataf = cape$capelin_m1[,c('year','maxarea')])
   
   cape$capelin_m1$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld_area$par), dataf = cape$capelin_m1[,c('year','maxarea')])
   
   
   
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
     geom_point(data = subset(cape$capelin_m1, year < 1991), aes(x = maxarea, y = logcapelin), shape=2, size=3) +
     geom_point(data = subset(cape$capelin_m1, year > 1991), aes(x = maxarea, y = logcapelin), shape=15, size=3) + 
     geom_errorbar(data = subset(cape$capelin_m1, year < 1991), aes(x = maxarea, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
     geom_errorbar(data = subset(cape$capelin_m1, year > 1991), aes(x = maxarea, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
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