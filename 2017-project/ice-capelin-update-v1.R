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
source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
library(plotrix)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
capelin <- read.csv('capelin-ice-2014.csv',header=T)

capelin_m1 <- read_csv("output-processing/capelin-m1.csv")
capelin_m2 <- read_csv("output-processing/capelin-m2.csv")
capelin_m3 <- read_csv("output-processing/capelin-m3.csv")
capelin_m4 <- read_csv("output-processing/capelin-m4.csv")
capelin_m5 <- read_csv("output-processing/capelin-m5.csv")
capelin_m6 <- read_csv("output-processing/capelin-m6.csv")

  
## Ale original-------------
capelin <- capelin[c("year", "maxarea", "minlat", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]

capelin_join <- capelin[c("year", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
## m1-----------
capelin_m1 <- left_join(capelin_m1, capelin_join, by = "year")

## m2-----------
capelin_m2 <- capelin_m2[c("year", "area", "minlats", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
capelin_m2 <- rename(capelin_m2, maxarea = area)
capelin_m2 <- rename(capelin_m2, minlat = minlats)

## m3-----------
capelin_m3 <- capelin_m3[c("year", "area", "minlats", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
capelin_m3 <- rename(capelin_m3, maxarea = area)
capelin_m3 <- rename(capelin_m3, minlat = minlats)

## m4-----------
capelin_m4 <- capelin_m4[c("year", "area", "minlats", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
capelin_m4 <- rename(capelin_m4, maxarea = area)
capelin_m4 <- rename(capelin_m4, minlat = minlats)

capelin$myear <- capelin$year
capelin$myear[45:46] <- c(1960,1961)


## Ale's original-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present
CapelinDomeFit <- optim(par = c(1,200,0.6),
                        dataf = capelin[which(capelin$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome, method=c("BFGS"))

# before 2010
CapelinDomeFitOld <- optim(par=c(0.2,180,0.6),
                        dataf = capelin[which(capelin$year < 2011 & 
                                                capelin$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome, 
                        method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par), dataf = capelin[,c('year','tice')])

capelin$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par), dataf = capelin[,c('year','tice')])



# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000), tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????


# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

# set values
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

yearInt <- seq(1982, 2014, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)
labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

windows()
optimGraphs(capelin, regime1, regime2, yearInt, lnbiomassInt)


## m1--------------
CapelinDomeFit1 <- optim(par = c(1,200,0.6),
                        dataf = capelin_m1[which(capelin_m1$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome1, 
                        method=c("BFGS"))

# before 2010
CapelinDomeFitOld1 <- optim(par=c(0.2,180,0.6),
                            dataf = capelin_m1[which(capelin_m1$year < 2011 & 
                                                       capelin_m1$logcapelin!='NA'),
                            c('year','tice','logcapelin')],
                            fn = SSQCapelinDome1, 
                            method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin_m1$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit1$par), dataf = capelin_m1[,c('year','tice')])

capelin_m1$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld1$par), dataf = capelin_m1[,c('year','tice')])


labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000),tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld1$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit1$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????


# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

# set values
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

yearInt <- seq(1982, 2018, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)
windows()

optimGraphs(capelin_m1, regime1, regime2, yearInt, lnbiomassInt)

## m2--------------
CapelinDomeFit2 <- optim(par = c(1,200,0.6),
                        dataf = capelin_m2[which(capelin_m2$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome1, method=c("BFGS"))

# before 2010
CapelinDomeFitOld2 <- optim(par=c(0.2,180,0.6),
                            dataf = capelin_m2[which(capelin_m2$year < 2011 &                                         capelin_m2$logcapelin!='NA'),
                            c('year','tice','logcapelin')],
                            fn = SSQCapelinDome1, method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin_m2$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit2$par), dataf         = capelin_m2[,c('year','tice')])

capelin_m2$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld2$par),          dataf = capelin_m2[,c('year','tice')])


labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000),tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld2$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit2$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

windows()

optimGraphs(capelin_m2, regime1, regime2, yearInt, lnbiomassInt)

## m3--------------
CapelinDomeFit3 <- optim(par = c(1,200,0.6),
                         dataf = capelin_m3[which(capelin_m3$logcapelin!='NA'),
                         c('year','tice','logcapelin')],
                         fn = SSQCapelinDome1, method=c("BFGS"))

# before 2010
CapelinDomeFitOld3 <- optim(par=c(0.2,180,0.6), 
                            dataf = capelin_m3[which(capelin_m3$year < 2011 & capelin_m3$logcapelin!='NA'),
                            c('year','tice','logcapelin')],
                            fn = SSQCapelinDome1, method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin_m3$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit3$par), dataf = capelin_m3[,c('year','tice')])

capelin_m3$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld3$par), dataf = capelin_m3[,c('year','tice')])


labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000),tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld3$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit3$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

# set values

windows()

optimGraphs(capelin_m3, regime1, regime2, yearInt, lnbiomassInt)

## m4--------------
CapelinDomeFit4 <- optim(par = c(1,200,0.6),
                         dataf = capelin_m4[which(capelin_m4$logcapelin!='NA'),
                                            c('year','tice','logcapelin')],
                         fn = SSQCapelinDome1, method=c("BFGS"))




# before 2010
CapelinDomeFitOld4 <- optim(par=c(0.2,180,0.6), 
                            dataf = capelin_m4[which(capelin_m4$year < 2011 & capelin_m4$logcapelin!='NA'),
                                               c('year','tice','logcapelin')],
                            fn = SSQCapelinDome1, method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin_m4$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit4$par), dataf = capelin_m4[,c('year','tice')])

capelin_m4$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld4$par), dataf = capelin_m4[,c('year','tice')])


labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000),tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld4$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit4$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????

regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]
# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

# set values

windows()

optimGraphs(capelin_m4, regime1, regime2, yearInt, lnbiomassInt)


# AIC for models
2*CapelinDomeFit$value+2*length(CapelinDomeFit$par)
2*CapelinDomeFit1$value+2*length(CapelinDomeFit1$par)
2*CapelinDomeFit2$value+2*length(CapelinDomeFit2$par)
2*CapelinDomeFit3$value+2*length(CapelinDomeFit3$par)
2*CapelinDomeFit4$value+2*length(CapelinDomeFit4$par)
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
   
   
## Ale's original changed from tice to area-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present
   CapelinDomeFit_area <- optim(par = c(1,500,0.6),
                           dataf = capelin[which(capelin$logcapelin!='NA'),
                                           c('year','maxarea','logcapelin')],
                           fn = SSQCapelinDome, method=c("BFGS"))
   
   # before 2010
   CapelinDomeFitOld_area <- optim(par=c(0.2,500,0.6),
                              dataf = capelin[which(capelin$year < 2011 & 
                                                      capelin$logcapelin!='NA'),
                                              c('year','maxarea','logcapelin')],
                              fn = SSQCapelinDome, 
                              method=c("BFGS"))
   
   ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
   capelin$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit_area$par), dataf = capelin[,c('year','maxarea')])
   
   capelin$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld_area$par), dataf = capelin[,c('year','maxarea')])
   
   
   
   # attach the optimization curves of capelin abundance to ice data
   xarea <- expand.grid(year = c(1990,2000), area = c(0:500,573.515,587.768))
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
   labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
   
   windows()
ggplot() +
     geom_line(data = regime1, aes(x = area, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
     geom_line(data = regime2, aes(x = area, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
     geom_line(data = regime1, aes(x = area, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
     geom_line(data = regime2, aes(x = area, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
     geom_point(data = subset(capelin, year < 1991), aes(x = maxarea, y = logcapelin), shape=2, size=3) +
     geom_point(data = subset(capelin, year > 1991), aes(x = maxarea, y = logcapelin), shape=15, size=3) + 
     geom_errorbar(data = subset(capelin, year < 1991), aes(x = maxarea, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
     geom_errorbar(data = subset(capelin, year > 1991), aes(x = maxarea, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
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
   summary(lm(logcapelin ~ maxarea, data=capelin_m1))
   summary(lm(logcapelin ~ maxarea, data=capelin_m2))
   summary(lm(logcapelin ~ maxarea, data=capelin_m3))
   summary(lm(logcapelin ~ maxarea, data=capelin_m4))
   
   
   
################################################
   data_path <- "output-processing"
   files <- list.files(data_path)
   pattern <- grep("capelin", files, value = T)
   
   files
   test <- files %>%
     map(~ read_csv(file.path(data_path, .))) %>%
     reduce(rbind)
   
   
   
   calcAreaVolLat <- function(y, ct = NULL, sa = NULL, sb = NULL) 
     #browser()
     #my.env <- new.env()
     year <- maxarea <- minlats <- tice <-  capelin <-  capelinlb <-  capelinub <-  logcapelin <-  logcapelinlb <- logcapelinub <- 
       areas <- rep(NA, length(pattern)) # create an empty object
     
     for (i in seq_along(pattern)) {
       ## load map data
       print(pattern[i])                   #start here when making single object for testing
     }
     