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

# Divide area by 1000
for(i in 1:length(cape)){
     cape[[i]]$max_area1000 <- cape[[i]]$max_area/1000
}

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

# Divide area by 1000
for(i in 1:length(med_cape_all)){
     med_cape_all[[i]]$med_area1000 <- med_cape_all[[i]]$med_area/1000
}


d5med_cape_all <- map2(d5med_cape_1992, d5med_cape_2017, bind_rows)
# Divide area by 1000
for(i in 1:length(med_cape_all)){
     d5med_cape_all[[i]]$d5med_area1000 <- d5med_cape_all[[i]]$d5med_area/1000
}

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


######################################################################
## OPTIMIZATION PLOTS
## common features to all plots------
# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

yearInt <- seq(1982, 2014, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)
yearLim <- c(1982, 2014)
labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures

#################################
## Ale's original-------------
## Optimization - produces lists of the optimization curve for use in figures below
# up to present

# tice
capelin_ale$myear <- capelin_ale$year
capelin_ale$myear[45:46] <- c(1960,1961)

range(capelin_ale$tice)

AleMaxTice <- calcFit(capelin_ale, var = "tice",
                     par = c(1, 200, 0.6), 
                     form1 = "Alpha*tmp*(1-(tmp/Beta))",
                     form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                     x_range = c(0:190,173.515,187.768)
                     )
#plot
optimGraphs(AleMaxTice$df, AleMaxTice$regime1, AleMaxTice$regime2, yearInt, lnbiomassInt, "AleMaxTice", "tice")
#ggsave("figs/optimization/AleMaxtice.pdf", width=10, height=8, units="in")

# summs of squares and convergence
AleMaxTice$cdf$value
AleMaxTice$cdf$convergence

# Max area
range(capelin_ale$maxarea)
AleMaxArea <- calcFit(capelin_ale, var = "maxarea",
                      par = c(1, 500, 0.6), 
                      form1 = "Alpha*tmp*(1-(tmp/Beta))",
                      form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                      x_range = c(0:650)
)

#plot
optimGraphs(AleMaxArea$df, AleMaxArea$regime1, AleMaxArea$regime2, yearInt, lnbiomassInt, "AleMaxArea", "maxarea")
ggsave("figs/optimization/AleMaxArea.pdf", width=10, height=8, units="in")

AleMaxArea$cdf$value
AleMaxArea$cdf$convergence
str(AleMaxArea, max.level = 2)
#######################################################
## Optimization for tice using maxarea from m1-m6----
titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")
yearInt <- seq(1982, 2017, by=4)
yearLim <- c(1982, 2017)

source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")

MaxTice <- calcFit_all(cape, titlenames, par = c(1, 200, 0.6), var = "tice", 
                         form1 = "Alpha*tmp*(1-(tmp/Beta))",
                         form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                         x_range = c(0:190,173.515,187.768))
str(MaxTice, max.level = 3)
#create and save graphs
for(i in 1:length(MaxTice$optim_ls)){
     df1 <- as.data.frame(MaxTice$optim_ls[[i]]$df)
     df2 <- as.data.frame(MaxTice$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(MaxTice$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "tice")
    ggsave(mm, filename = paste0("figs/optimization/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

# compare Sums of Squares
# $value reports the minimized residual sums of squares or minimized Likelihood
# http://www.magesblog.com/2013/03/how-to-use-optim-in-r.html
# http://plantecology.syr.edu/fridley/bio793/likelihood.html - shows how to do AIC
for(i in 1:length(MaxTice$optim_ls)){
     print(MaxTice$optim_ls[[i]]$cdf$value)
}

# print convergence - zero is the right answer
for(i in 1:length(MaxTice$optim_ls)){
     print(MaxTice$optim_ls[[i]]$cdf$convergence)
}

# create graphs in plotly to examine acutal data points.
mm <- optimGraphs(MaxTice$optim_ls$`MaxTice-m1`$df, MaxTice$optim_ls$`MaxTice-m1`$regime1, MaxTice$optim_ls$`MaxTice-m1`$regime2, yearInt, lnbiomassInt,  "MaxTice-m1", "tice")
library(plotly)
ggplotly(mm)  
testPlot <- function(df, reg1, reg2, yearInt, lnbiomassInt, title, var){
     browser()
     p3 <- ggplot() +
          geom_line(data = reg1, aes_string(x = var) + aes(y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
          geom_line(data = reg2, aes_string(x = var) + aes(y = ExpectedLogBiomass), colour="green", linetype=1, size=1.25) +
          geom_line(data = reg1, aes_string(x = var) + aes(y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
          geom_line(data = reg2, aes_string(x = var) + aes(y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
          geom_point(data = subset(df, year < 1991), aes_string(x = var) + aes(y = logcapelin), shape=2, size=3) +
          geom_point(data = subset(df, year > 1991), aes_string(x = var) + aes(y = logcapelin), shape=15, size=3) + 
          geom_errorbar(data = subset(df, year < 1991),  aes_string(x = var) + aes(ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
          geom_errorbar(data = subset(df, year > 1991), aes_string(x = var) + aes(ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
          xlab(paste(var)) +
          ylab("ln (Capelin biomass (ktons))") + 
          #ylim(0,9) +
          theme_bw()
}

test <- testPlot(MaxTice$optim_ls$`MaxTice-m1`$df, MaxTice$optim_ls$`MaxTice-m1`$regime1, MaxTice$optim_ls$`MaxTice-m1`$regime2, yearInt, lnbiomassInt,  "MaxTice-m1", "tice")
ggplotly(test)
View(MaxTice$optim_ls$`MaxTice-m1`$df)

#############################################################################################
## Optimization for max_area using maxarea----
titlenames <- c("MaxArea-m1", "MaxArea-m2", "MaxArea-m3", "MaxArea-m4", "MaxArea-m5", "MaxArea-m6")

range(cape$capelin_m1$max_area)
range(cape$capelin_m1$max_area1000)
MaxArea <- calcFit_all(cape, titlenames, par = c(50000, 400000, 0.6), 
                       var = "max_area", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))", 
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = seq(from = 50000, to = 400000, by = 500))

MaxArea <- calcFit_all(cape, titlenames, par = c(1, 400, 0.6), 
                       var = "max_area1000", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))", 
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:450))
#create and save graphs
for(i in 1:length(MaxArea$optim_ls)){
     df1 <- as.data.frame(MaxArea$optim_ls[[i]]$df)
     df2 <- as.data.frame(MaxArea$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(MaxArea$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "max_area1000")
     ggsave(mm, filename = paste0("figs/optimization/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

# compare Sums of Squares
for(i in 1:length(MaxArea$optim_ls)){
     print(MaxArea$optim_ls[[i]]$cdf$value)
}

# print convergence - zero is the right answer
for(i in 1:length(MaxArea$optim_ls)){
     print(MaxArea$optim_ls[[i]]$cdf$convergence)
}


######################################################################
## Optimization for med_area using maxarea----
titlenames <- c("MedArea-m1", "MedArea-m2", "MedArea-m3", "MedArea-m4", "MedArea-m5", "MedArea-m6")

range(med_cape_all$med_m1$med_area1000)

MedArea <- calcFit_all(med_cape_all, titlenames, par = c(1, 300, 0.6),                        var = "med_area1000", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))", 
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:210))

#create and save graphs
for(i in 1:length(MedArea$optim_ls)){
     df1 <- as.data.frame(MedArea$optim_ls[[i]]$df)
     df2 <- as.data.frame(MedArea$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(MedArea$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "med_area1000")
     ggsave(mm, filename = paste0("figs/optimization/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

# compare Sums of Squares
for(i in 1:length(MedArea$optim_ls)){
     print(MedArea$optim_ls[[i]]$cdf$value)
}

# print convergence - zero is the right answer
for(i in 1:length(MedArea$optim_ls)){
     print(MedArea$optim_ls[[i]]$cdf$convergence)
}

######################################################################
## Optimization for med_area using maxarea----
titlenames <- c("D5MedArea-m1", "D5MedArea-m2", "D5MedArea-m3", "D5MedArea-m4", "D5MedArea-m5", "D5MedArea-m6")

range(d5med_cape_all$d5med_m1$d5med_area1000)

d5MedArea <- calcFit_all(d5med_cape_all, titlenames, par = c(1, 400, 0.6),                        var = "d5med_area1000", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))", 
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:400))
#create and save graphs
for(i in 1:length(d5MedArea$optim_ls)){
     df1 <- as.data.frame(d5MedArea$optim_ls[[i]]$df)
     df2 <- as.data.frame(d5MedArea$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(d5MedArea$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "d5med_area1000")
     ggsave(mm, filename = paste0("figs/optimization/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

# compare Sums of Squares
for(i in 1:length(d5MedArea$optim_ls)){
     print(d5MedArea$optim_ls[[i]]$cdf$value)
}

# print convergence - zero is the right answer
for(i in 1:length(d5MedArea$optim_ls)){
     print(d5MedArea$optim_ls[[i]]$cdf$convergence)
}

#####################################################################################################  NOT SURE WHAT THE BELOW IF FOR

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
   
################################################################################################