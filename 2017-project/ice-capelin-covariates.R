################################################################
#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-08, R version 3.3.3 (2017-03-06)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-11.08 #
################################################################

# The purpose of this file is to:
# 1) Test additional covariates to improve the Tice model fit

rm(list=ls())
#setwd("D:/Keith/capelin/2017-project")

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(plotly)

## read in source code-----
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
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

# read in larval data
larvae <- read.csv('capelin_larval_indices.csv',header=T)
larvae$surface_tows_lag2 <- lag(larvae$surface_tows, 2)
capelin_join <- left_join(capelin_join, larvae, by = "year")

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
cape_2001 <- map(cape, ~filter(.x, year > 2001))

ggplot() + geom_point(data=cape_2001$capelin_m1, aes(x = surface_tows_lag2, y = logcapelin))
########################################################################

yearInt <- seq(2000, 2017, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")

MaxTice <- calcFit_all(cape_2001, titlenames, par = c(1, 200, 0.6), var = "tice*surface_tows_lag2", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))",
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:190,173.515,187.768))
str(MaxTice, max.level = 3)

source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
#create and save graphs
for(i in 1:length(MaxTice$optim_ls)){
     df1 <- as.data.frame(MaxTice$optim_ls[[i]]$df)
     df2 <- as.data.frame(MaxTice$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(MaxTice$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "tice")
     ggsave(mm, filename = paste0("figs/covariates/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}
