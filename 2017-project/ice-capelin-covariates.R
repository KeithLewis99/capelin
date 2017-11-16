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
cape_2001 <- map(cape, ~filter(.x, year > 2002 & year < 2017))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- cape_2001$capelin_m1[c("year", "tice", "logcapelin", "surface_tows_lag2")]

windows()
ggplot(data=cape_2001$capelin_m6) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=tice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)

########################################################################
# values for optim functions
yearInt <- seq(2000, 2017, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")

MaxTice <- calcFit_all(cape_2001, titlenames, par = c(1, 200, 0.6), var = "tice", 
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
     ggsave(mm, filename = paste0("figs/covariates/opt_", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

#################################
# IN DEVELOPMENT
#
# test of a generalization of the formula to multiple variables
source("D:/Keith/capelin/20117-project/ice-capelin-covariates-FUN.R")
test <- calcFit_all1(cape_2001, titlenames, par = c(5, 300, 5, 0.6), var1 = "tice", var2 = "surface_tows_lag2",
                       form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Delta*tmp2",
                       form2 = "Alpha*tmp1*(1-(tmp1/Beta))*Gamma + Delta*tmp2",
                       x_range = c(400:1500))

test <- calcFit_all(cape_2001, titlenames, par = c(1, 10000, 0.6), var = "surface_tows_lag2", 
                     form1 = "Alpha*tmp*(1-(tmp/Beta))",
                     form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                     x_range = c(400:4000))
str(MaxTice, max.level = 3)

Alpha <- 1
Beta <- 300
Delta <- 1
tmp1 <- cape_2001$capelin_m1$tice[2:15]
tmp2 <- cape_2001$capelin_m1$surface_tows_lag2[2:15]
Alpha*tmp1*(1-(tmp1/Beta)) + Delta*tmp2
View(cape_2001$capelin_m1)
#################################
## med area -----

# read the datasets and join them to capelin_join
pattern <- grep("sub2017-m+[^rsq].csv", files, value = T)
med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(med_cape_2017)){
     med_cape_2017[[i]]$med_area1000 <- med_cape_2017[[i]]$med_area/1000
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

#################################
## D5 med area----
pattern <- grep("iceMedD5p2017-m.a", files, value = T)

d5med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(d5med_cape_2017)){
     d5med_cape_2017[[i]]$d5med_area1000 <- d5med_cape_2017[[i]]$d5med_area/1000
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

titlenames <- c("D5MedArea-m1", "D5MedArea-m2", "D5MedArea-m3", "D5MedArea-m4", "D5MedArea-m5", "D5MedArea-m6")

for(i in 1:seq_along(d5med_cape_2001)){
mm <- ggplot(data=d5med_cape_2001[[i]]) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
          geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=d5tice), size=3) + 
          geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)
     ggsave(mm, filename = paste0("figs/covariates/scatter/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}
