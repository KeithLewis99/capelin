#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-31, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to:
# 1) Examine the Pseudocalanus data from H. Murphy (emailed 2017-11-07: file "Copy of PSEUSP27_1999_2016.xlsx" - no metadata - see Murhpy et al. 2017) - created metadata for xlsx and csv file
# - days of year, range of data, years, etc


rm(list=ls())

## libraries------
#library(plotrix)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
#library(purrr)
#library(plotly)
library(lubridate)

## read in source code-----
#source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
#source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")

## load original data----
pscal <- read_csv("data/pseudocal_1999_2016.csv")
str(pscal)

# Murhpy et al used mean Pseudocalanus spp density from 1 June to 1 October - use lubridate to determine these values
yday("2017-06-01") # 152
yday("2017-10-01") # 274

#Psuedocalanus - count of trawls on doy by year - didn't really want this
pscal %>%
     ggplot() +
     geom_histogram(aes(doy)) +
     facet_wrap(~year)

# means and sd by year
m1 <- pscal %>%
     group_by(year) %>%
     summarise(mean_tot = mean(total), sd_tot = sd(total))

# gather pscal
pscal_long <- pscal %>%
     gather(key = stage, value = density, 4:10) 
str(pscal_long)
range(pscal_long$density)
mean(pscal_long$density)
sd(pscal_long$density)
#spread(new1, key= stage, value=density)

# means and sd by year and stage
m2 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     group_by(year, stage) %>%
     summarise(mean_tot = mean(density), sd_tot = sd(density)) %>%
     select(c(1:3)) %>%
     spread(key = stage, value = mean_tot)

# means and sd by year and stage; total only
m3 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage == "total") %>%
     group_by(year) %>%
     summarise(mean_tot = mean(density), sd_tot = sd(density))



#plot of above - densities are similar among stages    
p1 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     group_by(year, stage) %>%
     summarise(mean_tot = mean(density), sd_tot = sd(density)) %>%
     ggplot() +
     geom_point(aes(x = stage, y = mean_tot)) + 
     geom_errorbar(aes(x = stage, ymin = mean_tot-sd_tot, ymax = mean_tot+sd_tot)) +
     facet_wrap(~year)

#plot of density by doy by stage - c1-c6 are similar
p2 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     ggplot() +
     geom_point(aes(x = doy, y = density)) + 
     facet_wrap(~stage)

#plot of density of total by doy
p3 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage == "total") %>%
     ggplot() +
     geom_point(aes(x = doy, y = density, colour = stage))

#plot of density by stage
p4 <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage != "total") %>%
     ggplot() +
     geom_point(aes(x = doy, y = density, colour = stage))
