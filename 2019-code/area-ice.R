#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-08, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Keith Lewis (2019-05-29) #

# The purpose of this file is to:

#1) to plot ice summaries (area and minlat) and ensure that the data are being properly summarized by ice-chart-processing-data
#2) examine the relationships between ice area, tice, and minlat 
#3) calculate tice for use in ice-capelin-jags_sequential.R
#4) produce *.csv files of the above, 


rm(list=ls())

#libraries----
library(ggplot2)
library(devtools)
library(ggmap)
library(sf)
library(tidyverse)
library(viridis)
library(rvest)
library(sp)
library(broom)
library(lubridate)
library(RArcInfo)
library(maptools)
library(rgdal)
library(rgeos)
library(data.table)
library(cowplot)
#devtools::install_github("tidyverse/ggplot2")
#library(raster)

setwd("D:\\Keith\\capelin\\2019-code")

# load code and data ----------
source("ice-chart-processing-function-v3.R")
source("area-ice-function.R")

#load("ice-trends-2017-m1-d3.Rdata")
#load data from ice-chart-processing-date-v3.R
load("output-processing/ice-trends-2019-m1-alla.Rdata")
# this file removes ice north of 55deg lat, Lake Melville, the Bays around NF, the Gulf of St. Lawrence, and the Scotian Shelf
load("output-processing/filters.Rdata")


# modify data---------
#get relevant data from trends.m1
#m1 <- trends.m1[c("date", "area", "volume")]

# add year to trends.m1
trends.m1$year <- year(trends.m1$date)

#Boxplots-----
#windows()
p1 <- ggplot(data = trends.m1, aes(x = year, y = area, group = year)) + 
  geom_boxplot()
p2 <- ggplot(data = trends.m1, aes(x = year, y = minlats, group = year)) + 
  geom_boxplot()

cowplot::plot_grid(p1, p2, labels = c("m1", ""), ncol=2)
ggsave("figs/1-area-minlatsBoxplots.pdf", width=10, height=8, units="in")


# yearly summaries of minlat and tice-------
# generate summaries by merging max area with minlat and tice and calculating the doy which gives tice
iceSum.m1 <- iceSummary(trends.m1)

## Fix 2013----
# note that 2013 was a very wonky ice year with 3-4 retreats and a large ice fragment on Apr 22 at < 49 deg latitude.  Because minlat was not a good proxy for tice in this year, I looked at the maps (D:\Keith\ice\Ice_Charts\Original_Ice_Data\GIF) and determined tice visually as Ale did in teh original work (Buren et al. 2014).  Note that it was a close call between February and April.  Basically, we decided that some low latitude ice hugging the coast in April but with little ice to the north in not the same as a lot of ice above the minlat in Feb and that the latter was more in accord with the spirit of Ale's work.

# from ice-chart-processing-date-v3.5 line207-220 
#confirm with AA that below is code is correct given changes taht have been made (the data and area are correct)
#minlats2013 <- trends.2013[8,]
#minlats2013
#    date     area       volume  minlats  minlongs
#8 2013-02-25 73609.46 23.49832 49.12323 -53.05903

# makes iceSum.m1 a list - this list and following loop were needed when there were 6 ice subsets.  Not necessary to do this way but it works so why change it?
iceSum_ls <- list(iceSum.m1=iceSum.m1)
# rm(ice_ls)

for(i in seq_along(iceSum_ls)){
     iceSum_ls[[i]]$date.y[45] <- "2013-02-25"
     iceSum_ls[[i]]$minlats[45] <- 49.12323
     iceSum_ls[[i]]$minlongs[45] <- -53.05903
     iceSum_ls[[i]]$jday.y[45] <- 105 #day of min lat
     iceSum_ls[[i]]$tice[45] <- yday(trends.m1$date[1313]) #day we chose as tice - see line 69 above
}

# puts the list "iceSum_ls into the global environment (we think)
list2env(iceSum_ls, .GlobalEnv)

##Plot year, tice, and area----
# plot year v tice and year v area for m1-m6

p1 <- ggplot(data = iceSum.m1, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m1")

p2 <- ggplot(data = iceSum.m1, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

cowplot::plot_grid(p1, p2, labels = c("m1", ""), ncol=2)
ggsave("figs/2-areaMaxScatterm1m3.pdf", width=10, height=8, units="in")


# ranges and ratios----------
#the following could be disgarded or brought together in a more useful manner, i.e. with the years for each object below.
area.range.low <- c(m1 = range(iceSum.m1$area)[1])
area.range.high <- c(m1 = range(iceSum.m1$area)[2])

# here is a stab at a more useful format - continue with tice and then join the objects based on year.
x <- subset(iceSum.m1, area == area.range.low)
y <- select(x, "year", "area")
x <- subset(iceSum.m1, area == area.range.high)
z <- select(x, "year", "area")
rbind(y,z)

tice.range.low <- c(m1 = range(iceSum.m1$tice)[1])
tice.range.high <- c(m1 = range(iceSum.m1$tice)[2])

ice.ranges <- as.data.frame(cbind(area.range.low, area.range.high, tice.range.low, tice.range.high))

area.ratio <- c(m1 = range(iceSum.m1$area)[2]/range(iceSum.m1$area)[1])
tice.ratio <- c(m1 = range(iceSum.m1$tice)[2]/range(iceSum.m1$tice)[1])

ice.ranges.ratios <- cbind(ice.ranges, area.ratio, tice.ratio)
write.csv(iceSum.m1, file = "output-processing/ice-ranges-ratios.csv", row.names=F, na="")


##Create maps of the minlats-----------
# plot water and minlongs/lats as points
# confirm that filters are working, i.e., this helps to check for tice and that there are no wonky minlats; check ice maps

# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

#m1 - plot for presentaiton only - save manually
par(mfrow=c(1,2))
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
# plot main map with filters; these are areas excluded from search for tice
plot(filters, border = "red", add = TRUE)  
points(iceSum.m1$minlongs, iceSum.m1$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
plot(filters, border = "red", add = TRUE)  # plot main map with filter
points(trends.m1$minlongs, trends.m1$minlats)


# Plot relationships between area, tice, and minlats FOR MAX_AREA-----
#m1
# plot tice against ice area
iceCorr.m1 <- iceScatterSummary(iceSum.m1, x = "area", x1 = "minlats", y = "tice", y1 = "minlats")
cowplot::plot_grid(iceCorr.m1$p1, iceCorr.m1$p2, iceCorr.m1$p3, labels = c("m1"), ncol=2)
ggsave("figs/3-relMax-m1.pdf", width=10, height=8, units="in")
#correlations for the above graph
iceCorr.m1.rsq <- lmoutputSummary(iceCorr.m1)
write.csv(iceCorr.m1.rsq, file = "output-processing/iceCorr-m1-rsq.csv", row.names=F, na="")


####################################################################################### OUTPUT of derived values for each year----
#3) calculate tice for use in ice-capelin-jags_sequential.R
#4) produce *.csv files of the above, 

iceSum.m1 <- iceSum.m1[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m1, file = "output-processing/capelin-m1.csv", row.names=F, na="")
