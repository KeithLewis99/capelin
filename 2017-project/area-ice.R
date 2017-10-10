################################################################
#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-08, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Keith Lewis (2017-09-22) #
################################################################

# The purpose of this file is to: 

#1) to map these summaries and ensure that the data are being properly summarized by ice-chart-processing-data
#2) examine the relationships between ice area, tice, and minlat using various subsets of the ice data
#3) explore more rigorous method of determining relationship between tice and minlat


rm(list=ls())

#libraries----
#library(lubridate)
library(ggplot2)
#library(dplyr)
library(devtools)
#devtools::install_github("tidyverse/ggplot2")
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
library(raster)
library(data.table)

# load code and data ----------
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
#load("ice-trends-2017-m1-d3.Rdata") 
load("output-processing/ice-trends-2017-m1-all.Rdata")
load("output-processing/ice-trends-2017-m2-all.Rdata")
load("output-processing/ice-trends-2017-m3-all.Rdata")
load("output-processing/ice-trends-2017-m4-all.Rdata")
load("output-processing/ice-trends-2017-m5-all.Rdata")
load("output-processing/ice-trends-2017-m6-all.Rdata")
#load("output-processing/ice-trends-2017-m1-subset.Rdata")
#load("output-processing/ice-trends-2017-m2-subset.Rdata")
load("output-processing/filters.Rdata")


##################################
# check data to ensure proper subset---------
# plots to determien if subsetting worked - alternate with minlats/minlongs
m1 <- trends.m1[c("date", "area", "volume")]
m2 <- trends.m2[c("date", "area", "volume")]
m3 <- trends.m3[c("date", "area", "volume")]
m4 <- trends.m4[c("date", "area", "volume")]
m5 <- trends.m4[c("date", "area", "volume")]
m6 <- trends.m4[c("date", "area", "volume")]

# no value should be above the slope of 1
subsetTestPlot(m1, m2, "date")
subsetTestPlot(m2, m3, "date")
subsetTestPlot(m3, m4, "date")
subsetTestPlot(m4, m5, "date")
subsetTestPlot(m5, m6, "date")

# summary graphs
trends.m1$year <- year(trends.m1$date)
trends.m2$year <- year(trends.m2$date)
trends.m3$year <- year(trends.m3$date)
trends.m4$year <- year(trends.m4$date)

#Boxplots
iceYearBox(df1 = trends.m1, df2 = trends.m2, df3 = trends.m3, df4 = trends.m4)


# yearly summaries of minlat and tice
# generate summaries by merging max area with minlat and tice-------
iceSum.m1a <- iceSummary(trends.m1)
iceSum.m2 <- iceSummary(trends.m2)
iceSum.m3 <- iceSummary(trends.m3)
iceSum.m4 <- iceSummary(trends.m4)

# this is the null model for tice and area
View(iceSum.m1a)

p1 <- ggplot(data = iceSum.m1, aes(x = year, y = tice)) + 
  geom_point() + geom_smooth(method=lm)

p2 <- ggplot(data = iceSum.m1, aes(x = year, y = area)) + 
  geom_point() + geom_smooth(method=lm)

p3 <- ggplot(data = iceSum.m2, aes(x = year, y = tice)) + 
  geom_point() + geom_smooth(method=lm)

p4 <- ggplot(data = iceSum.m2, aes(x = year, y = area)) + 
  geom_point() + geom_smooth(method=lm)

p5 <- ggplot(data = iceSum.m3, aes(x = year, y = tice)) + 
  geom_point() + geom_smooth(method=lm)

p6 <- ggplot(data = iceSum.m3, aes(x = year, y = area)) + 
  geom_point() + geom_smooth(method=lm)

p7 <- ggplot(data = iceSum.m4, aes(x = year, y = tice)) + 
  geom_point() + geom_smooth(method=lm)

p8 <- ggplot(data = iceSum.m4, aes(x = year, y = area)) + 
  geom_point() + geom_smooth(method=lm)

# filter so that 0 values are removed
multiplot(p1, p3, p5, p7, p2, p4, p6, p8, cols = 2)

iceSum.m3 %>%
  filter(year >= 1993) %>%
ggplot(aes(x = year, y = area)) + 
  geom_point() + geom_smooth(method=lm)


# ranges and ratios
area.range <- c(m1 = range(iceSum.m1$area), m2 = range(iceSum.m2$area), m3 = range(iceSum.m3$area), m4 = range(iceSum.m4$area))


tice.range <- c(m1 = range(iceSum.m1$tice), m2 = range(iceSum.m2$tice), m3 = range(iceSum.m3$tice), m4 = range(iceSum.m4$tice))

rbind(area.range, tice.range)

area.ratio <- c(m1 = range(iceSum.m1$area)[2]/range(iceSum.m1$area)[1],
  m2 = range(iceSum.m2$area)[2]/range(iceSum.m2$area)[1],
  m3 = range(iceSum.m3$area)[2]/range(iceSum.m3$area)[1],
  m4 = range(iceSum.m4$area)[2]/range(iceSum.m4$area)[1]
)

tice.ratio <- c(m1 = range(iceSum.m1$tice)[2]/range(iceSum.m1$tice)[1],
  m2 = range(iceSum.m2$tice)[2]/range(iceSum.m2$tice)[1],
  m3 = range(iceSum.m3$tice)[2]/range(iceSum.m3$tice)[1],
  m4 = range(iceSum.m4$tice)[2]/range(iceSum.m4$tice)[1]
)

# filter
cbind(area.ratio, tice.ratio)

#######################################################################
##Create maps of the minlats-----------
#m1-----------
# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
par(mfrow=c(4,2))

plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m1$minlongs, iceSum.m1$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m1$minlongs, trends.m1$minlats)

#m2-----------
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter 
  points(iceSum.m2$minlongs, iceSum.m2$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m2$minlongs, trends.m2$minlats)

#m3-----------
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter 
  points(iceSum.m3$minlongs, iceSum.m3$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m3$minlongs, trends.m3$minlats)

#m4-----------
  water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
  water <- gUnaryUnion(water)
  
  # plot water and minlongs/lats as points 
  # confirm that filters are working; check ice maps
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter 
  points(iceSum.m4$minlongs, iceSum.m4$minlats)
  
  # plot all points
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m4$minlongs, trends.m4$minlats)
  
par(mfrow=c(1,1))

###############################################################################
# plot relationships between area, tice, and minlats
#m1-----
# plot tice against ice area
iceCorr.m1 <- iceScatterSummary(iceSum.m1, "area", "minlats", "tice", "minlats")

iceCorr.m1$s1 
iceCorr.m1$s2
iceCorr.m1$s3

multiplot(iceCorr.m1$p1, iceCorr.m1$p3, iceCorr.m1$p2, cols=2)

#m2-------
# plot tice against ice area
iceCorr.m2 <- iceScatterSummary(iceSum.m2)

iceCorr.m2$s1 
iceCorr.m2$s2
iceCorr.m2$s3

multiplot(iceCorr.m2$p1, iceCorr.m2$p3, iceCorr.m2$p2, cols=2)

#m3-----
iceCorr.m3 <- iceScatterSummary(iceSum.m3)

iceCorr.m3$s1 
iceCorr.m3$s2
iceCorr.m3$s3

multiplot(iceCorr.m3$p1, iceCorr.m3$p3, iceCorr.m3$p2, cols=2)

#m4-----
iceCorr.m4 <- iceScatterSummary(iceSum.m4)

iceCorr.m4$s1 
iceCorr.m4$s2
iceCorr.m4$s3

multiplot(iceCorr.m4$p1, iceCorr.m4$p3, iceCorr.m4$p2, cols=2)

#####################################################################
# explore more rigorous method of determining relationship between tice and minlat-----

#check max area

str(trends.m1)
# calculate tice and year
trends.m1$tice <- yday(trends.m1$date)
trends.m1$year <- year(trends.m1$date)

trends.m2$tice <- yday(trends.m2$date)
trends.m2$year <- year(trends.m2$date)

trends.m3$tice <- yday(trends.m3$date)
trends.m3$year <- year(trends.m3$date)

trends.m4$tice <- yday(trends.m4$date)
trends.m4$year <- year(trends.m4$date)

# Medians of all values
#####1991-----
#m1
sub1991.m1 <- iceMedian(trends.m1, "year < 1992", "tice < 150", iceSum.m1)
str(sub1991.m1)
iceSummarylm(sub1991.m1)

# graph these for this year but only repeat if reasonable r-squared values appear
p1.m1  <- ggplot(data = iceSum.m1, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=iceSum.m1))

# plot tice against ice area
p2.m1 <- ggplot(data = iceSum.m1, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=iceSum.m1))

#plot ice area against minlats
p3.m1 <- ggplot(data = iceSum.m1, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=iceSum.m1))

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub1991.m1.plot <- iceDateScatter(sub1991.m1)
sub1991.m1.plot$p1
sub1991.m1.plot$p2
sub1991.m1.plot$p3

#m2--------
sub1991.m2 <- iceMedian(trends.m2, "year < 1992", "tice < 150", iceSum.m2)
iceSummarylm(sub1991.m2)
# scatter plot of date and minlats
sub1991.m2.plot <- iceDateScatter(sub1991.m2)
sub1991.m2.plot$p1
sub1991.m2.plot$p2
sub1991.m2.plot$p3

#m3--------
sub1991.m3 <- iceMedian(trends.m3, "year < 1992", "tice < 150", iceSum.m3)
iceSummarylm(sub1991.m3)

#m4--------
sub1991.m4 <- iceMedian(trends.m4, "year < 1992", "tice < 150", iceSum.m4)
iceSummarylm(sub1991.m4)

##############################################################
#####1992-2017------
##m1
sub2017.m1 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
iceSummarylm(sub2017.m1)
sub2017.m1.plot <- iceDateScatter(sub2017.m1)
sub2017.m1.plot$p1
sub2017.m1.plot$p2
sub2017.m1.plot$p3

##m2
sub2017.m2 <- iceMedian(trends.m2, "year > 1991", "tice < 150", iceSum.m2)
iceSummarylm(sub2017.m2)
# scatter plot of date and minlats
sub2017.m2.plot <- iceDateScatter(sub2017.m2, d = "d")
sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

##m3
sub2017.m3 <- iceMedian(trends.m3, "year > 1991", "tice < 150", iceSum.m3)
iceSummarylm(sub2017.m3)

##m4
sub2017.m4 <- iceMedian(trends.m4, "year > 1991", "tice < 150", iceSum.m4)
iceSummarylm(sub2017.m4)

###############################################
###############################################
###############################################
# Median of the top 5 values

count.year <- trends.m1 %>%
  group_by(year) %>%
  count()
View(count.year)

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
###1992-2017------
##m1
iceMedD5p1992.m1 <- iceMedianD5(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
iceSummarylm(iceMedD5p1992.m1, med = "d5", med1 = "d5")

# figures
plotMedD5p1992.m1 <- iceScatterSummary(iceMedD5p1992.m1$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
multiplot(plotMedD5p1992.m1$p1, plotMedD5p1992.m1$p3, plotMedD5p1992.m1$p2, cols=2)

# scatter plot of date and minlats
windows()
sub2017.m2.plot <- iceDateScatter(iceMedD5p1992.m1, d= "d5")

sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

###1992-2017------
##m2
iceMedD5p1992.m2 <- iceMedianD5(trends.m2, "year > 1991", "tice < 150", iceSum.m2)
iceSummarylm(iceMedD5p1992.m2, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub2017.m2.plot <- iceDateScatter(iceMedD5p1992.m2, d= "d5")
sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

###1992-2017------
##m3
iceMedD5p1992.m3 <- iceMedianD5(trends.m3, "year > 1991", "tice < 150", iceSum.m3)
iceSummarylm(iceMedD5p1992.m3, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5p1992.m3$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5p1992.m3$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5p1992.m3$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub2017.m3.plot <- iceDateScatter(iceMedD5p1992.m3, d= "d5")
sub2017.m3.plot$p1
sub2017.m3.plot$p2
sub2017.m3.plot$p3

###1992-2017------
##m4
iceMedD5.m4 <- iceMedianD5(trends.m4, "year > 1991", "tice < 150", iceSum.m4)
iceSummarylm(iceMedD5.m4, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5.m4$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5.m4$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5.m4$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub2017.m4.plot <- iceDateScatter(iceMedD5.m4, d= "d5")
sub2017.m4.plot$p1
sub2017.m4.plot$p2
sub2017.m4.plot$p3

###pre-1992------
##m1
iceMedD5p1992.m1 <- iceMedianD5(trends.m1, "year <= 1992", "tice < 150", iceSum.m1)
iceSummarylm(iceMedD5p1992.m1, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5.m1$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5.m1$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5.m1$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)



# scatter plot of date and minlats
sub1992.m2.plot <- iceDateScatter(iceMedD5.m1, d = "d5")
sub1992.m2.plot$p1
sub1992.m2.plot$p2
sub1992.m2.plot$p3

###1992-2017------
##m2
iceMedD5p1992.m2 <- iceMedianD5(trends.m2, "year <= 1991", "tice < 150", iceSum.m2)
iceSummarylm(iceMedD5p1992.m2, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5p1992.m2$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub1992.m2.plot <- iceDateScatter(iceMedD5p1992.m2, d= "d5")
sub1992.m2.plot$p1
sub1992.m2.plot$p2
sub1992.m2.plot$p3

###1992-2017------
##m3
iceMedD5.m3 <- iceMedianD5(trends.m3, "year <= 1991", "tice < 150", iceSum.m3)
iceSummarylm(iceMedD5.m3, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5.m3$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5.m3$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5.m3$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub1992.m3.plot <- iceDateScatter(iceMedD5.m3, d= "d5")
sub1992.m3.plot$p1
sub1992.m3.plot$p2
sub1992.m3.plot$p3

###1992-2017------
##m4
iceMedD5.m4 <- iceMedianD5(trends.m4, "year <= 1991", "tice < 150", iceSum.m4)
iceSummarylm(iceMedD5.m4, med = "d5", med1 = "d5")

# figures
p1.m1  <- ggplot(data = iceMedD5.m4$mall, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5.m4$mall, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5.m4$mall, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
sub1992.m4.plot <- iceDateScatter(iceMedD5.m4, d= "d5")
sub1992.m4.plot$p1
sub1992.m4.plot$p2
sub1992.m4.plot$p3


source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
#######################################################################################

iceSum.m1 <- iceSum.m1[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m1, file = "output-processing/capelin-m1.csv", row.names=F, na="")


iceSum.m2 <- iceSum.m2[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m2, file = "output-processing/capelin-m2.csv", row.names=F, na="")

iceSum.m3 <- iceSum.m3[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m3, file = "output-processing/capelin-m3.csv", row.names=F, na="")

iceSum.m4 <- iceSum.m4[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m4, file = "output-processing/capelin-m4.csv", row.names=F, na="")

#######Median values-----------

iceSum.m1 <- iceSum.m1[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m1, file = "output-processing/capelin-m1.csv", row.names=F, na="")


iceSum.m2 <- iceSum.m2[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m2, file = "output-processing/capelin-m2.csv", row.names=F, na="")

iceSum.m3 <- iceSum.m3[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m3, file = "output-processing/capelin-m3.csv", row.names=F, na="")

iceSum.m4 <- iceSum.m4[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m4, file = "output-processing/capelin-m4.csv", row.names=F, na="")

#######D5Median values-----------
iceMedD5b1992.m1 <- iceMedD5$med
#[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m1, file = "output-processing/capelin-m1.csv", row.names=F, na="")


iceSum.m2 <- iceSum.m2[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m2, file = "output-processing/capelin-m2.csv", row.names=F, na="")

iceSum.m3 <- iceSum.m3[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m3, file = "output-processing/capelin-m3.csv", row.names=F, na="")

iceSum.m4 <- iceSum.m4[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m4, file = "output-processing/capelin-m4.csv", row.names=F, na="")

#######################################################################################
#######################################################################################
###This is how to make a map with a SPDF - much harder in ggplot - not sure its worth the effort
# add to data a new column termed "id" composed of the rownames of data
# https://susanejohnston.wordpress.com/2012/07/03/creating-a-large-scale-map-using-ggplot2-a-step-by-step-guide/ for ggplot
# for converting SPDF to DF http://mazamascience.com/WorkingWithData/?p=1494
water@data$id <- rownames(water@data)
str(water, max.level = 3)
# create a data.frame from our spatial object
waterPoints <- tidy(water, region = "id")

# merge the "fortified" data with the data from our spatial object
waterDF <- merge(waterPoints, water@data, by = "id")
#water1DF <- full_join(water1Points, water1@data, by = "id")


ggWater <- ggplot(data = waterDF, aes(x=long, y=lat, group=group)) + 
  # map "tears" w/0 the group
  geom_polygon(fill="lightsteelblue") + 
  geom_path(col="black") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey"), panel.grid.major = element_line(colour = "grey90")) + theme_bw() 


guides(fill=FALSE)

ggWater

########################################################################################################