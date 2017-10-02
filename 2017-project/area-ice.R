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
load("output-processing/ice-trends-2017-m4-subset.Rdata")
#load("output-processing/ice-trends-2017-m1-subset.Rdata")
#load("output-processing/ice-trends-2017-m2-subset.Rdata")
load("output-processing/filters.Rdata")


##################################
# check data to ensure proper subset---------
# plots to determien if subsetting worked - alternate with minlats/minlongs
m1 <- trends.m1[c("date", "area", "volume")]
m2 <- trends.m2[c("date", "area", "volume")]
m3 <- trends.m3[c("date", "area", "volume")]

# no value should be above the slope of 1
subsetTestPlot(m1, m2, "date")
subsetTestPlot(m2, m3, "date")

# summary graphs
trends.m1$year <- year(trends.m1$date)
trends.m2$year <- year(trends.m2$date)
trends.m3$year <- year(trends.m3$date)

p1 <- ggplot(data = trends.m1, aes(x = year, y = area, group = year)) + 
  geom_boxplot()
p2 <- ggplot(data = trends.m1, aes(x = year, y = minlats, group = year)) + 
  geom_boxplot()
p3 <- ggplot(data = trends.m2, aes(x = year, y = area, group = year)) + 
  geom_boxplot()
p4 <- ggplot(data = trends.m2, aes(x = year, y = minlats, group = year)) + 
  geom_boxplot()
p5 <- ggplot(data = trends.m3, aes(x = year, y = area, group = year)) + 
  geom_boxplot()
p6 <- ggplot(data = trends.m3, aes(x = year, y = minlats, group = year)) + 
  geom_boxplot()

windows()

multiplot(p1, p3, p5, p2, p4, p6, cols=2)

# generate summaries by merging max area with minlat and tice-------
iceSum.m1 <- iceSummary(trends.m1)
iceSum.m2 <- iceSummary(trends.m2)
iceSum.m3 <- iceSummary(trends.m3)
View(iceSum.m1)

p1 <- ggplot(data = iceSum.m1, aes(x = year, y = tice)) + 
  geom_point()

p2 <- ggplot(data = iceSum.m1, aes(x = year, y = area)) + 
  geom_point()

p3 <- ggplot(data = iceSum.m2, aes(x = year, y = tice)) + 
  geom_point()

p4 <- ggplot(data = iceSum.m2, aes(x = year, y = area)) + 
  geom_point()

p5 <- ggplot(data = iceSum.m3, aes(x = year, y = tice)) + 
  geom_point()

p6 <- ggplot(data = iceSum.m3, aes(x = year, y = area)) + 
  geom_point()

multiplot(p1, p3, p5, p2, p4, p6, cols = 2)


area.range <- c(m1 = range(iceSum.m1$area), m2 = range(iceSum.m2$area), m3 = range(iceSum.m3$area))


tice.range <- c(m1 = range(iceSum.m1$tice), m2 = range(iceSum.m2$tice), m3 = range(iceSum.m3$tice))

rbind(area.range, tice.range)

area.ratio <- c(m1 = range(iceSum.m1$area)[2]/range(iceSum.m1$area)[1],
  m2 = range(iceSum.m2$area)[2]/range(iceSum.m2$area)[1],
  m3 = range(iceSum.m3$area)[2]/range(iceSum.m3$area)[1]
)

tice.ratio <- c(m1 = range(iceSum.m1$tice)[2]/range(iceSum.m1$tice)[1],
  m2 = range(iceSum.m2$tice)[2]/range(iceSum.m2$tice)[1],
  m3 = range(iceSum.m3$tice)[2]/range(iceSum.m3$tice)[1]
)

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
par(mfrow=c(3,2))

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
  
par(mfrow=c(1,1))

###############################################################################
# plot relationships between area, tice, and minlats
#m1-----
# plot tice against ice area
p1.m1  <- ggplot(data = iceSum.m1, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=iceSum.m1))

# plot tice against ice area
p2.m1 <- ggplot(data = iceSum.m1, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=iceSum.m1))

#plot ice area against minlats
p3.m1 <- ggplot(data = iceSum.m1, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=iceSum.m1))

multiplot(p1.m1, p3.m1, p2.m1, cols=2)

#m2-------
# plot tice against ice area
p1.m2  <- ggplot(data = iceSum.m2, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=iceSum.m2))

# plot tice against ice area
p2.m2 <- ggplot(data = iceSum.m2, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=iceSum.m2))

#plot ice area against minlats
p3.m2 <- ggplot(data = iceSum.m2, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=iceSum.m2))

multiplot(p1.m2, p3.m2, p2.m2, cols=2)

#m3-----
p1.m3  <- ggplot(data = iceSum.m3, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=iceSum.m3))

# plot tice against ice area
p2.m3 <- ggplot(data = iceSum.m3, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=iceSum.m3))

#plot ice area against minlats
p3.m3 <- ggplot(data = iceSum.m3, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=iceSum.m3))

multiplot(p1.m3, p3.m3, p2.m3, cols=2)

#####################################################################
# explore more rigorous method of determining relationship between tice and minlat-----

#check max area

str(trends.m1)
# calculate tice and year
trends.m1$tice <- yday(trends.m1$date)
trends.m1$year <- year(trends.m1$date)

trends.m2$tice <- yday(trends.m2$date)
trends.m2$year <- year(trends.m2$date)

#####1991-----
#m1
sub1991.m1 <- iceMedian(trends.m1, "year < 1992", "tice < 150", iceSum.m1)
str(sub1991)
iceSummarylm(sub1991)

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
sub1991.m1.plot <- datePlots(sub1991.m1)
sub1991.m1.plot$p1
sub1991.m1.plot$p2
sub1991.m1.plot$p3

#m2--------
sub1991.m2 <- iceMedian(trends.m2, "year < 1992", "tice < 150", iceSum.m1)
iceSummarylm(sub1991.m2)
# scatter plot of date and minlats
sub1991.m2.plot <- datePlots(sub1991.m2)
sub1991.m2.plot$p1
sub1991.m2.plot$p2
sub1991.m2.plot$p3

##############################################################
#####1992-2017------
##m1
sub2017.m1 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
iceSummarylm(sub2017.m1)
sub2017.m1.plot <- datePlots(sub2017.m1)
sub2017.m1.plot$p1
sub2017.m1.plot$p2
sub2017.m1.plot$p3

##m2
sub2017.m2 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
iceSummarylm(sub2017.m2)
# scatter plot of date and minlats
sub2017.m2.plot <- datePlots(sub2017.m2)
sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

###############################################
###############################################
###############################################
# all data - I haven't made this function subsettable yet
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

count.year <- trends.m1 %>%
  group_by(year) %>%
  count()
View(count.year)


iceMedD5.m1 <- iceMedianD5(trends.m1)

summary(lm(d5area~d5tice, data=iceMedD5.m1))
summary(lm(d5minlats~d5tice, data=iceMedD5.m1))
summary(lm(d5minlats~d5area, data=iceMedD5.m1))
#iceSummarylm(iceMedD5.m1)

p1.m1  <- ggplot(data = iceMedD5.m1, aes(x = d5area, y = d5tice)) + geom_point() + geom_smooth(method=lm)
#, colour = year < 1992)
# plot tice against ice area
p2.m1 <- ggplot(data = iceMedD5.m1, aes(x = d5minlats, y = d5tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = iceMedD5.m1, aes(x = d5area, y = d5minlats)) + geom_point() + geom_smooth(method=lm)

multiplot(p1.m1, p3.m1, p2.m1, cols=2)


# scatter plot of date and minlats
ggplot(data = trends.m1, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = trends.m1, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = iceMedD5.m1, aes(xintercept = d5minlats), colour = "red")

ggplot(data = trends.m1, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = iceMedD5.m1, aes(xintercept = d5area), colour = "red")

# testing
subset(trends.m1, year == 2012)

test <- trends.m1 %>%
  group_by(year) %>%
  arrange(minlats) %>%
  slice(1:5) %>%
  summarize(d5tice=median(tice))
View(test)


trends.m1 %>%
  group_by(year) %>%
  arrange(desc(area)) %>%
  slice(1:5) %>%
  summarize(d5area=median(area))


########################################################################################################
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
# filter by year



trends.m1.1993 <- filter(trends.m1, year < 1994 & tice <150)
str(trends.m1.1993)

m1.mminlat <- trends.m1.1993 %>%
  group_by(year) %>%
  summarize(dmedian = median(minlats))

m1.mtice <- trends.m1.1993 %>%
  group_by(year) %>%
  summarize(dtice = median(tice))

m1.marea <- trends.m1.1993 %>%
  group_by(year) %>%
  summarize(darea = median(area))

m1.med <- left_join(m1.mminlat, m1.mtice, by = "year")

iceSum.m1.1993 <- filter(iceSum.m1, year < 1994)
m1.all <- left_join(iceSum.m1.1993, m1.med, by = "year")

yr <- 1994
ti <- 150
df <- trends.m1
rm(yr, ti, df, df1)
operator1 <- "<"
do.call(operator1, list(get(t1), yr))

########################################################################################################

# create a new data set with max area of ice and minlat - to be plotted and compared to raw maps
# change model name for ease
trends <- trends.m1
lookAt(trends)



##modify trends
# add year
trends$year <- year(trends$date)
trend.check <- filter(trends, year=="1992")
test <- trend.check[order(trend.check$area), ]
lookAt(trend.check)


# calculate the maximum area for ice by year
trends1 <- trends[c("date", "area", "volume", "year")]
temp1 <- trends1 %>%
  group_by(year) %>%
  slice(which.max(area)) 
temp1 


# calculate the minimum latitude for ice by year
trends2 <- trends[c("date", "minlats", "minlongs", "year")]
temp2 <- trends2 %>%
  group_by(year) %>%
  slice(which.min(minlats)) 

# merge summarized dataframes
iceSum.m1 <- full_join(temp1, temp2, by = "year")

#convert date.y to doy
iceSum.m1$tice <- yday(iceSum.m1$date.y)
lookAt(iceSum.m1)
View(iceSum.m1)

