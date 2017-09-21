
rm(list=ls())
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
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

#load("ice-trends-2017-m1-d3.Rdata") 
load("output-processing/ice-trends-2017-m1-all.Rdata")
load("output-processing/ice-trends-2017-m2-all.Rdata")
load("output-processing/ice-trends-2017-m1-subset.Rdata")
load("output-processing/ice-trends-2017-m2-subset.Rdata")
load("output-processing/ice_trends_subset-test.Rdata")

load("output-processing/filters.Rdata")

##################################
# plots to determien if subsetting worked - alternate with minlats/minlongs
m1 <- trends.m1[c("date", "area", "volume")]
m2 <- trends.m2[c("date", "area", "volume")]

# no value should be above the slope of 1
mtest <-merge(m1, m2, by ="date")
head(mtest)
plot(mtest$area.x, mtest$area.y)
abline(a=0, b = 1, col="red", lwd=3)

#################
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
maxIce.minLat <- full_join(temp1, temp2, by = "year")

#convert date.y to doy
maxIce.minLat$tice <- yday(maxIce.minLat$date.y)
lookAt(maxIce.minLat)
View(maxIce.minLat)

#########################################################################Create maps of the minlats
# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
plot(filters, border = "red", add = TRUE) # plot main map with filter
points(maxIce.minLat$minlongs, maxIce.minLat$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
plot(filters, border = "red", add = TRUE) # plot main map with filter
points(trends.m1$minlongs, trends.m1$minlats)

# plot tice against ice area
p1  <- ggplot(data = maxIce.minLat, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=maxIce.minLat))

# plot tice against ice area
p2 <- ggplot(data = maxIce.minLat, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=maxIce.minLat))

#plot ice area against minlats
p3 <- ggplot(data = maxIce.minLat, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=maxIce.minLat))

multiplot(p1, p3, p2, cols=2)

m1.all
p1  <- ggplot(data = m1.all, aes(x = area, y = dtice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=maxIce.minLat))

#####################################################################
# explore more rigorous method of determining relationship between tice and minlat

#check max area

str(trends.m1)
# calculate tice and year
trends.m1$tice <- yday(trends.m1$date)
trends.m1$year <- year(trends.m1$date)

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")


#####1991------
sub1991 <- iceMedian(trends.m1, "year < 1992", "tice < 150", maxIce.minLat)

summary(lm(dminlats~darea, data=sub1991$mall))
summary(lm(dminlats~dtice, data=sub1991$mall))
summary(lm(darea~dtice, data=sub1991$mall))

# scatter plot of date and minlats
ggplot(data = trends.m1, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub1991$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub1991$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991$mall, aes(xintercept = darea), colour = "red")
##############################################################
#####1992-2017------
sub2017 <- iceMedian(trends.m1, "year > 1991", "tice < 150", maxIce.minLat)
# scatter plot of date and minlats
ggplot(data = trends.m1, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub2017$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub2017$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017$mall, aes(xintercept = darea), colour = "red")

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

maxIce.minLat.1993 <- filter(maxIce.minLat, year < 1994)
m1.all <- left_join(maxIce.minLat.1993, m1.med, by = "year")

yr <- 1994
ti <- 150
df <- trends.m1
rm(yr, ti, df, df1)
operator1 <- "<"
do.call(operator1, list(get(t1), yr))