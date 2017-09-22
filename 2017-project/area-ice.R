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


# load code and data ----------
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
#load("ice-trends-2017-m1-d3.Rdata") 
load("output-processing/ice-trends-2017-m1-all.Rdata")
load("output-processing/ice-trends-2017-m2-all.Rdata")
#load("output-processing/ice-trends-2017-m1-subset.Rdata")
#load("output-processing/ice-trends-2017-m2-subset.Rdata")
load("output-processing/filters.Rdata")


##################################
# check data to ensure proper subset---------
# plots to determien if subsetting worked - alternate with minlats/minlongs
m1 <- trends.m1[c("date", "area", "volume")]
m2 <- trends.m2[c("date", "area", "volume")]

# no value should be above the slope of 1
mtest <-merge(m1, m2, by ="date")
head(mtest)
plot(mtest$area.x, mtest$area.y)
abline(a=0, b = 1, col="red", lwd=3)

# generate summaries-------
iceSum.m1 <- iceSummary(trends.m1)
iceSum.m2 <- iceSummary(trends.m2)

#######################################################################
##Create maps of the minlats-----------
#m1-----------
# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
q1 <- plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
q1 + plot(filters, border = "red", add = TRUE) # plot main map with filter
q1 + points(iceSum.m1$minlongs, iceSum.m1$minlats)

# plot all points
q2 <- plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
q2 + plot(filters, border = "red", add = TRUE) # plot main map with filter
q2 + points(trends.m1$minlongs, trends.m1$minlats)

#m2-----------
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
# confirm that filters are working; check ice maps
q1 <- plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
q1 + plot(filters, border = "red", add = TRUE) # plot main map with filter
q1 + points(iceSum.m2$minlongs, iceSum.m2$minlats)

# plot all points
q2 <- plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA) 
q2 + plot(filters, border = "red", add = TRUE) # plot main map with filter
q2 + points(trends.m2$minlongs, trends.m2$minlats)

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
#####################################################################
# explore more rigorous method of determining relationship between tice and minlat-----

#check max area

str(trends.m1)
# calculate tice and year
trends.m1$tice <- yday(trends.m1$date)
trends.m1$year <- year(trends.m1$date)

trends.m2$tice <- yday(trends.m2$date)
trends.m2$year <- year(trends.m2$date)

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")


#####1991-----
#m1
sub1991 <- iceMedian(trends.m1, "year < 1992", "tice < 150", iceSum.m1)
str(sub1991)
summary(lm(dminlats~darea, data=sub1991$mall))
summary(lm(dminlats~dtice, data=sub1991$mall))
summary(lm(darea~dtice, data=sub1991$mall))

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
ggplot(data = sub1991$data, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub1991$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub1991$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991$mall, aes(xintercept = darea), colour = "red")

#m2--------
sub1991.m2 <- iceMedian(trends.m2, "year < 1992", "tice < 150", iceSum.m1)

summary(lm(dminlats~darea, data=sub1991.m2$mall))
summary(lm(dminlats~dtice, data=sub1991.m2$mall))
summary(lm(darea~dtice, data=sub1991.m2$mall))

# scatter plot of date and minlats
ggplot(data = sub1991.m2$data, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub1991.m2$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991.m2$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub1991.m2$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub1991.m2$mall, aes(xintercept = darea), colour = "red")

##############################################################
#####1992-2017------
##m1
sub2017.m1 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)

summary(lm(dminlats~darea, data=sub2017.m1$mall))
summary(lm(dminlats~dtice, data=sub2017.m1$mall))
summary(lm(darea~dtice, data=sub2017.m1$mall))

# scatter plot of date and minlats------
ggplot(data = sub2017.m1$data, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub2017.m1$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017.m1$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub2017$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017.m1$mall, aes(xintercept = darea), colour = "red")

##m2
sub2017.m2 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)

summary(lm(dminlats~darea, data=sub2017.m2$mall))
summary(lm(dminlats~dtice, data=sub2017.m2$mall))
summary(lm(darea~dtice, data=sub2017.m2$mall))

# scatter plot of date and minlats
ggplot(data = sub2017.m2$data, aes(x = date, y = minlats)) + geom_point() + geom_smooth(method=lm)

ggplot(data = sub2017.m2$data, aes(minlats)) + 
  geom_histogram() + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017.m2$mall, aes(xintercept = dminlats), colour = "red")

ggplot(data = sub2017.m2$data, aes(area)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ year) +
  geom_vline(data = sub2017.m2$mall, aes(xintercept = darea), colour = "red")

###############################################
# all data - I haven't made this function subsettable yet

count.year <- trends.m1 %>%
  group_by(year) %>%
  count()
View(count.year)

iceMedD5.m1 <- iceMedianD5(trends.m1)

summary(lm(d5area~d5tice, data=iceMedD5.m1))
summary(lm(d5minlats~d5tice, data=iceMedD5.m1))
summary(lm(d5minlats~d5area, data=iceMedD5.m1))


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



## testing
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

