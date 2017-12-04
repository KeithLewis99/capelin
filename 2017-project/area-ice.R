################################################################
#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-08, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Keith Lewis (2017-09-22) #
################################################################

# The purpose of this file is to:

#1) to map these summaries and ensure that the data are being properly summarized by ice-chart-processing-data
#2) examine the relationships between ice area, tice, and minlat using various subsets of the ice data
#3) explore more rigorous method of determining relationship between tice and minlat such as median and median of the top 5 values
#4) produce *.csv files of the above


rm(list=ls())

#libraries----
library(ggplot2)
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
#library(raster)
library(data.table)
#library(cowplot)

# load code and data ----------
#source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
source("D:/Keith/capelin/2017-project/area-ice-function.R")

#load("ice-trends-2017-m1-d3.Rdata")
load("output-processing/ice-trends-2017-m1-alla.Rdata")
load("output-processing/ice-trends-2017-m2-all.Rdata")
load("output-processing/ice-trends-2017-m3-all.Rdata")
load("output-processing/ice-trends-2017-m4-all.Rdata")
load("output-processing/ice-trends-2017-m5-all.Rdata")
load("output-processing/ice-trends-2017-m6-all.Rdata")
#load("output-processing/ice-trends-2017-m1-subset.Rdata")
#load("output-processing/ice-trends-2017-m2-subset.Rdata")
load("output-processing/filters.Rdata")

##################################
# modify data---------
m1 <- trends.m1[c("date", "area", "volume")]
m2 <- trends.m2[c("date", "area", "volume")]
m3 <- trends.m3[c("date", "area", "volume")]
m4 <- trends.m4[c("date", "area", "volume")]
m5 <- trends.m5[c("date", "area", "volume")]
m6 <- trends.m6[c("date", "area", "volume")]

# summary graphs
trends.m1$year <- year(trends.m1$date)
trends.m2$year <- year(trends.m2$date)
trends.m3$year <- year(trends.m3$date)
trends.m4$year <- year(trends.m4$date)
trends.m5$year <- year(trends.m5$date)
trends.m6$year <- year(trends.m6$date)


# check data to ensure proper subset---------
# no value should be above the slope of 1
# m4 and m5 area not subsets of each other so did not plot
ss1 <- subsetTestPlot(m1, m2, "date")
ss2 <- subsetTestPlot(m2, m3, "date")
ss3 <- subsetTestPlot(m3, m4, "date")
ss4 <- subsetTestPlot(m3, m5, "date")#
#subsetTestPlot(m4, m5, "date")
ss5 <- subsetTestPlot(m5, m6, "date")
ss6 <- subsetTestPlot(m4, m6, "date")

# value is the x axis first.  So m1-m2 is x=m1 and y = m2
testSubset <- cowplot::plot_grid(ss1, ss2, ss3, ss4, ss5, ss6, labels = c("m1-2", "m2-3", "m3-4", "m3-5", "m5-6", "m4-6"), ncol=2)
ggsave("figs/1-testSubset.pdf", width=10, height=8, units="in")

#Boxplots-----
#windows()
iceYearBox(df1 = trends.m1, df2 = trends.m2, df3 = trends.m3)
ggsave("figs/2-allyearsBox-m1m3.pdf", width=10, height=8, units="in")
iceYearBoxm4m6 <- iceYearBox(df1 = trends.m4, df2 = trends.m5, df3 = trends.m6, title = "m4-m6")
ggsave("figs/3-allyearsBox-m4m6.pdf", width=10, height=8, units="in")

# yearly summaries of minlat and tice-------
# generate summaries by merging max area with minlat and tice
iceSum.m1 <- iceSummary(trends.m1)
iceSum.m2 <- iceSummary(trends.m2)
iceSum.m3 <- iceSummary(trends.m3)
iceSum.m4 <- iceSummary(trends.m4)
iceSum.m5 <- iceSummary(trends.m5)
iceSum.m6 <- iceSummary(trends.m6)
# this is the null model for tice and area

# Fix 2013
# from ice-capelin-data processing
#minlats2013 <- trends.2013[8,]
#minlats2013
#    date     area       volume  minlats  minlongs
#8 2013-02-25 73609.46 23.49832 49.12323 -53.05903

iceSum_ls <- list(iceSum.m1=iceSum.m1, iceSum.m2=iceSum.m2, iceSum.m3=iceSum.m3, iceSum.m4=iceSum.m4, iceSum.m5=iceSum.m5, iceSum.m6=iceSum.m6)
rm(ice_ls)

for(i in seq_along(iceSum_ls)){
     iceSum_ls[[i]]$date.y[45] <- "2013-02-25"
     iceSum_ls[[i]]$minlats[45] <- 49.12323
     iceSum_ls[[i]]$minlongs[45] <- -53.05903
     iceSum_ls[[i]]$jday.y[45] <- 105
     iceSum_ls[[i]]$tice[45] <- yday(trends.m1$date[1313])
}

list2env(iceSum_ls, .GlobalEnv)

# plot year v tice and year v area for m1-m6
p1 <- ggplot(data = iceSum.m1, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m1")

p2 <- ggplot(data = iceSum.m1, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

p3 <- ggplot(data = iceSum.m2, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m2")

p4 <- ggplot(data = iceSum.m2, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

p5 <- ggplot(data = iceSum.m3, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m3")

p6 <- ggplot(data = iceSum.m3, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

p7 <- ggplot(data = iceSum.m4, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m4")

p8 <- ggplot(data = iceSum.m4, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

p9 <- ggplot(data = iceSum.m5, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m5")

p10 <- ggplot(data = iceSum.m5, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

p11 <- ggplot(data = iceSum.m6, aes(x = year, y = tice)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m6")

p12 <- ggplot(data = iceSum.m6, aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm)

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, labels = c("m1", "", "m2", "", "m3", ""), ncol=2)
ggsave("figs/4-areaMaxScatterm1m3.pdf", width=10, height=8, units="in")
cowplot::plot_grid(p7, p8, p9, p10, p11, p12, ncol=2)
ggsave("figs/5-areaMaxScatterm4m5.pdf", width=10, height=8, units="in")

# filter so that 0 values are removed
si1 <- iceSum.m4 %>%
  filter(year >= 1983) %>%
ggplot(aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m4 > 1983")

si2 <- iceSum.m6 %>%
  filter(year >= 1983) %>%
  ggplot(aes(x = year, y = area)) +
  geom_point() + geom_smooth(method=lm) + ggtitle("m6 > 1983")

areaMaxScatterm4m6 <- cowplot::plot_grid(si1, si2, ncol=2)
ggsave("figs/6-areaMaxScatterm4m6.pdf", width=10, height=8, units="in")

# ranges and ratios----------
area.range.low <- c(m1 = range(iceSum.m1$area)[1],
                m2 = range(iceSum.m2$area)[1],
                m3 = range(iceSum.m3$area)[1],
                m4 = range(iceSum.m4$area)[1],
                m5 = range(iceSum.m5$area)[1],
                m6 = range(iceSum.m6$area)[1])

area.range.high <- c(m1 = range(iceSum.m1$area)[2],
                    m2 = range(iceSum.m2$area)[2],
                    m3 = range(iceSum.m3$area)[2],
                    m4 = range(iceSum.m4$area)[2],
                    m5 = range(iceSum.m5$area)[2],
                    m6 = range(iceSum.m6$area)[2])

tice.range.low <- c(m1 = range(iceSum.m1$tice)[1],
                m2 = range(iceSum.m2$tice)[1],
                m3 = range(iceSum.m3$tice)[1],
                m4 = range(iceSum.m4$tice)[1],
                m5 = range(iceSum.m5$tice)[1],
                m6 = range(iceSum.m6$tice)[1])
tice.range.high <- c(m1 = range(iceSum.m1$tice)[2],
                    m2 = range(iceSum.m2$tice)[2],
                    m3 = range(iceSum.m3$tice)[2],
                    m4 = range(iceSum.m4$tice)[2],
                    m5 = range(iceSum.m5$tice)[2],
                    m6 = range(iceSum.m6$tice)[2])

ice.ranges <- as.data.frame(cbind(area.range.low, area.range.high, tice.range.low, tice.range.high))

area.ratio <- c(m1 = range(iceSum.m1$area)[2]/range(iceSum.m1$area)[1],
  m2 = range(iceSum.m2$area)[2]/range(iceSum.m2$area)[1],
  m3 = range(iceSum.m3$area)[2]/range(iceSum.m3$area)[1],
  m4 = range(iceSum.m4$area)[2]/range(iceSum.m4$area)[1],
  m5 = range(iceSum.m5$area)[2]/range(iceSum.m5$area)[1],
  m6 = range(iceSum.m6$area)[2]/range(iceSum.m6$area)[1]
)

tice.ratio <- c(m1 = range(iceSum.m1$tice)[2]/range(iceSum.m1$tice)[1],
  m2 = range(iceSum.m2$tice)[2]/range(iceSum.m2$tice)[1],
  m3 = range(iceSum.m3$tice)[2]/range(iceSum.m3$tice)[1],
  m4 = range(iceSum.m4$tice)[2]/range(iceSum.m4$tice)[1],
  m5 = range(iceSum.m5$tice)[2]/range(iceSum.m5$tice)[1],
  m6 = range(iceSum.m6$tice)[2]/range(iceSum.m6$tice)[1]
)

ice.ranges.ratios <- cbind(ice.ranges, area.ratio, tice.ratio)
write.csv(iceSum.m1, file = "output-processing/ice-ranges-ratios.csv", row.names=F, na="")

#######################################################################
##Create maps of the minlats-----------
# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)


#m1 - plot for presentaiton only - save manually
par(mfrow=c(1,2))
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
plot(filters, border = "red", add = TRUE)  # plot main map with filter
points(iceSum.m1$minlongs, iceSum.m1$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
plot(filters, border = "red", add = TRUE)  # plot main map with filter
points(trends.m1$minlongs, trends.m1$minlats)

# plot water and minlongs/lats as points
# confirm that filters are working; check ice maps
par(mfrow=c(3,2))

#m1
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m1$minlongs, iceSum.m1$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m1$minlongs, trends.m1$minlats)

#m2
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m2$minlongs, iceSum.m2$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m2$minlongs, trends.m2$minlats)

#m3
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m3$minlongs, iceSum.m3$minlats)

# plot all points
plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m3$minlongs, trends.m3$minlats)

par(mfrow=c(3,2))
#m4
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m4$minlongs, iceSum.m4$minlats)

  # plot all points
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m4$minlongs, trends.m4$minlats)

#m5
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m5$minlongs, iceSum.m5$minlats)

  # plot all points
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m5$minlongs, trends.m5$minlats)

#m6
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(iceSum.m6$minlongs, iceSum.m6$minlats)

  # plot all points
  plot(water, xlim = c(-70, -48), ylim = c(40, 55), col = "lightblue", border = NA)
  plot(filters, border = "red", add = TRUE)  # plot main map with filter
  points(trends.m6$minlongs, trends.m6$minlats)

par(mfrow=c(1,1))

###############################################################################
# plot relationships between area, tice, and minlats FOR MAX_AREA-----

source("D:/Keith/capelin/2017-project/area-ice-function.R")
#m1
# plot tice against ice area
iceCorr.m1 <- iceScatterSummary(iceSum.m1, x = "area", x1 = "minlats", y = "tice", y1 = "minlats")
cowplot::plot_grid(iceCorr.m1$p1, iceCorr.m1$p2, iceCorr.m1$p3, labels = c("m1"), ncol=2)
ggsave("figs/8-relMax-m1.pdf", width=10, height=8, units="in")
iceCorr.m1.rsq <- lmoutputSummary(iceCorr.m1)
write.csv(iceCorr.m1.rsq, file = "output-processing/iceCorr-m1-rsq.csv", row.names=F, na="")

#m2
# plot tice against ice area
iceCorr.m2 <- iceScatterSummary(iceSum.m2, "area", "minlats", "tice", "minlats")
cowplot::plot_grid(iceCorr.m2$p1, iceCorr.m2$p2, iceCorr.m2$p3, labels = c("m2"), ncol=2)
ggsave("figs/9-relMax-m2.pdf", width=10, height=8, units="in")
iceCorr.m2.rsq <- lmoutputSummary(iceCorr.m2)
write.csv(iceCorr.m2.rsq, file = "output-processing/iceCorr-m2-rsq.csv", row.names=F, na="")

#m3
iceCorr.m3 <- iceScatterSummary(iceSum.m3, "area", "minlats", "tice", "minlats")
cowplot::plot_grid(iceCorr.m3$p1, iceCorr.m3$p2, iceCorr.m3$p3, labels = c("m3"), ncol=2)
ggsave("figs/10-relMax-m3.pdf", width=10, height=8, units="in")
iceCorr.m3.rsq <- lmoutputSummary(iceCorr.m3)
write.csv(iceCorr.m3.rsq, file = "output-processing/iceCorr-m3-rsq.csv", row.names=F, na="")

#m4
iceCorr.m4 <- iceScatterSummary(iceSum.m4, "area", "minlats", "tice", "minlats")
cowplot::plot_grid(iceCorr.m4$p1, iceCorr.m4$p2, iceCorr.m4$p3, labels = c("m4"), ncol=2)
ggsave("figs/11-relMax-m4.pdf", width=10, height=8, units="in")
iceCorr.m4.rsq <- lmoutputSummary(iceCorr.m4)
write.csv(iceCorr.m4.rsq, file = "output-processing/iceCorr-m4-rsq.csv", row.names=F, na="")

#m5
iceCorr.m5 <- iceScatterSummary(iceSum.m5, "area", "minlats", "tice", "minlats")
cowplot::plot_grid(iceCorr.m5$p1, iceCorr.m5$p2, iceCorr.m5$p3, labels = c("m5"), ncol=2)
ggsave("figs/12-relMax-m5.pdf", width=10, height=8, units="in")
iceCorr.m5.rsq <- lmoutputSummary(iceCorr.m5)
write.csv(iceCorr.m5.rsq, file = "output-processing/iceCorr-m5-rsq.csv", row.names=F, na="")

#m6
iceCorr.m6 <- iceScatterSummary(iceSum.m6, "area", "minlats", "tice", "minlats")
cowplot::plot_grid(iceCorr.m6$p1, iceCorr.m6$p2, iceCorr.m6$p3, labels = c("m6"), ncol=2)
ggsave("figs/13-relMax-m6.pdf", width=10, height=8, units="in")
iceCorr.m6.rsq <- lmoutputSummary(iceCorr.m6)
write.csv(iceCorr.m6.rsq, file = "output-processing/iceCorr-m6-rsq.csv", row.names=F, na="")


#####################################################################
# explore more rigorous method of determining relationship between tice and minlat-----
# calculate tice and year

str(trends.m1)

trends.m1$tice <- yday(trends.m1$date)
trends.m1$year <- year(trends.m1$date)

trends.m2$tice <- yday(trends.m2$date)
trends.m2$year <- year(trends.m2$date)

trends.m3$tice <- yday(trends.m3$date)
trends.m3$year <- year(trends.m3$date)

trends.m4$tice <- yday(trends.m4$date)
trends.m4$year <- year(trends.m4$date)

trends.m5$tice <- yday(trends.m5$date)
trends.m5$year <- year(trends.m5$date)

trends.m6$tice <- yday(trends.m6$date)
trends.m6$year <- year(trends.m6$date)

# Medians of all values
#####1991-----
#m1
sub1991.m1 <- iceMedian(trends.m1, "year < 1992", "tice < 150", iceSum.m1) # list of three dataframes
sub1991.m1.sum <- iceSummarylm(sub1991.m1) # lm Summaries of yearly median values
sub1991.m1.rsq <- lmoutputSummary(sub1991.m1.sum) # rsqared values only
write.csv(sub1991.m1.rsq, file = "output-processing/sub1991-m1-rsq.csv", row.names=F, na="")

# graph these for this model but only repeat if reasonable r-squared values appear
p1.m1  <- ggplot(data = sub1991.m1$mall, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~area, data=sub1991.m1$mall))

# plot tice against ice area
p2.m1 <- ggplot(data = sub1991.m1$mall, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)
summary(lm(tice~minlats, data=sub1991.m1$mall))

#plot ice area against minlats
p3.m1 <- ggplot(data = sub1991.m1$mall, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)
summary(lm(minlats~area, data=sub1991.m1$mall))

# plot of median values
cowplot::plot_grid(p1.m1, p2.m1, p3.m1, labels = c("medians-1969-1991-m1"), ncol=2)
ggsave("figs/14-relMeds1991.pdf", width=10, height=8, units="in")

# scatter plots of all dates
sub1991.m1.plot <- iceDateScatter(sub1991.m1)
sub1991.m1.plot$p1
ggsave("figs/15-minlatMedYear1991m1.pdf", width=10, height=8, units="in")
sub1991.m1.plot$p2
ggsave("figs/16-minlatMedViolin991m1.pdf", width=10, height=8, units="in")
sub1991.m1.plot$p3
ggsave("figs/17-areaMedViolin1991m1.pdf", width=10, height=8, units="in")

#m2
sub1991.m2 <- iceMedian(trends.m2, "year < 1992", "tice < 150", iceSum.m2)
sub1991.m2.sum <- iceSummarylm(sub1991.m2)
sub1991.m2.rsq <- lmoutputSummary(sub1991.m2.sum)
write.csv(sub1991.m2.rsq, file = "output-processing/sub1991-m2-rsq.csv", row.names=F, na="")

# scatter plot of date and minlats
sub1991.m2.plot <- iceDateScatter(sub1991.m2)
sub1991.m2.plot$p1
sub1991.m2.plot$p2
sub1991.m2.plot$p3

#m3
sub1991.m3 <- iceMedian(trends.m3, "year < 1992", "tice < 150", iceSum.m3)
sub1991.m3.sum <- iceSummarylm(sub1991.m3)
sub1991.m3.rsq <- lmoutputSummary(sub1991.m3.sum)
write.csv(sub1991.m3.rsq, file = "output-processing/sub1991-m3-rsq.csv", row.names=F, na="")

#m4
sub1991.m4 <- iceMedian(trends.m4, "year < 1992", "tice < 150", iceSum.m4)
sub1991.m4.sum <- iceSummarylm(sub1991.m4)
sub1991.m4.rsq <- lmoutputSummary(sub1991.m4.sum)
write.csv(sub1991.m4.rsq, file = "output-processing/sub1991-m4-rsq.csv", row.names=F, na="")

#m5
sub1991.m5 <- iceMedian(trends.m5, "year < 1992", "tice < 150", iceSum.m5)
sub1991.m5.sum <- iceSummarylm(sub1991.m5)
sub1991.m5.rsq <- lmoutputSummary(sub1991.m5.sum)
write.csv(sub1991.m5.rsq, file = "output-processing/sub1991-m5-rsq.csv", row.names=F, na="")

#m6
sub1991.m6 <- iceMedian(trends.m6, "year < 1992", "tice < 150", iceSum.m6)
sub1991.m6.sum <- iceSummarylm(sub1991.m6)
sub1991.m6.rsq <- lmoutputSummary(sub1991.m6.sum)
write.csv(sub1991.m6.rsq, file = "output-processing/sub1991-m6-rsq.csv", row.names=F, na="")

##############################################################
#####1992-2017------
##m1
sub2017.m1 <- iceMedian(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
sub2017.m1.sum <- iceSummarylm(sub2017.m1)
sub2017.m1.rsq <- lmoutputSummary(sub2017.m1.sum)
write.csv(sub2017.m1.rsq, file = "output-processing/sub2017-m1-rsq.csv", row.names=F, na="")

write.csv(trends.m1, file = "output-processing/trends_m1.csv", row.names=F, na="")

p1.m1  <- ggplot(data = sub2017.m1$mall, aes(x = area, y = tice)) + geom_point() + geom_smooth(method=lm)

# plot tice against ice area
p2.m1 <- ggplot(data = sub2017.m1$mall, aes(x = minlats, y = tice)) + geom_point() + geom_smooth(method=lm)

#plot ice area against minlats
p3.m1 <- ggplot(data = sub2017.m1$mall, aes(x = area, y = minlats)) + geom_point() + geom_smooth(method=lm)

# plot of median values
cowplot::plot_grid(p1.m1, p2.m1, p3.m1,  labels = c("medians-1992-2017-m1"), ncol=2)

ggsave("figs/18-relMeds2017m1.pdf", width=10, height=8, units="in")

sub2017.m1.plot <- iceDateScatter(sub2017.m1)
sub2017.m1.plot$p1
sub2017.m1.plot$p2
sub2017.m1.plot$p3

cowplot::plot_grid(p1.m1, p2.m1, p3.m1, labels = c("medians-m1"), ncol=2)

##m2
sub2017.m2 <- iceMedian(trends.m2, "year > 1991", "tice < 150", iceSum.m2)
sub2017.m2.sum <- iceSummarylm(sub2017.m2)
sub2017.m2.rsq <- lmoutputSummary(sub2017.m2.sum)
write.csv(sub2017.m2.rsq, file = "output-processing/sub2017-m2-rsq.csv", row.names=F, na="")

# scatter plot of date and minlats
sub2017.m2.plot <- iceDateScatter(sub2017.m2, d = "d")
sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

##m3
sub2017.m3 <- iceMedian(trends.m3, "year > 1991", "tice < 150", iceSum.m3)
sub2017.m3.sum <- iceSummarylm(sub2017.m3)
sub2017.m3.rsq <- lmoutputSummary(sub2017.m3.sum)
write.csv(sub2017.m3.rsq, file = "output-processing/sub2017-m3-rsq.csv", row.names=F, na="")

##m4
sub2017.m4 <- iceMedian(trends.m4, "year > 1991", "tice < 150", iceSum.m4)
sub2017.m4.sum <- iceSummarylm(sub2017.m4)
sub2017.m4.rsq <- lmoutputSummary(sub2017.m4.sum)
write.csv(sub2017.m4.rsq, file = "output-processing/sub2017-m4-rsq.csv", row.names=F, na="")

##m5
sub2017.m5 <- iceMedian(trends.m5, "year > 1991", "tice < 150", iceSum.m5)
sub2017.m5.sum <- iceSummarylm(sub2017.m5)
sub2017.m5.rsq <- lmoutputSummary(sub2017.m5.sum)
write.csv(sub2017.m5.rsq, file = "output-processing/sub2017-m5-rsq.csv", row.names=F, na="")

##m6
sub2017.m6 <- iceMedian(trends.m6, "year > 1991", "tice < 150", iceSum.m6)
sub2017.m6.sum <- iceSummarylm(sub2017.m6)
sub2017.m6.rsq <- lmoutputSummary(sub2017.m6.sum)
write.csv(sub2017.m6.rsq, file = "output-processing/sub2017-m6-rsq.csv", row.names=F, na="")

###############################################
###############################################
# Median of the top 5 values
count.year <- trends.m1 %>%
  group_by(year) %>%
  count()
View(count.year)

###1992-2017------
##m1
iceMedD5p2017.m1 <- iceMedianD5(trends.m1, "year > 1991", "tice < 150", iceSum.m1)
iceMedD5p2017.m1.sum <-iceSummarylm(iceMedD5p2017.m1, med = "d5", med1 = "d5")
iceMedD5p2017.m1.rsq <- lmoutputSummary(iceMedD5p2017.m1.sum)
write.csv(iceMedD5p2017.m1.rsq, file = "output-processing/iceMedD5p2017-m1-rsq.csv", row.names=F, na="")

# figures
plotMedD5p2017.m1 <- iceScatterSummary(iceMedD5p2017.m1$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p2017.m1$p1, plotMedD5p2017.m1$p2, plotMedD5p2017.m1$p3,  labels = c("medians-m1"), ncol=2)
ggsave("figs/19-relD5Med2017m1.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
#windows()
sub2017.m1.plot <- iceDateScatter(iceMedD5p2017.m1, d= "d5")
sub2017.m1.plot$p1
sub2017.m1.plot$p2
sub2017.m1.plot$p3

###1992-2017------
##m2
iceMedD5p2017.m2 <- iceMedianD5(trends.m2, "year > 1991", "tice < 150", iceSum.m2)
iceMedD5p2017.m2.sum <- iceSummarylm(iceMedD5p2017.m2, med = "d5", med1 = "d5")
iceMedD5p2017.m2.rsq <- lmoutputSummary(iceMedD5p2017.m2.sum)
write.csv(iceMedD5p2017.m2.rsq, file = "output-processing/iceMedD5p2017-m2-rsq.csv", row.names=F, na="")

# figures
plotMedD5p2017.m2 <- iceScatterSummary(iceMedD5p2017.m2$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p2017.m2$p1, plotMedD5p2017.m2$p2, plotMedD5p2017.m2$p3,  labels = c("medians-m1"), ncol=2)
ggsave("figs/20-relD5Med2017m2.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
sub2017.m2.plot <- iceDateScatter(iceMedD5p2017.m2, d= "d5")
sub2017.m2.plot$p1
sub2017.m2.plot$p2
sub2017.m2.plot$p3

###1992-2017------
##m3
iceMedD5p2017.m3 <- iceMedianD5(trends.m3, "year > 1991", "tice < 150", iceSum.m3)
iceMedD5p2017.m3.sum <- iceSummarylm(iceMedD5p2017.m3, med = "d5", med1 = "d5")
iceMedD5p2017.m3.rsq <- lmoutputSummary(iceMedD5p2017.m3.sum)
write.csv(iceMedD5p2017.m3.rsq, file = "output-processing/iceMedD5p2017-m3-rsq.csv", row.names=F, na="")

# figures
plotMedD5p2017.m3 <- iceScatterSummary(iceMedD5p2017.m3$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p2017.m3$p1, plotMedD5p2017.m3$p2, plotMedD5p2017.m3$p3,  labels = c("medians-m1"), ncol=2)
ggsave("figs/21-relD5Med2017m3.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
sub2017.m3.plot <- iceDateScatter(iceMedD5p2017.m3, d= "d5")
sub2017.m3.plot$p1
sub2017.m3.plot$p2
sub2017.m3.plot$p3

###1992-2017------
##m4
iceMedD5p2017.m4 <- iceMedianD5(trends.m4, "year > 1991", "tice < 150", iceSum.m4)
iceMedD5p2017.m4.sum <- iceSummarylm(iceMedD5p2017.m4, med = "d5", med1 = "d5")
iceMedD5p2017.m4.rsq <- lmoutputSummary(iceMedD5p2017.m4.sum)
write.csv(iceMedD5p2017.m4.rsq, file = "output-processing/iceMedD5p2017-m4-rsq.csv", row.names=F, na="")

# figures
plotMedD5p2017.m4 <- iceScatterSummary(iceMedD5p2017.m4$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p2017.m4$p1, plotMedD5p2017.m4$p2, plotMedD5p2017.m4$p3,  labels = c("medians-m1"), ncol=2)
ggsave("figs/22-relD5Med2017m4.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
sub2017.m4.plot <- iceDateScatter(iceMedD5p2017.m4, d= "d5")
sub2017.m4.plot$p1
sub2017.m4.plot$p2
sub2017.m4.plot$p3

##m5
iceMedD5p2017.m5 <- iceMedianD5(trends.m5, "year > 1991", "tice < 150", iceSum.m5)
iceMedD5p2017.m5.sum <- iceSummarylm(iceMedD5p2017.m5, med = "d5", med1 = "d5")
iceMedD5p2017.m5.rsq <- lmoutputSummary(iceMedD5p2017.m5.sum)
write.csv(iceMedD5p2017.m5.rsq, file = "output-processing/iceMedD5p2017-m5-rsq.csv", row.names=F, na="")

##m6
iceMedD5p2017.m6 <- iceMedianD5(trends.m6, "year > 1991", "tice < 150", iceSum.m6)
iceMedD5p2017.m6.sum <- iceSummarylm(iceMedD5p2017.m6, med = "d5", med1 = "d5")
iceMedD5p2017.m6.rsq <- lmoutputSummary(iceMedD5p2017.m6.sum)
write.csv(iceMedD5p2017.m6.rsq, file = "output-processing/iceMedD5p2017-m6-rsq.csv", row.names=F, na="")

##########################################################################################
###pre-1992------
##m1
iceMedD5p1992.m1 <- iceMedianD5(trends.m1, "year <= 1991", "tice < 150", iceSum.m1)
iceMedD5p1992.m1.sum <- iceSummarylm(iceMedD5p1992.m1, med = "d5", med1 = "d5")
iceMedD5p1992.m1.sum.rsq <- lmoutputSummary(iceMedD5p1992.m1.sum)
write.csv(iceMedD5p1992.m1.sum.rsq, file = "output-processing/iceMedD5p1992-m1-sum.csv", row.names=F, na="")

# figures
plotMedD5p1992.m1 <- iceScatterSummary(iceMedD5p1992.m1$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p1992.m1$p1, plotMedD5p1992.m1$p2, plotMedD5p1992.m1$p3,  labels = c("d5medians1960-1992-m1"), ncol=2)
ggsave("figs/23-relD5Med1992m1.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
sub1992.m1.plot <- iceDateScatter(iceMedD5p1992.m1, d = "d5")
sub1992.m1.plot$p1
sub1992.m1.plot$p2
sub1992.m1.plot$p3

###1992-2017------
##m2
iceMedD5p1992.m2 <- iceMedianD5(trends.m2, "year <= 1991", "tice < 150", iceSum.m2)
iceMedD5p1992.m2.sum <- iceSummarylm(iceMedD5p1992.m2, med = "d5", med1 = "d5")
iceMedD5p1992.m2.sum.rsq <- lmoutputSummary(iceMedD5p1992.m2.sum)
write.csv(iceMedD5p1992.m2.sum.rsq, file = "output-processing/iceMedD5p1992-m2-sum.csv", row.names=F, na="")

# figures
plotMedD5p1992.m2 <- iceScatterSummary(iceMedD5p1992.m2$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
cowplot::plot_grid(plotMedD5p1992.m2$p1, plotMedD5p1992.m2$p2, plotMedD5p1992.m2$p3,  labels = c("d5medians1960-1992-m1"), ncol=2)
ggsave("figs/24-relD5Med1992m2.pdf", width=10, height=8, units="in")

# scatter plot of date and minlats
sub1992.m2.plot <- iceDateScatter(iceMedD5p1992.m2, d= "d5")
sub1992.m2.plot$p1
sub1992.m2.plot$p2
sub1992.m2.plot$p3

###1992-2017------
##m3
iceMedD5p1992.m3 <- iceMedianD5(trends.m3, "year <= 1991", "tice < 150", iceSum.m3)
iceMedD5p1992.m3.sum <- iceSummarylm(iceMedD5p1992.m3, med = "d5", med1 = "d5")
iceMedD5p1992.m3.sum.rsq <- lmoutputSummary(iceMedD5p1992.m3.sum)
write.csv(iceMedD5p1992.m3.sum.rsq, file = "output-processing/iceMedD5p1992-m3-sum.csv", row.names=F, na="")

# figures
plotMedD5p1992.m3 <- iceScatterSummary(iceMedD5p1992.m3$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
multiplot(plotMedD5p1992.m3$p1, plotMedD5p1992.m3$p3, plotMedD5p1992.m2$p2, cols=2)

# scatter plot of date and minlats
sub1992.m3.plot <- iceDateScatter(iceMedD5p1992.m3, d= "d5")
sub1992.m3.plot$p1
sub1992.m3.plot$p2
sub1992.m3.plot$p3

###1992-2017------
##m4
iceMedD5p1992.m4 <- iceMedianD5(trends.m4, "year <= 1991", "tice < 150", iceSum.m4)
iceMedD5p1992.m4.sum <- iceSummarylm(iceMedD5p1992.m4, med = "d5", med1 = "d5")
iceMedD5p1992.m4.sum.rsq <- lmoutputSummary(iceMedD5p1992.m4.sum)
write.csv(iceMedD5p1992.m4.sum.rsq, file = "output-processing/iceMedD5p1992-m4-sum.csv", row.names=F, na="")

# figures
plotMedD5p1992.m4 <- iceScatterSummary(iceMedD5p1992.m4$mall, "d5area", "d5minlats", "d5tice", "d5minlats")
multiplot(plotMedD5p1992.m4$p1, plotMedD5p1992.m4$p3, plotMedD5p1992.m4$p3, cols=2)

# scatter plot of date and minlats
sub1992.m4.plot <- iceDateScatter(iceMedD5p1992.m4, d= "d5")
sub1992.m4.plot$p1
sub1992.m4.plot$p2
sub1992.m4.plot$p3

##m5
iceMedD5p1992.m5 <- iceMedianD5(trends.m5, "year <= 1991", "tice < 150", iceSum.m5)
iceMedD5p1992.m5.sum <- iceSummarylm(iceMedD5p1992.m5, med = "d5", med1 = "d5")
iceMedD5p1992.m5.sum.rsq <- lmoutputSummary(iceMedD5p1992.m5.sum)
write.csv(iceMedD5p1992.m5.sum.rsq, file = "output-processing/iceMedD5p1992-m5-sum.csv", row.names=F, na="")

##m6
iceMedD5p1992.m6 <- iceMedianD5(trends.m6, "year <= 1991", "tice < 150", iceSum.m6)
iceMedD5p1992.m6.sum <- iceSummarylm(iceMedD5p1992.m6, med = "d5", med1 = "d5")
iceMedD5p1992.m6.sum.rsq <- lmoutputSummary(iceMedD5p1992.m6.sum)
write.csv(iceMedD5p1992.m6.sum.rsq, file = "output-processing/iceMedD5p1992-m6-sum.csv", row.names=F, na="")

###########################################################
# extra graphs
iceMed.m1 <- iceMedian(trends.m1, "year < 2018", "tice < 150", iceSum.m1)
iceMed.m2 <- iceMedian(trends.m2, "year < 2018", "tice < 150", iceSum.m2)
iceMed.m3 <- iceMedian(trends.m3, "year < 2018", "tice < 150", iceSum.m3)
iceMed.m4 <- iceMedian(trends.m4, "year < 2018", "tice < 150", iceSum.m4)
iceMed.m5 <- iceMedian(trends.m5, "year < 2018", "tice < 150", iceSum.m5)
iceMed.m6 <- iceMedian(trends.m6, "year < 2018", "tice < 150", iceSum.m6)

iceMedD5.m1 <- iceMedianD5(trends.m1, "year < 2018", "tice < 150", iceSum.m1)
iceMedD5.m2 <- iceMedianD5(trends.m2, "year < 2018", "tice < 150", iceSum.m2)
iceMedD5.m3 <- iceMedianD5(trends.m3, "year < 2018", "tice < 150", iceSum.m3)
iceMedD5.m4 <- iceMedianD5(trends.m4, "year < 2018", "tice < 150", iceSum.m4)
iceMedD5.m5 <- iceMedianD5(trends.m5, "year < 2018", "tice < 150", iceSum.m5)
iceMedD5.m6 <- iceMedianD5(trends.m6, "year < 2018", "tice < 150", iceSum.m6)

source("D:/Keith/capelin/2017-project/area-ice-function.R")
#m1
atScat.m1 <- area_ticeScatter(df1 = iceSum.m1, df2 = iceMed.m1, df3 = iceMedD5.m1, x = "tice", y1 = "area", y2 = "darea",   y3 = "d5area", sub_val = "m1")
ggsave("figs/area_v_tice/41-areaticeScat-m1.pdf", width=10, height=8, units="in")

#m2
atScat.m2 <- area_ticeScatter(df1 = iceSum.m2, df2 = iceMed.m2, df3 = iceMedD5.m2, x="tice", y="area", y2 = "darea",   y3 = "d5area", sub_val = "m2")
ggsave("figs/area_v_tice/42-areaticeScat-m2.pdf", width=10, height=8, units="in")

#m3
atScat.m3 <- area_ticeScatter(df1 = iceSum.m3, df2 = iceMed.m3, df3 = iceMedD5.m3, x = "tice", y1 = "area", y2 = "darea",   y3 = "d5area", sub_val = "m3")
ggsave("figs/area_v_tice/43-areaticeScat-m3.pdf", width=10, height=8, units="in")

#m4
atScat.m4 <- area_ticeScatter(iceSum.m4, iceMed.m4, iceMedD5.m4, x = "tice", y1 = "area", y2 = "darea",   y3 = "d5area", sub_val = "m4")
ggsave("figs/area_v_tice/44-areaticeScat-m4.pdf", width=10, height=8, units="in")

 #m5
atScat.m5 <- area_ticeScatter(iceSum.m5, iceMed.m5, iceMedD5.m5, x = "tice", y1 = "area", y2 = "darea",   y3 = "d5area", sub_val = "m5")
ggsave("figs/area_v_tice/45-areaticeScat-m5.pdf", width=10, height=8, units="in")

#m6
atScat.m6 <- area_ticeScatter(iceSum.m6, iceMed.m6, iceMedD5.m6, x = "tice", y1 = "area", y2 = "darea",   y3 = "d5area", sub_val = "m6")
ggsave("figs/area_v_tice/46-areaticeScat-m6.pdf", width=10, height=8, units="in")

source("D:/Keith/capelin/2017-project/area-ice-function.R")

## Correlations
area_tice_corr_1 <- area_ticeRsq(df1 = iceSum.m1, df2 = iceMed.m1, df3 = iceMedD5.m1, x = "tice", y = "area")
write.csv(area_tice_corr_1, file = "output-processing/area_tice_corr_1.csv", row.names=F, na="")

area_tice_corr_2 <- area_ticeRsq(df1 = iceSum.m2, df2 = iceMed.m2, df3 = iceMedD5.m2, x = "tice", y = "area")
write.csv(area_tice_corr_2, file = "output-processing/area_tice_corr_2.csv", row.names=F, na="")

area_tice_corr_3 <- area_ticeRsq(df1 = iceSum.m3, df2 = iceMed.m3, df3 = iceMedD5.m3, x = "tice", y = "area")
write.csv(area_tice_corr_3, file = "output-processing/area_tice_corr_3.csv", row.names=F, na="")

area_tice_corr_4 <- area_ticeRsq(df1 = iceSum.m4, df2 = iceMed.m4, df3 = iceMedD5.m4, x = "tice", y = "area")
write.csv(area_tice_corr_4, file = "output-processing/area_tice_corr_4.csv", row.names=F, na="")

area_tice_corr_5 <- area_ticeRsq(df1 = iceSum.m5, df2 = iceMed.m5, df3 = iceMedD5.m5, x = "tice", y = "area")
write.csv(area_tice_corr_5, file = "output-processing/area_tice_corr_5.csv", row.names=F, na="")

area_tice_corr_6 <- area_ticeRsq(df1 = iceSum.m6, df2 = iceMed.m6, df3 = iceMedD5.m6, x = "tice", y = "area")
write.csv(area_tice_corr_6, file = "output-processing/area_tice_corr_6.csv", row.names=F, na="")

####################################################################################### OUTPUT of derived values for each year
#######Max values-----------

# need to change all fo these values from here to below
iceSum.m1 <- iceSum.m1[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m1, file = "output-processing/capelin-m1.csv", row.names=F, na="")

iceSum.m2 <- iceSum.m2[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m2, file = "output-processing/capelin-m2.csv", row.names=F, na="")

iceSum.m3 <- iceSum.m3[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m3, file = "output-processing/capelin-m3.csv", row.names=F, na="")

iceSum.m4 <- iceSum.m4[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m4, file = "output-processing/capelin-m4.csv", row.names=F, na="")

iceSum.m5 <- iceSum.m5[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m5, file = "output-processing/capelin-m5.csv", row.names=F, na="")

iceSum.m6 <- iceSum.m6[c("year", "area", "minlats", "tice")]
write.csv(iceSum.m6, file = "output-processing/capelin-m6.csv", row.names=F, na="")

#######################################################################################
#######Median values-----------
## < 1992
sub1991.m1a <- sub1991.m1$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m1a, file = "output-processing/sub1991-m1.csv", row.names=F, na="")

sub1991.m2a <- sub1991.m2$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m2a, file = "output-processing/sub1991-m2.csv", row.names=F, na="")

sub1991.m3a <- sub1991.m3$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m3a, file = "output-processing/sub1991-m3.csv", row.names=F, na="")

sub1991.m4a <- sub1991.m4$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m4a, file = "output-processing/sub1991-m4.csv", row.names=F, na="")

sub1991.m5a <- sub1991.m5$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m5a, file = "output-processing/sub1991-m5.csv", row.names=F, na="")

sub1991.m6a <- sub1991.m6$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub1991.m6a, file = "output-processing/sub1991-m6.csv", row.names=F, na="")

## > 1992
sub2017.m1a <- sub2017.m1$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m1a, file = "output-processing/sub2017-m1.csv", row.names=F, na="")

sub2017.m2a <- sub2017.m2$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m2a, file = "output-processing/sub2017-m2.csv", row.names=F, na="")

sub2017.m3a <- sub2017.m3$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m3a, file = "output-processing/sub2017-m3.csv", row.names=F, na="")

sub2017.m4a <- sub2017.m4$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m4a, file = "output-processing/sub2017-m4.csv", row.names=F, na="")

sub2017.m5a <- sub2017.m5$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m5a, file = "output-processing/sub2017-m5.csv", row.names=F, na="")

sub2017.m6a <- sub2017.m6$mall[c("year", "darea", "dminlats", "dtice", "tice")]
write.csv(sub2017.m6a, file = "output-processing/sub2017-m6.csv", row.names=F, na="")

#######D5Median values-----------
# <1992
iceMedD5p1992.m1a <- iceMedD5p1992.m1$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m1a, file = "output-processing/iceMedD5p1992-m1a.csv", row.names=F, na="")

iceMedD5p1992.m2a <- iceMedD5p1992.m2$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m2a, file = "output-processing/iceMedD5p1992-m2a.csv", row.names=F, na="")

iceMedD5p1992.m3a <- iceMedD5p1992.m3$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m3a, file = "output-processing/iceMedD5p1992-m3a.csv", row.names=F, na="")

iceMedD5p1992.m4a <- iceMedD5p1992.m4$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m4a, file = "output-processing/iceMedD5p1992-m4a.csv", row.names=F, na="")

iceMedD5p1992.m5a <- iceMedD5p1992.m5$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m5a, file = "output-processing/iceMedD5p1992-m5a.csv", row.names=F, na="")

iceMedD5p1992.m6a <- iceMedD5p1992.m6$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p1992.m6a, file = "output-processing/iceMedD5p1992-m6a.csv", row.names=F, na="")

# >1992
iceMedD5p2017.m1a <- iceMedD5p2017.m1$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m1a, file = "output-processing/iceMedD5p2017-m1a.csv", row.names=F, na="")

iceMedD5p2017.m2a <- iceMedD5p2017.m2$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m2a, file = "output-processing/iceMedD5p2017-m2a.csv", row.names=F, na="")

iceMedD5p2017.m3a <- iceMedD5p2017.m3$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m3a, file = "output-processing/iceMedD5p2017-m3a.csv", row.names=F, na="")

iceMedD5p2017.m4a <- iceMedD5p2017.m4$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m4a, file = "output-processing/iceMedD5p2017-m4a.csv", row.names=F, na="")

iceMedD5p2017.m5a <- iceMedD5p2017.m5$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m5a, file = "output-processing/iceMedD5p2017-m5a.csv", row.names=F, na="")

iceMedD5p2017.m6a <- iceMedD5p2017.m6$mall[c("year", "d5area", "d5minlats", "d5tice", "tice")]
write.csv(iceMedD5p2017.m6a, file = "output-processing/iceMedD5p2017-m6a.csv", row.names=F, na="")

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

trends.m6 %>%
     filter(year > 2002 & tice < 200) %>%
     ggplot() + 
     geom_line(aes(x=tice, y=minlats)) + 
     facet_wrap(~year)
