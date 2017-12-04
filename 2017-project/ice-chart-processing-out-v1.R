################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# The purpose of this file is to:
#1) calculate non-Gaussian curves that Paul examined.  File deprecated for now.

## Area analysis ---------------------------------------------------------------

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
library(data.table)
library(ggplot2)

load("ice_trends_2017.Rdata") # this is date, area, and volume
trends <- data.table(trends) # shows head and tail of "trends"
trends$doy <- as.POSIXlt(trends$date)$yday
trends$doy <- ifelse(trends$doy > 305, trends$doy - 365, trends$doy) # force doy to negative if between 305 to 365 (nov-dec)
trends$day <- as.integer(format(trends$date, "%d")) 
trends$month <- as.integer(format(trends$date, "%m"))
trends$year <- as.integer(format(trends$date, "%Y"))
trends$year <- ifelse(trends$month %in% 11:12, trends$year + 1, trends$year) # if nov or dec, plus 1 to year
trends$psudo.date <- as.Date(ISOdate(ifelse(trends$month %in% 11:12, 2015, 2016), 
                                     trends$month, trends$day)) # for plotting
lookAt(trends)


############create a bunch of figures
# area v. date
ggplot() + geom_line(aes(x = date, y = area), trends, size = 1)

# graph in three year chunks
p1 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 1990:1992], size = 1)
p1
p2 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 1995:1997], size = 1)

p3 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 2006:2008], size = 1)

p4 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 2007:2010], size = 1)

multiplot(p1, p3, p2, p4, cols=2)

# facet graph over 11 years
ggplot() + geom_line(aes(x = psudo.date, y = area), 
                     trends[year %in% 2005:2016], size = 1) + facet_wrap(~ year)


## Not sure if volume is appropriate...too big of an approximation maybe?

#########################################################################
# this prints out ~ graphs plotting the predicted values of two differnt models as well as the data from 1970-2016!!!!!!!!!!!!!!!!!
# a = area/1000
# t = doy
# models are mathematical rather than biological - meant to be like growing degree days but for ice. Can determine ice retreat by calculating the max of each curve.


years <- unique(trends$year) # create a vector of years
adv.mag <- ret.mag <- amax <- tstart <- tend <- tmax <-
  raw.amax <- raw.tmax <- rep(NA, length(years)) # create a series of vectors filled with NA

modelGraphs(years) # makes all the graphs

# need something to calculate max of each curve and return the date
# this is in the code - just need to return it
# 
###################################################
# the following are line graphs that show various measures of ice at 3 levels, e.g. changes in the timing of the start, max, and end of ice.  Important exploratory work.

pars <- data.table(year = years, adv.mag, ret.mag, amax, tstart, 
                   tend, tmax, raw.amax, raw.tmax)
pars$tdur <- pars$tend - pars$tstart
pars$ret.dur <- pars$tend - pars$tmax
pars$adv.dur <- pars$tmax - pars$tstart
pars$mag <- pars$adv.mag + pars$ret.mag

save(pars, file = "ice_pars.Rdata")


## PLot make Area v year
ggplot(pars, aes (x = year, y = amax)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Plot of various measures of "timing" by year
tice <- melt(pars, id.vars = "year", measure.vars = c("tstart", "tmax", "tend"),
             variable.name = "timing", value.name = "doy")
ggplot(tice, aes(x = year, y = doy, col = timing)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Duration
dice <- melt(pars, id.vars = "year", measure.vars = c("tdur", "adv.dur", "ret.dur"),
             variable.name = "type", value.name = "duration")
ggplot(dice, aes(x = year, y = duration, col = type)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Maginitude
mice <- melt(pars, id.vars = "year", measure.vars = c("mag", "adv.mag", "ret.mag"),
             variable.name = "type", value.name = "magnitude")
ggplot(mice, aes(x = year, y = magnitude, col = type)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

