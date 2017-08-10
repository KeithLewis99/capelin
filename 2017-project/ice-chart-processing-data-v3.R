################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# The purpose of this file is to:
  #1) Set-up: bring in appropriate source code, libraries, and set up directories
  #2) Metadata: set up metadata for Environment Canada (EC) ice codes
  #3) download the EC ice map files
  #4) convert e00 files
  #5) Cross-check coastline (map projection)
  #6) Make filter polygons
  #7) Area and volume calculations

# Main issues (see README):
  #1) Few errors in downloading files from the internet or converting e00 to SpatialPolygonDataframe
  #2) Too many connections open when converting files
  #3) Lots of long loops that could be made into functions to simplify code. Done.
  #4) can date computations be simplified with lubridate?
  #5) circular problem with ice_trends data. See explanation below
  #6) how to get into archive so that we can see what has been missed and why. Can get into archive.
  #7) make sure that we can extract the same values as in Ale's paper to ensure this is accurately working (maxIce, tice)

# rest to output file

# http://ice-glaces.ec.gc.ca//www_archive/AOI_12/Coverages/rgc_a12_19690117_XXXXXX.e00
# http://ice-glaces.ec.gc.ca//www_archive/AOI_12/Coverages/rgc_a12_19690124_XXXXXX.e00
# http://ice-glaces.ec.gc.ca//www_archive/AOI_12/Coverages/rgc_a12_19770213_XXXXXX.e00
# http://ice-glaces.ec.gc.ca//www_archive/AOI_12/Coverages/rgc_a12_20050606_EXPREC.e00
# http://ice-glaces.ec.gc.ca//www_archive/AOI_12/Coverages/rgc_a12_20160509_CEXPREC.e00


rm(list=ls()) # Clear the workspace

# Files downloaded from http://iceweb1.cis.ec.gc.ca/Archive/page1.xhtml?grp=Guest&mn=&lang=en
## Set-up ----------------------------------------------------------------------

# devtools::install_github("eblondel/cleangeo") # for cleaning SpatialPolygon objects; also on CRAN
options(stringsAsFactors = FALSE)

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v2a.R")

#setwd("C:/Users/Paul/Documents/DFO/ice")
setwd("D:/Keith/capelin/2017-project")
library(RArcInfo)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cleangeo)
library(dplyr)
library(tidyr)
library(broom)
#library(doParallel)

## create dirs for data storage
if(!dir.exists("avc_data")) dir.create("avc_data")
if(!dir.exists("e00_data")) dir.create("e00_data")
if(!dir.exists("sp_data")) dir.create("sp_data")


## Metadata --------------------------------------------------------------------

## stage codes (https://ec.gc.ca/glaces-ice/default.asp?lang=En&n=D5F7EA14-1&offset=1&toc=show)
stages <- fread("description	thickness	code
                New ice	< 10 centimetres	1
                Nilas, Ice rind	< 10 centimetres	2
                Young Ice	10 - 30 centimetres	3
                Grey Ice	10 - 15 centimetres	4
                Grey-white ice	15 - 30 centimetres	5
                First-year ice	>= 30 centimetres	6
                Thin first-year ice	30 - 70 centimetres	7
                First stage thin first-year	30 - 50 centimetres	8
                Second stage thin first-year	50 - 70 centimetres	9
                Medium first-year ice	70 - 120 centimetres	1?
                Thick first-year ice	> 120 centimetres	4?")
stages$thickness <- gsub(" centimetres", "", stages$thickness) #get rid of centimeters in stages
trange <- strsplit(stages$thickness, " - ") # create a list of thickness values in stages
stages$thickness <- sapply(seq_along(trange), function(i) {
  if(length(trange[[i]]) == 2) {
    median(as.numeric(trange[[i]]))
  } else {
    as.numeric(gsub("[^0-9]", "", trange[[i]]))
  }
}) # creates a single numeric value based on thickness (mean for ranges, same value for <>)


## Download e00 files ----------------------------------------------------------
# this is very robust.  Seems to be only 1 error in almost 1500 files - however, if e00 files are deleted, the code does not always seem to look backwards.  Use with caution and check if there is reason to believe an e00 file has been deleted.

downloaded <- list.files("e00_data/", pattern = ".e00") # this code assumes at least one e00 file in directory, i.e., seed the directory
downloaded <- downloaded[downloaded != "test.e00"]
downloaded <- strptime(downloaded, "%Y%m%d.e00")
dates <- seq(ISOdate(1969, 1, 17), ISOdate(2017, 7, 11), by = "day") # update the second ISOdate to present as desired - however, this does not work well if not updated: also, it reports a numer of erros, one for each day since last e00 file - these should be checked but can generally be ignored.

dates <- dates[format(dates, "%m") %in% c("11", "12", "01", "02", "03", "04", "05", "06", "07")]

####################################################################
#############STOP
####################################################################
#run one of next two lines
dates <- dates[difftime(max(downloaded), dates, units = "days") < -1] #this line is for simply updating - assumes all past files are downloaded and only looks for new files on EC website past the terminal data supplied in line 92 above
#dates <- dates[!as.Date(dates) %in% as.Date(downloaded)] #this line is for downloading all data from the EC website - takes A LONG TIME!!!!!
dates <- format(dates, "%Y%m%d")
dates


# downloads all files not already in e00_data folder
# error messages indicate maps that have not been or may never be created.  Can be ignored. See previous comment
e00Download(dates)

## e00 file conversion to Spatial Polygons Dataframe---------------------------------------------------------
# e00 to avc_data (coverages) and then to SpatialPolygonsDAtaframe in sp_data

## note: had to run this across multiple sessions...too many open files.
##       don't know how to close the connections open by RArcInfo..grrrrr

# this is just to refersh the download file to ensure that all dates accounted for
downloaded <- list.files("e00_data/", pattern = ".e00") # this code assumes at least one e00 file in directory, i.e., seed the directory
downloaded <- downloaded[downloaded != "test.e00"]
downloaded <- strptime(downloaded, "%Y%m%d.e00")
converted <- list.files("sp_data/", pattern = ".Rdata")
converted <- strptime(converted, "%Y%m%d.Rdata")
dates1 <- downloaded[!downloaded %in% converted] # files left to convert ###########not sure this is working right
dates1 <- format(dates1, "%Y%m%d")
dates1

# helper function to minimize the code
e00_to_SpatialPolygonDataframe(dates1)

# note that there is no corresponding e00 file on the EC website.  PR suspects that it may be hidden.  More imortantly, there is one from 20070507 (one day before).  Eliminate file
file.remove("e00_data/20070508.e00")
file.remove("sp_data/20070508.Rdata")
unlink("avc_data/20070508")


# refresh converted and calculate percentage of files that did not convert to a SPDF
converted <- list.files("sp_data/", pattern = ".Rdata")
converted <- strptime(converted, "%Y%m%d.Rdata")
dates2 <- downloaded[!downloaded %in% converted] # files left to convert
length(dates2)/length(converted) # conversion didn't work out for a small number of e00 files


## Cross-check coastline (map projection) --------------------------------------

# ## load trusted map data to cross-check coastline against
# can <- getData("GADM", country = "CAN", level = 1)
# nl <- can[can$NAME_1 == "Newfoundland and Labrador", ]
# plot(nl)
# nf <- crop(nl, extent(c(-60.12000, -51.24861, 46.06806, 51.94517))) # Newfoundland
# plot(nf)
# bell <- crop(nl, extent(c(-53.03932, -52.89715, 47.58573, 47.66841))) # Bell Island
# plot(bell)
# 
# ## cycle through a subset of the converted ice charts and Bell Island coastline
# converted <- list.files("sp_data/", pattern = ".Rdata")
# dates <- strptime(converted, "%Y%m%d.Rdata")
# #dates <- dates[round(seq(1, length(dates), length.out = 20))] # size up a sub-sample of the files
# yrs <- as.numeric(format(dates, "%Y"))
# #dates <- dates[2005 <= yrs & yrs <= 2006]
# dates <- dates[yrs == 2007]
# for(i in seq_along(dates)) {
#   load(format(dates[i], "sp_data/%Y%m%d.Rdata"))
#   chart.nl <- ice[ice$A_LEGEND == "Land", ]
#   chart.nl <- spTransform(chart.nl, CRS(proj4string(nl)))
#   par(mfcol = c(1, 2), mar = c(3, 3, 0, 0), oma = c(0, 0, 3, 1))
#   plot(nf, col = "grey", border = NA, axes = TRUE)
#   plot(chart.nl, border = "black", add = TRUE)
#   box()
#   plot(bell, col = "grey", border = NA, axes = TRUE)
#   plot(chart.nl, border = "black", add = TRUE)
#   title(main = dates[i], outer = TRUE)
#   box()
# }
# # looks like they improved their map data (and/or changed projections) 20070716...


## Make filter polygons --------------------------------------------------------
# purpose here is too .....filter certain areas where ice in bays and inlets (e.g. Lake Melville) might throw off the calculations

## Any polygon will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84")) # this is a SpatialPolygondDataFrame
water <- gUnaryUnion(water) # this makes a Spatial Polygon

## Lake Melvelle
plot(water, xlim = c(-62, -56), ylim = c(53.05, 54.52), col = "lightblue", border = NA)
x <- c(-58.26572, -58.18660, -57.68811, -59.42886,
       -61.16169, -61.39115, -59.93525, -58.70091) # captured using locator function
y <- c(54.25903, 54.21695, 54.12346, 52.94540, 
       52.92670, 53.93646, 54.45069, 54.58626)
p1 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
proj4string(p1) <- proj4string(water) # set coordinates reference system
melville <- gIntersection(water, p1, byid = TRUE) # intesection of water (ocean surronding NL) and Lake Melville
melville <- gBuffer(melville, width = 0.1) # add a small buffer
plot(melville, border = "red", add = TRUE) # plots Lake Melville with border

# plot North Atlantic (Gulf of Maine to Labrador Sea)
plot(water, col = "lightblue", border = NA) 
north55 <- extent(water) # this makes the extent of water an object
north55@ymin <- 55      # this uses the extent to make a rectangle that starts at 55 lat
north55 <- as(north55, "SpatialPolygons") # makes north55 a SpatialPolygon
proj4string(north55) <- proj4string(water)
plot(north55, border = "red", add = TRUE) # plot red rectangle over norther Lab Sea (Latitude 55 deg???)


# not needed
# intesection of water (ocean surronding NL) and Lake Melville
#eastcoast <- gIntersection(water, p1, byid = TRUE) 
#eastcoast2 <- gIntersection(water, p2, byid = TRUE) 
#eastcoast <- gBuffer(eastcoast, width = 0.1) # add a small buffer
#plot(eastcoast, border = "red", add = TRUE) # plots Lake Melville with border

#Eastern Newfoundland
plot(water, xlim = c(-53.25, -53.23), ylim = c(46.75, 49.5), col = "lightblue", border = NA)
x <- c(-52.79640, -52.91794, -53.07998, -53.24203, -53.52560, -53.93072,
       -53.99148, -53.98136, -54.13327, -54.30545, -54.02187, -54.16366, 
       -54.32570, -53.64714, -53.09011, -52.91794, -52.79640) 

x<- c(-53.49140, -53.01775, -52.88416, -52.71413, -53.57641, -54.20794, -55.33742, -54.54800, -54.03792)# captured using locator()
y <- c(47.75702, 47.52098, 47.39285, 47.51424, 47.50075, 47.72330,
       47.92562, 48.21561, 48.33025, 48.40444, 48.60001, 48.72814,
       48.78884, 49.20022, 48.64047, 48.11445, 47.76376)
y <- c(49.25103, 48.65924, 48.12419, 47.75128, 46.66497, 46.81900, 47.00545, 47.71885, 49.42938)
p2 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
proj4string(p2) <- proj4string(water)
eastNF <- gIntersection(water, p2, byid = TRUE) # intesection of water (ocean surronding NL) and Lake Melville
eastNF <- gBuffer(eastNF, width = 0.1) # add a small buffer
plot(eastNF, border = "red", add = TRUE) # plots Trinity/Bonavista with border

# Gulf of St. Lawrence/Scotian Shelf Filter
plot(water, xlim = c(-70, -62), ylim = c(40, 52), col = "lightblue", border = NA) # plot Gulf of St. Lawrense
x <- c(-55.56948, -55.44441, -73.62120, -73.12092, -66.49223, -56.69511, -56.57004, -57.07032, -54.56892, -55.61117) # captured using locator()
y <- c(46.91320, 42.37340, 42.37340, 47.35539, 51.06977, 51.42352, 50.39175, 49.21258, 47.82706, 46.97216)
p3 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
proj4string(p3) <- proj4string(water)
gulf <- gIntersection(water, p3, byid = TRUE) # intesection of water (ocean surronding NL) and Gulf of St. Lawrence
gulf <- gBuffer(gulf, width = 0.1) # add a small buffer
plot(gulf, border = "red", add = TRUE) # plots gulf with border

# plot North Atlantic (Gulf of Maine to Labrador Sea)
plot(water, col = "lightblue", border = NA) 
west <- extent(water)           # this makes the rectangle
west@ymin <- 41                 # this lowers the rectangle to 55 lat
west@ymax <- 52                 # this lowers the rectangle to 55 lat
west@xmin <- -72                # this lowers the rectangle to 55 lat
west@xmax <- -56                # this lowers the rectangle to 55 lat
west <- as(west, "SpatialPolygons")   # makes west a SpatialPolygon
proj4string(west) <- proj4string(water)
plot(west, border = "red", add = TRUE) # plot red rec


## Ignoring Lake Melville made little difference
filters <- gUnaryUnion(rbind(melville, north55))
filters <- gUnaryUnion(rbind(filters, eastNF))
#filters <- gUnaryUnion(rbind(filters, gulf))
filters <- gUnaryUnion(rbind(filters, gulf))
filters <- gUnaryUnion(rbind(filters, west, makeUniqueIDs = T)) # there were nonunique IDs in this line - the addition of makeUniqueIDs seems to work.
slot(filters@polygons[[1]], "ID") <- "- filters"
plot(water, col = "lightblue", border = NA)
plot(filters, border = "red", add = TRUE) # plot main map with filter


## Area and volume calculations ------------------------------------------------
# calculate the area of the sea ice and .....

converted <- list.files("sp_data/", pattern = ".Rdata")
load("ice_trends.Rdata")  # this is circular but is working under the assumption that code has been run in past years and that there is already and ice_trends.Rdata file.  So just load it and continue on with the code - it should only be a year out of date.
trends.calculated <- trends$date # this is from ice_trends.Rdata
dates3 <- strptime(converted, "%Y%m%d.Rdata") # use this for the full run
dates3 <- dates3[!format(dates3, "%Y%m%d") %in% format(trends.calculated, "%Y%m%d")] # these are dates that are not in ice_trends.Rdata- hence it is small


## subset the ice data
# create "model sets"

ct <- c(1:10)
sa <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "1.", "4.", "7.", "8.", "9.") # don't have "ice of land origin" or "undetermined or unknown"
m1 <- list(ct=ct, sa=sa)
m1

ct <- c(3:10)
sa <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
m2 <- list(ct=ct, sa=sa)
m2

ct <- c(8, 9)
sa <- c("5")
m3 <- list(ct=ct, sa=sa)
m3


# use function to calc Area/Volumes of ice
# works OK
# this is just to test on small batches
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v2a.R")

trends_update_subset <- calcAreaVolume(dates3[-c(23)], ct=m1$ct, sa=m1$sa)
trends <- rbind(trends, na.omit(data.frame(date = dates3[1:10], 
                                  area = trends_update_subset$areas, 
                                  volume = trends_update_subset$volumes,
                                  minlats = trends_update_subset$minlats,
                                  minlongs = trends_update_subset$minlongs
                                  )))  # it makes no sense to bind minlat to this because the original ice-trends.Rdata does not have this value - therefore, updating it makes no sense.

trends_update_subset <- calcAreaVolume(dates3[1:10], ct=m1$ct, sa=m1$sa)
trends <- rbind(data.frame(date = dates3[1:10], 
                             area = trends_update_subset$areas, 
                             volume = trends_update_subset$volumes,
                             minlats = trends_update_subset$minlats,
                             minlongs = trends_update_subset$minlongs))  # it makes no sense to bind minlat to th
save(trends, file = "ice_trends_subset-test.Rdata")

# test with NULL
# works OK but seems to break if run more than once
trends_update2 <- calcAreaVolume(dates3, ct=m1$ct, sa=m1$sa)
trends.m1 <- rbind(data.frame(date = dates3, 
                           area = trends_update2$areas, 
                           volume = trends_update2$volumes,
                           minlats = trends_update2$minlats,
                           minlongs = trends_update2$minlongs)) 

# this gets errors if you subset on dates3 as in some of the tests here.  Works fine if calcAreaVolume done on full dates3 vector
lookAt(trends.m1)
trends.m1 <- trends.m1[order(trends.m1$date), ]
save(trends.m1, file = "ice-trends-2017-m1-all.Rdata")

#source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v2.R")

# test with all of date3
# problem with 2002-07-15
trends_update3 <- calcAreaVolume(dates3, ct=m2$ct, sa=m2$sa)
trends.m2 <- rbind(data.frame(date = dates3, 
                           area = trends_update3$areas, 
                           volume = trends_update3$volumes, 
                           minlats = trends_update3$minlats,
                           minlongs = trends_update3$minlongs )) 
trends.m2 <- trends.m2[order(trends.m2$date), ]
save(trends.m2, file = "ice-trends-2017-m2-all.Rdata")


# test with other values
# works OK
trends_update3 <- calcAreaVolume(dates3, ct=m3$ct, sa=m3$sa)




str(dates3)
## Annual maps -----------------------------------------------------------------

# yr <- 2003
# converted <- list.files("sp_data/", pattern = ".Rdata")
# dates <- strptime(converted, "%Y%m%d.Rdata")
# dates <- dates[format(dates, "%Y") == as.character(yr)]
# for(i in seq_along(dates)) {
#   load(format(dates[i], "sp_data/%Y%m%d.Rdata"))
#   plotIce(ice, main = dates[i])
# }
# not used

## Manual discard --------------------------------------------------------------
dates3[-c(23)]

trends$area[format(trends$date, "%Y%m%d") == "19910401"] <- NA # polygons were miss-classified
trends <- trends[trends$date > ISOdate(1969, 11, 1), ] # 1969 was a strange year - difficult to assess patterns/fit curve - unreliable data?
trends <- na.omit(trends)
save(trends, file = "ice_trends_2017.Rdata")
