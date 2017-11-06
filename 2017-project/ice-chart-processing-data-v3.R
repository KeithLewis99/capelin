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
  #6) Area and volume calculations

# Main issues (see README):
  #1) Few errors in downloading files from the internet or converting e00 to SpatialPolygonDataframe
  #2) Too many connections open when converting files
  #3) Lots of long loops that could be made into functions to simplify code. Done.
  #4) can date computations be simplified with lubridate? yes but why change code now?
  #5) circular problem with ice_trends data. See explanation below
  #6) how to get into archive so that we can see what has been missed and why. Can get into archive.
  #7) make sure that we can extract the same values as in Ale's paper to ensure this is accurately working (maxIce, tice) - done

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


#setwd("C:/Users/Paul/Documents/DFO/ice")
#setwd("D:/Keith/capelin/2017-project")

library(RArcInfo)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
#devtools::install_github("tidyverse/ggplot2")
library(RColorBrewer)
library(cleangeo)
#library(ggplot2)
#library(dplyr)
#library(tidyr)
library(broom)
library(tidyverse)
#library(doParallel)
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

# from here to "subset ice data" can be carved off into another file

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
dates <- dates[!as.Date(dates) %in% as.Date(downloaded)] #this line is for downloading all data from the EC website - takes A LONG TIME!!!!!
dates <- format(dates, "%Y%m%d")
dates

# downloads all files not already in e00_data folder
# error messages indicate maps that have not been or may never be created.  Can be ignored. See previous comment
e00Download(dates)

## remove problem files at this stage
file.remove("e00_data/20070508.e00")

## e00 file conversion to Spatial Polygons Dataframe-------------------
# e00 to avc_data (coverages) and then to SpatialPolygonsDAtaframe in sp_data

## note: had to run this across multiple sessions...too many open files.
##       don't know how to close the connections open by RArcInfo..grrrrr


## convert e00 files to RData--------------
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
# see "code-map-check.xlsx: unconverted maps"

# refresh converted and calculate percentage of files that did not convert to a SPDF
converted <- list.files("sp_data/", pattern = ".Rdata")
converted <- strptime(converted, "%Y%m%d.Rdata")
dates2 <- downloaded[!downloaded %in% converted] # files left to convert
length(dates2)/length(converted) # conversion didn't work out for a small number of e00 files  - See "code-map-check.xlsx:unconverted maps" for details


## remove problem files at this stage
# note that there is no corresponding e00 file on the EC website.  PR suspects that it may be hidden.  More imortantly, there is one from 20070507 (one day before).  Eliminate file
file.remove("sp_data/20070508.Rdata")
unlink("avc_data/20070508") # not sure why we removed this folder but it is a fragment and not needed

# remove these dates before running due to fragments, i.e., ice fragments that came very low or where there was some question what tice was.  See "code-map-check.xlsx:area-tice" for details
#"1991-02-18 NST", "1992-03-19 NST", "2001-04-09 NDT", "2007-04-23 NDT", "2011-04-04", "2013-04-15", "2015-04-27 NDT"
file.remove("sp_data/19910218.Rdata")
file.remove("sp_data/19920319.Rdata")
file.remove("sp_data/20010409.Rdata")
file.remove("sp_data/20070423.Rdata")
file.remove("sp_data/20110404.Rdata")
#file.remove("sp_data/20130415.Rdata") tried removing dates but this is such a wonky year that it doesn't really work (too many chunks of ice in the bays or just outside of the filters).  Went over the 
file.remove("sp_data/20150427.Rdata")

# second round of removals
file.remove("sp_data/19920402.Rdata")
file.remove("sp_data/20070507.Rdata")
#file.remove("sp_data/20130422.Rdata")
file.remove("sp_data/20150323.Rdata")
file.remove("sp_data/20110214.Rdata")

# third round of removals
#file.remove("sp_data/20130128.Rdata")
## Area and volume calculations ------------------------------------------------
# calculate the area of the sea ice and .....

converted <- list.files("sp_data/", pattern = ".Rdata")
load("ice_trends.Rdata")  # this is circular but is working under the assumption that code has been run in past years and that there is already and ice_trends.Rdata file.  So just load it and continue on with the code - it should only be a year out of date.
trends.calculated <- trends$date # this is from ice_trends.Rdata
dates3 <- strptime(converted, "%Y%m%d.Rdata") # use this for the full run
save(dates3, file = "output-processing/dates3all.Rdata")
dates3 <- dates3[!format(dates3, "%Y%m%d") %in% format(trends.calculated, "%Y%m%d")] # these are dates that are not in ice_trends.Rdata- hence it is small
save(dates3, file = "output-processing/dates3.Rdata")

##########################################################################
## subset the ice data
#
# create "model sets"
#load("output-processing/dates3.Rdata")
load("output-processing/dates3all.Rdata")
load("output-processing/subset-lists.Rdata")
load("output-processing/filters.Rdata")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

# calculate proper dates for 2013 and manually insert below
trends_2013 <- calcAreaVolLat(dates3[1306:1330], ct=m1$ct, sa=m1$sa, sb=m1$sb)
trends.2013 <- iceOutput(trends_2013, dates3[1306:1330])
max(trends.2013$area) 
maxarea2013 <- trends.2013[5,] # this is the proper value :: keep
min(trends.2013$minlats) # the minlats is actually 2013-04-15 but bc this was a fragment
minlats2013 <- trends.2013[8,]
minlats2013 <- trends.2013[-15, -16] #remove the fragments
lookAt(trends.m1)
#View(trends.m1)
save(trends.m1, file = "output-processing/ice-trends-2013-m1.Rdata")

# it makes no sense to bind minlat to this because the original ice-trends.Rdata does not have this value - therefore, updating it makes no sense.

# this gets errors if you subset on dates3 as in some of the tests here.  Works fine if calcAreaVolume done on full dates3 vector

## m1
trends_update.m1 <- calcAreaVolLat(dates3, ct=m1$ct, sa=m1$sa, sb=m1$sb)
trends.m1 <- iceOutput(trends_update.m1, dates3)
lookAt(trends.m1)
#View(trends.m1)
save(trends.m1, file = "output-processing/ice-trends-2017-m1-alla.Rdata")
#save(trends.m1, file = "output-processing/ice-trends-2017-m1-subset.Rdata")

# check to insure that subset is less than full set
x <- trends_update.m1$areas - trends_update3$areas
subset(x, x < 0)

## m2
trends_update.m2 <- calcAreaVolLat(dates3, ct=m2$ct, sa=m2$sa,  sb=m2$sb)
trends.m2 <- iceOutput(trends_update.m2, dates3)
lookAt(trends.m2)
save(trends.m2, file = "output-processing/ice-trends-2017-m2-all.Rdata")
#save(trends.m2, file = "output-processing/ice-trends-2017-m2-subset.Rdata")

## m3
trends_update.m3 <- calcAreaVolLat(dates3, ct=m3$ct, sa=m3$sa,  sb=m3$sb)
trends.m3 <- iceOutput(trends_update.m3, dates3)
lookAt(trends.m3)
save(trends.m3, file = "output-processing/ice-trends-2017-m3-all.Rdata")
#save(trends.m3, file = "output-processing/ice-trends-2017-m3-subset.Rdata")

## m4
trends_update.m4 <- calcAreaVolLat(dates3, ct=m4$ct, sa=m4$sa,  sb=m4$sb)
trends.m4 <- iceOutput(trends_update.m4, dates3)
lookAt(trends.m4)
#save(trends.m4, file = "output-processing/ice-trends-2017-m4-subset.Rdata")
save(trends.m4, file = "output-processing/ice-trends-2017-m4-all.Rdata")

## m5
trends_update.m5 <- calcAreaVolLat(dates3, ct=m5$ct, sa=m5$sa,  sb=m5$sb)
trends.m5 <- iceOutput(trends_update.m5, dates3)
lookAt(trends.m5)
#save(trends.m4, file = "output-processing/ice-trends-2017-m4-subset.Rdata")
save(trends.m5, file = "output-processing/ice-trends-2017-m5-all.Rdata")

## m6
trends_update.m6 <- calcAreaVolLat(dates3, ct=m6$ct, sa=m6$sa,  sb=m6$sb)
trends.m6 <- iceOutput(trends_update.m6, dates3)
lookAt(trends.m6)
#save(trends.m4, file = "output-processing/ice-trends-2017-m4-subset.Rdata")
save(trends.m6, file = "output-processing/ice-trends-2017-m6-all.Rdata")

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

dates3[-c(618, 656, 933, 1125, 1257, 1328, 1404)]

dates3[1328]
trends$area[format(trends$date, "%Y%m%d") == "19910401"] <- NA # polygons were miss-classified
trends <- trends[trends$date > ISOdate(1969, 11, 1), ] # 1969 was a strange year - difficult to assess patterns/fit curve - unreliable data?
trends <- na.omit(trends)
save(trends, file = "\output-processing\ice_trends_2017.Rdata")
