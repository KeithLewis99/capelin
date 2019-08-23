################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2019-05-29 #
#  Project taken over by Aaron Adamack in 2019; Lewis helped modify scripts from original 2017 project
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
  #8) Major modification from 2017 project to clean up unneeded code and streamline the project for future stock assessments

# all of the downloading below could probably been done more efficiently with an API!


# rest to output file
# Note that the ice files (the 'egg' data) come with different file names (below are several examples)
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


library(RArcInfo)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(RColorBrewer)
library(cleangeo)
library(broom)
library(tidyverse)

## read in source code (Lewis functions)-----
source("ice-chart-processing-function-v3.R")
load("output-processing/subset-lists.Rdata") # this is now only needed to bring in "m1" which is all of the ice egg data.  There were subsets m2-m6 which removed ice stages and types.
# this file removes ice north of 55deg lat, Lake Melville, the Bays around NF, the Gulf of St. Lawrence, and the Scotian Shelf
load("output-processing/filters.Rdata")

## create dirs for data storage
if(!dir.exists("avc_data")) dir.create("avc_data")
if(!dir.exists("e00_data")) dir.create("e00_data")
if(!dir.exists("sp_data")) dir.create("sp_data")
if(!dir.exists("output-processing")) dir.create("output-processing")
if(!dir.exists("figs")) dir.create("figs")
if(!dir.exists("data")) dir.create("data")

## Metadata --------------------------------------------------------------------
## For stage codes, see: (https://ec.gc.ca/glaces-ice/default.asp?lang=En&n=D5F7EA14-1&offset=1&toc=show) or https://www.canada.ca/en/environment-climate-change/services/ice-forecasts-observations/latest-conditions/archive-overview/information-about-data.html
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
#get rid of centimeters in stages
stages$thickness <- gsub(" centimetres", "", stages$thickness) 

# create a list of thickness values in stages
trange <- strsplit(stages$thickness, " - ") 

# creates a single numeric value based on thickness (mean for ranges, same value for <>)
stages$thickness <- sapply(seq_along(trange), function(i) {
  if(length(trange[[i]]) == 2) {
    median(as.numeric(trange[[i]]))
  } else {
    as.numeric(gsub("[^0-9]", "", trange[[i]]))
  }
}) 


## Download e00 files ----------------------------------------------------------
# this is very robust.  Seems to be only 1 error in almost 1500 files - however, if e00 files are deleted, the code does not always seem to look backwards.  Use with caution and check if there is reason to believe an e00 file has been deleted.

# this code assumes at least one e00 file in directory, i.e., seed the directory by downloading a file from https://iceweb1.cis.ec.gc.ca/Archive/page1.xhtml (choose e00 data in black and white for the NL region)
downloaded <- list.files("e00_data/", pattern = ".e00") #lists the e00 files
downloaded <- downloaded[downloaded != "test.e00"] # removes the test.e00
downloaded <- strptime(downloaded, "%Y%m%d.e00") # takes the date from the e00 files

###############################################################################
##### HEADS UP, FIX THE DATE TO CURRENT DATE!!!!!!!!!---- ##############################
####################################################################################
#this creates a POSIXct vector of all possible dates between the dates specified i.e 50 years!!
dates <- seq(ISOdate(1969, 1, 17), ISOdate(2019, 6, 19), by = "day") # update the second ISOdate to present as desired - however, this does not work well if not updated: also, it reports a numer of erros, one for each day since last e00 file - these should be checked but can generally be ignored.

#this makes dates just for any months that could have ice
dates <- dates[format(dates, "%m") %in% c("11", "12", "01", "02", "03", "04", "05", "06", "07")]

####################################################################
#############STOP
#################################################################### See N. Soontiens - 
#run one of next two lines
#this line is for simply updating - assumes all past files are downloaded and only looks for new files on EC website past the terminal data supplied in line 103 above; basically just looks for dates not downloaded by 1) finding the max date in downloaded, 2) finding the difference between "downloaded" and "dates" and keeping values < -1 ie only the most recent dates that would need to downloaded
dates <- dates[difftime(max(downloaded), dates, units = "days") < -1] 
#this line is for downloading all data from the EC website - takes A LONG TIME!!!!! #dates <- dates[!as.Date(dates) %in% as.Date(downloaded)] 
dates <- format(dates, "%Y%m%d")

# downloads all files not already in e00_data folder
# error messages indicate maps that have not been or may never be created.  Can be ignored. See previous comment.
# Note that if doing from 1969 to present, this will take a long time. Most will say error but that is just because this function looks at ~4 different naming conventions and all possible dates.  So expect to see "Error in download.file" but don't worry - its doing what it is supposed to.
e00Download(dates)


## e00 file conversion to Spatial Polygons Dataframe-------------------
# e00 to avc_data (coverages) and then to SpatialPolygonsDAtaframe in sp_data

## note: originally, I had to run this across multiple sessions...too many open files.
##       don't know how to close the connections open by RArcInfo..grrrrr


## convert e00 files to RData--------------
# this is just to refresh the download file to ensure that all dates accounted for - see lines 105-107 for an explanation
downloaded <- list.files("e00_data/", pattern = ".e00") # this code assumes at least one e00 file in directory, i.e., seed the directory
downloaded <- downloaded[downloaded != "test.e00"]
downloaded <- strptime(downloaded, "%Y%m%d.e00")
# this does as in lines 105-107 for the .Rdata  files
converted <- list.files("sp_data/", pattern = ".Rdata")
converted <- strptime(converted, "%Y%m%d.Rdata")
dates1 <- downloaded[!downloaded %in% converted] # files left to convert ###########not sure this is working right # KL to check these dates against GIF files!!! Note that length dates1 = 26 but dates2 =16
dates1 <- format(dates1, "%Y%m%d")
dates1

#"19910218" "19920319" "19920402" "20010409" "20070423" "20070507" "20110214" "20110404" 20150323" "20150427" "20160208"

# helper function to minimize the code
# e00 to avc_data (coverages) and then to SpatialPolygonsDAtaframe in sp_data
e00_to_SpatialPolygonDataframe(dates1)
# see "code-map-check.xlsx: unconverted maps"

# refresh converted and calculate percentage of files that did not convert to a SPDF - not neccessary for analysis but indicates if there is a problem with the code
converted <- list.files("sp_data/", pattern = ".Rdata")
converted <- strptime(converted, "%Y%m%d.Rdata")
dates2 <- downloaded[!downloaded %in% converted] # files left to convert
length(dates2)/length(converted) # conversion didn't work out for a small number of e00 files (usually about 1%) - See "code-map-check.xlsx:unconverted maps" for details

# [1] "1978-01-22 NST" "1984-03-13 NST" "1990-02-26 NST"  "2016-02-08 NST" - checked all of these dates; some were close but none were tice; 
# "2015-02-16 NST" there is no Rdata file or GIF for this date
#"1993-12-13 NST" "2011-07-18 NDT"[6] "2011-11-14 NST" "2011-11-21 NST" "2012-07-23 NDT" "2012-07-30 NDT" "2012-11-26 NST"[11] "2013-12-23 NST"  "2015-11-09 NST"  "2016-12-05 NST"[16] "2017-07-31 NDT" - none of these dates can possible influence analysis so did not check

## remove problem files at this stage
# note that there is no corresponding e00 file on the EC website.  PR suspects that it may be hidden.  More imortantly, there is one from 20070507 (one day before).  Eliminate file
file.remove("sp_data/20070508.Rdata")
unlink("avc_data/20070508") # not sure why we removed this folder but it is a fragment and not needed

# remove these dates before running due to fragments, i.e., ice fragments that came very low or where there was some question what tice was.  See "code-map-check.xlsx:area-tice" in conjunction with ice maps for details
#"1991-02-18 NST", "1992-03-19 NST", "2001-04-09 NDT", "2007-04-23 NDT", "2011-04-04", "2013-04-15", "2015-04-27 NDT"
file.remove("sp_data/19910218.Rdata")
file.remove("sp_data/19920319.Rdata")
file.remove("sp_data/20010409.Rdata")
file.remove("sp_data/20070423.Rdata")
file.remove("sp_data/20110404.Rdata")
#file.remove("sp_data/20130415.Rdata") tried removing dates but this is such a wonky year that it doesn't really work (too many chunks of ice in the bays or just outside of the filters).  So did it manually - see line 205 below 
file.remove("sp_data/20150427.Rdata")

# second round of removals - as above: see "code-map-check.xlsx: area-tice" in conjuction with icemaps for details.
file.remove("sp_data/19920402.Rdata")
file.remove("sp_data/20070507.Rdata")
#file.remove("sp_data/20130422.Rdata") - see explantion above for hashing out 2013
file.remove("sp_data/20150323.Rdata")
file.remove("sp_data/20110214.Rdata")

# third round of removals
#file.remove("sp_data/20130128.Rdata") - see explantion above for hashing out 2013


## Area and volume calculations ------------------------------------------------
# calculate the area of the sea ice and .....
#list the .Rdata files and put in object "converted"
converted <- list.files("sp_data/", pattern = ".Rdata")
dates3 <- strptime(converted, "%Y%m%d.Rdata") # use this for the full run

##########################################################################
# leave the loads in ##  because sometimes its good to have these when troubleshooting and you don't want to rerun code

# calculate proper dates for 2013 and manually in area-ice.R
# note that 2013 was a very wonky ice year with 3-4 retreats and a large ice fragment on Apr 22 at < 49 deg latitude.  Because minlat was not a good proxy for tice in this year, I looked at the maps (D:\Keith\ice\Ice_Charts\Original_Ice_Data\GIF) and determined tice visually as Ale did in teh original work (Buren et al. 2014).  Note that it was a close call between February and April.  Basically, we decided that some low latitude ice hugging the coast in April but with little ice to the north in not the same as a lot of ice above the minlat in Feb and that the latter was more in accord with the spirit of Ale's work.
trends_2013 <- calcAreaVolLat(dates3[1306:1330], ct=m1$ct, sa=m1$sa, sb=m1$sb)
trends.2013 <- iceOutput(trends_2013, dates3[1306:1330]) #makes trends_2013 a dataframe
max(trends.2013$area) 
maxarea2013 <- trends.2013[5,] # this is the proper value based on a visual assessment of the ice maps :: keep
min(trends.2013$minlats) # the minlats is actually 2013-04-15 but this was bc of a fragment hugging the coast
minlats2013 <- trends.2013[8,] #based on the description above, I chose this day for tice, 2013-02-25 and made teh corrections in area-ice.R on line 71-76
minlats2013 <- trends.2013[c(-15, -16), ] #remove the dates with ice fragments
lookAt(trends.m1)
lookAt(trends.2013)
#View(trends.m1)
save(trends.2013, file = "output-processing/ice-trends-2013-m1.Rdata") #this Rdata is not really used but again, saves having to run the above code

## m1
trends_update.m1 <- calcAreaVolLat(dates3, ct=m1$ct, sa=m1$sa, sb=m1$sb)
trends.m1 <- iceOutput(trends_update.m1, dates3)
lookAt(trends.m1)
#View(trends.m1)
save(trends.m1, file = "output-processing/ice-trends-2019-m1-alla.Rdata")
#save(trends.m1, file = "output-processing/ice-trends-2017-m1-subset.Rdata")

## Annual maps -----------------------------------------------------------------
# this is for viewing the flow of ice for a given year - provides an animation.  Not needed for the analysis and is similar to going through the gif files
# yr <- 2003
# converted <- list.files("sp_data/", pattern = ".Rdata")
# dates <- strptime(converted, "%Y%m%d.Rdata")
# dates <- dates[format(dates, "%Y") == as.character(yr)]
# for(i in seq_along(dates)) {
#   load(format(dates[i], "sp_data/%Y%m%d.Rdata"))
#   plotIce(ice, main = dates[i])
# }


## Manual discard --------------------------------------------------------------
#Not sure what below if for but it should have no influence on the analysis
#dates3[-c(23)]
#23 = "1970-01-16 NST" - no egg data on the ice charts for this day


#dates3[-c(618, 656, 933, 1125, 1257, 1328, 1404)]
#618 = "1991-02-25 NST" - not tice; there is an .Rdata file - not sure why this was selected
#656 = 1992-04-09 NDT" - as above
#933 = "2001-05-07 NDT" - as above
#1125 = "2007-06-04 NDT" - as above
#1257="2011-05-30 NDT" - as above
#1328 = "2013-06-10 NDT" - as above
#1404 = "2015-07-06 NDT" - as above

# the below doesn't seem to do anything.  There is no date or Rdata but there is a GIF file for 19910401.  Therefore, there is no place for the NA to go and there appear to be no NAs to omit.  DEPRECATE
#load("ice_trends.Rdata")
#trends$area[format(trends$date, "%Y%m%d") == "19910401"] <- NA # polygons were miss-classified; there is no Rdata file for this date
#trends <- trends[trends$date > ISOdate(1969, 11, 1), ] # 1969 was a strange year - difficult to assess patterns/fit curve - unreliable data?
#trends <- na.omit(trends)
#save(trends, file = "output-processing/ice_trends_2019.Rdata")
