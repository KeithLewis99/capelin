
################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 2017-09-29, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-09-29v #
################################################################

# The purpose of this file is to:
#1) Make filter polygons

# Main issues:
# none for now

rm(list=ls()) # Clear the workspace

## Set-up ----------------------------------------------------------------------

# devtools::install_github("eblondel/cleangeo") # for cleaning SpatialPolygon objects; also on CRAN
options(stringsAsFactors = FALSE)

#setwd("C:/Users/Paul/Documents/DFO/ice")
setwd("D:/Keith/capelin/2017-project")

# this may be overkill on libraries but not sure what does what
library(RArcInfo)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(RColorBrewer)
library(cleangeo)

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

# Newfoundland North Coast and Bays
plot(water, xlim = c(-53.45, -53.35), ylim = c(47.8, 50.5), col = "lightblue", border = NA)
x <- c(-54.00649, -54.53953, -55.27395, -55.57009, -55.35687, -56.16237, -56.30451, -55.57009, -54.29078, -54.06572)
y <- c(49.44365, 49.18022, 49.04850, 49.07949, 49.30418, 49.35842, 49.61410, 49.95501, 49.59860, 49.44365)
p4 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
proj4string(p4) <- proj4string(water)
northNF <- gIntersection(water, p4, byid = TRUE) # intesection of water (ocean surronding NL) and Lake Melville
northNF <- gBuffer(northNF, width = 0.1) # add a small buffer
plot(northNF, border = "red", add = TRUE) # plots Trinity/Bonavista with border


# Cape Freels to Musgrave Harbour
plot(water, xlim = c(-53.45, -53.35), ylim = c(47.8, 50.5), col = "lightblue", border = NA)
x <- c(-54.12189, -53.10497, -53.46905, -53.90846, -54.48597, -54.57385, -54.13444)
y <- c(49.74124, 49.27318, 48.92007, 48.92007, 49.15821, 49.47026, 49.70018)
p41 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
proj4string(p41) <- proj4string(water)
capeFreels <- gIntersection(water, p41, byid = TRUE) # intesection of water (ocean surronding NL) and Lake Melville
capeFreels <- gBuffer(capeFreels, width = 0.1) # add a small buffer
plot(capeFreels, border = "red", add = TRUE) # plots Trinity/Bonavista with border

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
filters <- gUnaryUnion(rbind(filters, gulf))
filters <- gUnaryUnion(rbind(filters, gulf))
filters <- gUnaryUnion(rbind(filters, northNF))
filters <- gUnaryUnion(rbind(filters, capeFreels))
filters <- gUnaryUnion(rbind(filters, west, makeUniqueIDs = T)) # there were nonunique IDs in this line - the addition of makeUniqueIDs seems to work.
slot(filters@polygons[[1]], "ID") <- "- filters"
plot(water, col = "lightblue", border = NA)
plot(filters, border = "red", add = TRUE) # plot main map with filter

save(filters, file = "output-processing/filters.Rdata")

