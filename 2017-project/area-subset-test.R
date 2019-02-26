rm(list=ls())

test1 <- lookAtIce(dates3, 25, ct=m1$ct, sa=m1$sa, sb=m1$sb)

test1$a
egg <- subsetProject(test1$a)    # ice is not in the local environment
sub.egg <- filterEgg(egg)

sub.egg.lcc <- sp::spTransform(sub.egg, CRS(proj4string(test1$a)))

area <- sapply(slot(sub.egg, "polygons"), function(x) slot(x, "area"))
area
sub.egg@data$AREA/1000^2

area.lcc <- sapply(slot(sub.egg.lcc, "polygons"), function(x) slot(x, "area"))
area.lcc/1000^2
sub.egg.lcc@data$AREA

a <- try(gArea(sub.egg.lcc, byid = TRUE))
sub.egg.lcc@data$AREA <- a * 1e-6 # replace polygon area (use square km)
sub.egg.lcc@data$AREAice <- sub.egg.lcc@data$AREA*as.numeric(sub.egg@data$N_CT)/10
sub.egg.lcc@data$AREA

#plot
zz <- SpatialPolygons((sub.egg@polygons))
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)
plot(water, xlim = c(-60, -53), ylim = c(50, 56), col = "lightblue", border = NA) 
plot(water, xlim = c(-54.5, -52.5), ylim = c(51, 51.5), col = "lightblue", border = NA) 

plot(sub.egg, add = T)
plot(zz[1], border = "cyan", add = T)  # C filter?
plot(zz[2], border = "chartreuse", add = T, lwd=4) #
plot(zz[3], border = "black", add = T, lwd=4)
plot(zz[4], border = "red", add = T, lwd=4)
plot(zz[5], border = "red", add = T, lwd=4)
plot(zz[6], border = "red", add = T, lwd=4)
plot(zz[7], border = "red", add = T, lwd=4)
plot(zz[8], border = "red", add = T, lwd=4)
plot(zz[9], border = "red", add = T, lwd=4)
plot(zz[10], border = "red", add = T, lwd=4)
plot(zz[11], border = "red", add = T, lwd=4)
plot(zz[12], border = "red", add = T, lwd=4)
plot(zz[13], border = "red", add = T, lwd=4)
plot(zz[14], border = "cyan", add = T, lwd=4)
plot(zz[15], border = "green", add = T, lwd=4)
plot(zz[16], border = "orange", add = T, lwd=4)

#subset
lattest <- sub.egg.lcc[sub.egg.lcc@data$AREAice > 100, ]
sum(lattest@data$AREAice)/sum(sub.egg.lcc@data$AREAice)

#plot subsetted polygons
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)
plot(water, xlim = c(-60, -53), ylim = c(50, 56), col = "lightblue", border = NA) 

yy <- SpatialPolygons((sub.egg.lcc@polygons))
x <- spTransform(lattest, CRS("+proj=longlat +datum=WGS84"))
plot(x, add = T)


# trying to get the polygons to form a union is not straight forward and its not clear to me how to subset the most southerly one.
test <- unionSpatialPolygons(sub.egg, 16)
#browser()
#print("inside IT")
#print(environment())
#print(ls())
d <- tidy(sub.egg)
tid <- d[which.min(d$lat),]
return(tid)


# try this with a union of polygons
nc1 <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],proj4string=CRS("+proj=longlat +datum=NAD27"))
lps <- coordinates(nc1)
ID <- cut(lps[,1], quantile(lps[,1]), include.lowest=TRUE)
reg4 <- unionSpatialPolygons(nc1, ID)
row.names(reg4)


lps <- coordinates(sub.egg)
ID <- cut(lps[,1], quantile(lps[,1]), include.lowest=TRUE)
test <- unionSpatialPolygons(sub.egg, ID)



plot(test)
lps <- coordinates(test)
ID <- cut(lps[,1], quantile(lps[,1]), include.lowest=TRUE)
test1 <- unionSpatialPolygons(test, ID)
plot(test1)

plot(test)
plot(test[1], border = "cyan", add = T)  # C filter?
plot(test[2], border = "chartreuse", add = T, lwd=4) #
plot(test[3], border = "black", add = T, lwd=4)
plot(test[4], border = "red", add = T, lwd=4)
plot(sub.egg)


test <- gUnion(sub.egg, id = row.names(sub.egg@data))
plot(test)

row.names(sub.egg@data) <- row.names(sub.egg@data)
sub.egg <- spChFIDs(sub.egg, row.names(sub.egg@data))
test <- gUnaryUnion(sub.egg, id = row.names(sub.egg@data))
plot(test)


setwd("D:/Keith/R/spatial-test/states")
#
# script demonstrares polygon dissolve using maptools
# and area calculation using PBSmapping package.
#
# first, load the R packages that we will need
#
library(maptools)   # for geospatial services; also loads foreign and sp
library(gpclib)     # General Polygon Clipping library 
library(rgdal)      # for map projection work; also loads sp
library(PBSmapping) # for GIS_like geospatial object manipulation / anslysis including poly
# 
# Part 1:
# First example, derived from the maptools unionSpatialPolygon() method documentation in the R 
# maptools user manual: We calculate and sum the area of each polygon in the North Carolina map. 
# Then, use unionSpatialPolygons() to dissolve the county polygon boundries, and assign each
# polygon's area to one of four regions, based on longitude thresholds.
#
print("start of demo. Hit key...")
browser()
print("reading and transforming North Carolina Shapefile...")
#
NorthCaroBase <- readShapePoly("sids.shp",proj4string=CRS("+proj=aea +ellps=GRS80 +datum=WGS84"))
#
# Transform the polygons (which were read in as unprojected geographic coordinates)
# to an Albers Equal Area projection.
#
NorthCaroProj = spTransform(NorthCaroBase,CRS("+proj=aea +ellps=GRS80 +datum=WGS84"))
#
# Convert to a PolygonSet for compatability with PBSmapping package routines.
#
NorthCaroProjPS = SpatialPolygons2PolySet(NorthCaroProj)
#
plotPolys(NorthCaroProjPS, proj = TRUE,col="wheat1",xlab="longitude",ylab="latitude")
#
# polygon area calculations
#
print("Calculating North Carolina polygon areas...")
attr(NorthCaroProjPS, "projection") <- "LL"
NCPolyAreas = calcArea(NorthCaroProjPS,rollup=1)
#
# Compute the area of state by summing the area of the counties.
#
numCountyPolys = length(NCPolyAreas[,1])
NCArea = sum(NCPolyAreas[1:numCountyPolys,2])
print(sprintf("North Carolina (%d county polygons) area: %g sq km. Hit key to continue",
              numCountyPolys,NCArea))
browser()
#
# create a set of 'breakpoints' that unionSpatialPolygons() method will use to 
# place the dissolved polygons into one of four longitude 'zones' on the output map: 
# Each county/polygon's x/y label coordinates gives the location of that polygon's center.
# 
print("Dissolving North Carolina polygons...")

lps <- getSpPPolygonsLabptSlots(NorthCaroProj)
# 
# Assign each county to one of four longitudinal 'bins': 
#
IDFourBins <- cut(lps[,1], quantile(lps[,1]), include.lowest=TRUE)
#
# Dissolve operations: result is a SpatialPolygons object.
# Convert to PolySet to display new polygons and calculate areas.0
#
NcDissolve   <- unionSpatialPolygons(NorthCaroProj ,IDFourBins)
NcDissolvePSFour <- SpatialPolygons2PolySet(NcDissolve)
plotPolys(NcDissolvePSFour, proj = TRUE,col="wheat1",xlab="longitude",ylab="latitude")
#
# projecton attribute must be "UTM" or "LL" for areas in km ** 2
#
attr(NcDissolvePSFour, "projection") <- "LL"
NCDissolvePolyAreas = calcArea(NcDissolvePSFour ,rollup=1)
#
# sum the areas of all polygons in the poly set - 
# the expression 'length(NCDissolvePolyAreas [,1]' returns the number of polygons
#
NCDissolveArea = sum(NCDissolvePolyAreas [1:length(NCDissolvePolyAreas [,1]),2])
#
print(sprintf("North Carolina (four regions) area: %g sq km. Hit key to continue",
              NCDissolveArea))
browser()
#
# next, lets put all county polygons into a single 'bin' to get just a county outline: 
# replace the quantile() call with one to range(), that returns just the min/max of the lps[,1] vector.
#
print("Dissolving North Carolina polygons (again)...")
IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
NcDissolve   <- unionSpatialPolygons(NorthCaroProj ,IDOneBin)
NcDissolvePSOne <- SpatialPolygons2PolySet(NcDissolve)

plotPolys(NcDissolvePSOne, proj = TRUE,col="wheat1",xlab="longitude",ylab="latitude")
#
# projecton attribute must be "UTM" or "LL" for areas in km ** 2
#
print("Calculating North Carolina polygon area (again)...")
attr(NcDissolvePSOne, "projection") <- "LL"
NCDissolvePolyAreas = calcArea(NcDissolvePSOne ,rollup=1)
print(sprintf("North Carolina (1 region) area: %g sq km. Hit key to continue",NCDissolveArea))
#
browser()
#
# close all the graphics devices
#