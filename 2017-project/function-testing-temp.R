# The question here is, what is the right area? calcAreaVolLat() faithfully reproduces an areas calculation but is it right?  I have tried to take this function apart to reproduce the values.  Unfortunately, its always different.  Why?  
# things seem to be changing from day to day - things that were the same yesterday aren't today.  Could be an environment thing??


rm(list=ls())
load("output-processing/dates3.Rdata")
load("output-processing/filters.Rdata")
load("output-processing/dates3all.Rdata")
load("output-processing/subset-lists.Rdata")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
options(stringsAsFactors = FALSE)
setwd("D:/Keith/capelin/2017-project")
options(scipen=999) 

dates3[1180]
rm(test)
rm(sub)
rm(test1)
rm(orig)

z <- dates3
i <- 1
z <- sub.egg1
x <- sub.egg1
y <- a

my.polygon@bbox <- as.matrix(extent(my.raster))
extent(sub.egg)
extent(sub.egg1)

## Compare the value of AREA and a = gArea------------
# extract essence of sub.egg1
test <- extractSubegg1(dates3, 25, ct=m1$ct, sa=m1$sa, sb=m1$sb)
test <- extractSubegg1(dates3, 25, ct=m2$ct, sa=m2$sa, sb=m2$sb)
test <- extractSubegg1(dates3, 25, ct=m3$ct, sa=m3$sa, sb=m3$sb)


test1 <- extractSPDFfinal(dates3, 1180, ct=m1$ct, sa=m1$sa, sb=m1$sb)
test2 <- extractSPDFfinal(dates3, 1180, ct=m2$ct, sa=m2$sa, sb=m2$sb)
test3 <- extractSPDFfinal(dates3, 1180, ct=m3$ct, sa=m3$sa, sb=m3$sb)

test[[1]]@data
sub <- test[[1]]@data[, c("AREA",  "AREAice", "A_LEGEND", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "SE","AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")] #
#extract all sub.egg - remember that ice is unfiltered and will have lots more polygons
head(sub, 20)
sapply(slot(test[[1]], "polygons"), function(x) slot(x, "ID"))

test1 <- lookAtIce(dates3, 592, ct=m1$ct, sa=m1$sa, sb=m1$sb)
test2 <- lookAtIce(dates3, 1, ct=m2$ct, sa=m2$sa, sb=m2$sb)
test3 <- lookAtIce(dates3, 1, ct=m2$ct, sa=m2$sa, sb=m2$sb)


orig <- test1$c@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "SE")]
head(orig, 20)
str(orig)
ice.orig <- test1$a@data[, c("AREA", "PERIMETER", "BIN#")]

# calculate total area from AREA for sub.egg1
AREAice <- sum(test[[1]]@data$AREAice)/1000000

# calculate total area from AREA of each ice type
A <- sum(test[[1]]@data$AREA_SA)/1000000
B <- sum(test[[1]]@data$AREA_SB)/1000000
C <- sum(test[[1]]@data$AREA_SC)/1000000
D <- sum(test[[1]]@data$AREA_SD)/1000000

TOT <- sum(A, B, C, D)

AREAice/TOT # this compares favourably to AREA 

# look at AREA and AREA[SX]
compare <- test[[1]]@data[, c("AREA", "AREAice", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")]
# sum the columns of AREA[SX]
sumcols <- (test[[1]]@data$AREA_SA + test[[1]]@data$AREA_SB + test[[1]]@data$AREA_SC + test[[1]]@data$AREA_SD)

# compare the values of AREA and sumcols - should be equal 
# note that this will change - sumcols should be a percentage of AREA corresponding to CT
as.data.frame(cbind(compare, sumcols))
as.data.frame(cbind(compare[,c("AREA", "AREAice")], sumcols))
# These match suggesting that the AREA_SX is summing correctly (internally consistent)
sum(sumcols)
sum(compare[, c("AREA")])
sum(compare[, c("AREAice")])


#why does sum of AREA not equal areas???? Externally inconsistent
#pick up where iceSubset leaves off
# this calculates area
sub.poly <- calcPolyA(test1$c) # calcPoly calculates the area of each polygon manually using either 1) gArea or 2) mannualy calculating the area from all the slot @area in each polygon

length(test1$c)
length(sub.poly)
sub.trend <- trendsCalc(test1$c, sub.poly) # converts the AREAS in the SPDF and sums all areas
test1$c
sub.poly

# this further breaks down the code to run it line by line
z <- test1$c

#z <- test[[1]]

  a <- try(gArea(z, byid = TRUE)) 
  str(a)
# from iceAREA  
  z@data$AREA <- a * 1e-6
  z@data$AREAice <- z@data$AREA*as.numeric(z@data$CT)/10
  subarea <- iceArea(z, ct=m1$ct, sa=m1$sa, sb=m1$sb)
  options(scipen=999)  

subarea[[1]][, c("AREA", "AREAice", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")]  
# sub.trend and subarea are equal    
  
  # calculation of areas does NOT match subarea (has been through all)
  #calcAreaVolLat(dates3[25], ct=m1$ct, sa=m1$sa, sb=m1$sb)  
  calcAreaVolLat(dates3[25], ct=m1$ct, sa=m1$sa, sb=m1$sb)  
  calcAreaVolLat(dates3[25], ct=m2$ct, sa=m2$sa, sb=m2$sb)  
  sub.poly1 <- calcPolyA(test[[1]]) 
  sub.trend1 <- trendsCalc(test[[1]], sub.poly1)

  # sum of AREA (has been through iceArea but not gArea)
  sum(test[[1]]@data$AREAice) # this value is teh same as AREA - good!!!!! 14 values
  sum(test[[1]]$AREA_SA + test[[1]]$AREA_SB + test[[1]]$AREA_SC + test[[1]]$AREA_SD)
  test[[2]]
  
  # sum of AREAS[SX] - has been through all 
  sum(subarea[[1]]$AREA_SA + subarea[[1]]$AREA_SB + subarea[[1]]$AREA_SC + subarea[[1]]$AREA_SD)
  sum(subarea[[1]]$AREAice)                  # 19 values
  sub.trend1
  
  # this matches subarea perfectly - good
  
# subtrend and subarea and calcAreaVolLat and sum of AREA_S[X] give the same results (internal consistency)  

  # compare AREA and a ########## THESE were VERY, VERY different but now eqaul.  SB = 0 was the problem
  cbind(AREAice=test[[1]]@data$AREAice, a=subarea[[1]]$AREAice) # shows differences between approaches by polygon 
##################################################################

# compare data sets using the polygons
  # How is AREA calculated - its not using @area
str(test1$c, max.level = 5)

# dates[9]-----
# not gArea
    test1$c@polygons[[1]]@Polygons[[1]]@area/1000000
    str(test1[[3]]@polygons)
# not gArea
    test[[1]]@polygons[[1]]@Polygons[[1]]@area/1000000
    str(test[[1]]@polygons)

# this perfectly replicates the area for subarea    
    p1 <- z@polygons[[1]]@Polygons[[1]]@area/1000000
    p2 <- z@polygons[[1]]@Polygons[[2]]@area/1000000
    p3 <- z@polygons[[1]]@Polygons[[3]]@area/1000000
p1 + p2 + p3    # value for z@data$AREA[1]

test1$c@data$AREA/1000000
test[[1]]@data$AREA/1000000
z@data$AREA


#dates[25]-------
#test[[1]]@polygons[[4]]

# raw data
str(test1$a@polygons[[16]])
test1$a@polygons[[16]]@Polygons[[1]]@area/1000000
test1$a$AREA[16]/1000^2

# polygons have been subsetted but AREA has not been corrected
test1$c@polygons[[4]]@Polygons[[1]]@area/1000000
test1$c@polygons[[4]]@Polygons[[2]]@area/1000000
str(test1$c@polygons[[4]])
test1$c@data$AREA[[4]]/1000^2

# polygons have been subsetted but AREA has not been corrected
test[[1]]@polygons[[4]]@Polygons[[1]]@area/1000000
test[[1]]@polygons[[4]]@Polygons[[2]]@area/1000000
str(test[[1]]@polygons[[4]])
test[[1]]@data$AREA[[4]]/1000^2

# gArea - subsetted and corrected wtih gArea
# this perfectly replicates the area for subarea    
p1 <- z@polygons[[4]]@Polygons[[1]]@area/1000000
p2 <- z@polygons[[4]]@Polygons[[2]]@area/1000000
p1 + p2    # value for 
z@data$AREA[4]

# check to make sure ID numbers match
sapply(slot(test1$c, "polygons"), function(x) slot(x, "ID"))

test1$a@polygons[[16]]@ID
test1$c@polygons[[4]]@ID
test[[1]]@polygons[[4]]@ID
z@polygons[[4]]@ID

#test[[1]]@polygons[[5]]------
# raw data
str(test1$a@polygons[[20]])
test1$a@polygons[[20]]@Polygons[[1]]@area/1000000
test1$a@polygons[[20]]@Polygons[[2]]@area/1000000
test1$a$AREA[20]/1000^2
test1$c$AREA[5]/1000^2

# polygons have been subsetted but AREA has not been corrected
test1$c@polygons[[5]]@Polygons[[1]]@area/1000000
test1$c@polygons[[5]]@Polygons[[2]]@area/1000000
str(test1$c@polygons[[5]])
test1$c@data$AREA[[5]]/1000^2

# polygons have been subsetted but AREA has not been corrected
test[[1]]@polygons[[5]]@Polygons[[1]]@area/1000000
test[[1]]@polygons[[5]]@Polygons[[2]]@area/1000000
str(test[[1]]@polygons)
str(test[[1]]@polygons[[5]])
test[[1]]@data$AREA[[5]]/1000^2

# gArea - subsetted and corrected wtih gArea
# this perfectly replicates the area for subarea 
q1 <- z@polygons[[5]]@Polygons[[1]]@area/1000000
q2 <- z@polygons[[5]]@Polygons[[2]]@area/1000000 # this is a hole
z@polygons[[5]]@Polygons[[2]]@hole
q1 - q2    # value for z@data$AREA[1]

# check to make sure ID numbers match
sapply(slot(test1$a, "polygons"), function(x) slot(x, "ID"))
test1$a@polygons[[20]]@ID # ID 
test1$c@polygons[[5]]@ID
test[[1]]@polygons[[5]]@ID
z@polygons[[5]]@ID


### confirm that subsetting is accurately represented on a map
#PLOT------------------

# load water SPDF and convert to SP so that plot will work
windows()
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)


# LARGE SCALE - plot water and minlongs/lats as points ----
# confirm that filters are working; check ice maps
plot(water, xlim = c(-70, -48), ylim = c(40, 60), col = "lightblue", border = NA) 
plot(filters, border = "red", add = TRUE) # plot main map with filter

yy <- spTransform(test1$a, CRS("+proj=longlat +datum=WGS84")) 

plot(yy, border = "yellow", add = T)

zz <- spTransform(z, CRS("+proj=longlat +datum=WGS84")) 
zzz <- SpatialPolygons((zz@polygons))
str(zz, max.level = 3)
plot(zz, border = "black", add = T)
plot(zzz[4], border = "orange", add = T)

# SMALL SCALE - plot water and minlongs/lats as points ----
windows()
plot(water, xlim = c(-61, -54), ylim = c(49, 57), col = "lightblue", border = NA) 
plot(filters, border = "red", add = TRUE) # plot main map with filter
# letters do not match order in the polygons nor does egg data match what is on the map
plot(yy, border = "black", add = T)
plot(zz, border = "white", add = T)

plot(zzz[1], border = "cyan", add = T)  # C filter?
plot(zzz[2], border = "chartreuse", add = T, lwd=4)  # G  deep 
plot(zzz[3], border = "deeppink", add = T)  # E filter 1
plot(zzz[4], border = "orange", add = T) # F - northern
plot(zzz[5], border = "blue", add = T)  # I
plot(zzz[6], border = "dark green", add = T)  # H
plot(zzz[7], border = "seagreen1", add = T)  # F - southern
plot(zzz[8], border = "maroon1", add = T, lwd=2)  # H
plot(zzz[9], border = "olivedrab", add = T, lwd=2)  # K
plot(zzz[10], border = "yellow", add = T)  # J
plot(zzz[11], border = "yellow", add = T)  # L
plot(zzz[12], border = "cyan", add = T) 
plot(zzz[13], border = "chartreuse", add = T) 
plot(zzz[14], border = "deeppink", add = T) # Northern Penn
plot(zzz[15], border = "cyan", add = T, lwd=2) # north coast
plot(zzz[16], border = "blue", add = T, lwd=4) # 
plot(zzz[17], border = "yellow", add = T, lwd=4) 
plot(zzz[18], border = "deeppink", add = T, lwd=4) 
plot(zzz[19], border = "seagreen1", add = T, lwd=4) 
plot(zzz[20], border = "orange", add = T, lwd=4) 


plot(water, xlim = c(-59, -55), ylim = c(49, 53), col = "lightblue", border = NA)
# tried to plot the holes but it didn't work

# plot polygon with specific ID
class(water)
class(filters)
class(zz)
methods(class = "Spatialpolygons")
methods("plot")
#################################### 
# show all areas-----

# show all areas
a <- sapply(slot(test1$a, "polygons"), function(x) slot(x, "area"))/1000^2

b <- sapply(slot(test1$c, "polygons"), function(x) slot(x, "area"))/1000^2

c <- sapply(slot(z, "polygons"), function(x) slot(x, "area"))/1000^2


bc <- as.data.frame(cbind(b, c))

abc <- left_join(a, b)

str(test1$a@polygons[[9]])
test1$a@polygons[[9]]@Polygons[[1]]@area/1000000
test1$a$AREA[9]/1000^2

# polygons have been subsetted but AREA has not been corrected
test1$c@polygons[[4]]@Polygons[[1]]@area/1000000
test1$c@polygons[[4]]@Polygons[[2]]@area/1000000
str(test1$c@polygons[[1]])
test1$c@data$AREA[[4]]/1000^2

# polygons have been subsetted but AREA has not been corrected
test[[1]]@polygons[[4]]@Polygons[[1]]@area/1000000
test[[1]]@polygons[[4]]@Polygons[[2]]@area/1000000
str(test[[1]]@polygons[[4]])
test[[1]]@data$AREA[[4]]/1000^2

# gArea - subsetted and corrected wtih gArea
# this perfectly replicates the area for subarea    
p1 <- z@polygons[[4]]@Polygons[[1]]@area/1000000
p2 <- z@polygons[[4]]@Polygons[[2]]@area/1000000
p1 + p2    # value for 
z@data$AREA[4]

# check to make sure ID numbers match
sapply(slot(test1$a, "polygons"), function(x) slot(x, "ID"))
test1$a@polygons[[9]]@ID
test1$c@polygons[[1]]@ID
test[[1]]@polygons[[1]]@ID
z@polygons[[1]]@ID

###############################
# Eric's code----
want_bin = test[[1]]$`BIN#`[6:10]
> want_bin
length(test1$c)
[1] 14
> length(want_poly)
plot(test1$c[want_poly,])
> plot(water, xlim = c(-70, -48), ylim = c(40, 60), col = "lightblue", border = NA)
> plot(test1$c[want_poly,], add=T)
> plot(test1$c[want_poly,], add=T,col="red")
> a = test1$c[want_poly,]
> a$AREA
[1]   43918119  130360938  404375186 3606216830   31709556
> a$AREA/1000^2
[1]   43.91812  130.36094  404.37519 3606.21683   31.70956
> plot(test1$c[want_poly,], add=T,col="red")
> plot(a, add=T)

############################################
# compare AREA and area in dataframe and polygons while working through the funcitons.  This is what showed that SB was being subsetted improperly------
z <- dates3
i <- 9
dates3[25]

rm(ice)
rm(egg)
rm(sub.egg)
rm(sub.egg1)
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3-test.R")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

# 9604
ice@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
ice@data[17, c("AREA", "BIN#", "A_LEGEND", "ID")]
#ice@data$AREA[9]/1000^2
str(ice@polygons[17])
ice@data[17,]

egg@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
egg@data[3, c("AREA", "BIN#", "A_LEGEND", "ID")]
#egg@data$AREA[6]/1000^2
str(egg@polygons[3])

# reduction in polygons is due to logical issues with SB - OK for now
sub.egg@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
sub.egg@data[3, c("AREA", "BIN#", "A_LEGEND", "ID")]
#sub.egg@data$AREA[3]/1000^2
str(sub.egg@polygons[3])

#subset
sub.egg1@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
sub.egg1@data[1, c("AREA", "BIN#", "A_LEGEND", "ID")]
#sub.egg1@data$AREA[3]/1000^2
str(sub.egg1@polygons[1])

#transform
sub.egg1@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
sub.egg1@data$AREA[3]/1000^2
str(sub.egg1@polygons[3])

z@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
z@data$AREA[1]
str(z@polygons[1])

# eggATTR.cols
rm(x)
rm(y)
x <- sub.egg1

x@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
x@data$AREA[3]
str(x@polygons[3])
x@data$AREA[1]
str(x@polygons[1])

#eggAttr.query - here is where the problems start!!! 10-filters becomaes 5-filters
# look within query
# had to work through these with the browser but it showed that some filters were being subsetted and that SB=0 was the problem
x@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
x@data$AREA[1]
str(x@polygons[1])

# within eggAttr.query - but this gives the right answer
x@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
x@data$AREA[1]
str(x@polygons[1])


y@data[c("AREA", "BIN#", "A_LEGEND", "ID")]
y@data$AREA[1]
str(y@polygons[1])

x <- sub.egg1

##########################################################################
rm(spdf)
spdf <- extractSPDFfinal(dates3, 25, ct=m1$ct, sa=m1$sa, sb=m1$sb)

AREA <- sum(spdf[[1]]$AREA)/1000000
iceAREA <- sum(spdf[[1]]$AREAice)/1000000

# calculate total area from AREA of each ice type
A <- sum(spdf[[1]]$AREA_SA)/1000000
B <- sum(spdf[[1]]$AREA_SB)/1000000
C <- sum(spdf[[1]]$AREA_SC)/1000000
D <- sum(spdf[[1]]$AREA_SD)/1000000

TOT <- sum(A, B, C, D)

iceAREA/TOT # this compares favourably to AREA 

# look at AREA and AREA[SX]
compare <- spdf[[1]][, c("AREA", "AREAice", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")]
# sum the columns of AREA[SX]
sumcols <- (spdf[[1]]$AREA_SA + spdf[[1]]$AREA_SB + spdf[[1]]$AREA_SC + spdf[[1]]$AREA_SD)

# compare the values of AREA and sumcols - should be equal 
# note that this will change - sumcols should be a percentage of AREA corresponding to CT
as.data.frame(cbind(compare, sumcols))
as.data.frame(cbind(compare[,c("AREA", "AREAice")], sumcols))
# These match suggesting that the AREA_SX is summing correctly (internally consistent)
sum(sumcols)
sum(compare[, c("AREA")])
sum(compare[, c("AREAice")])
