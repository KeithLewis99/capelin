# The question here is, what is the right area? calcAreaVolLat() faithfully reproduces an areas calculation but is it right?  I have tried to take this function apart to reproduce the values.  Unfortunately, its always different.  Why?  
# things seem to be changing from day to day - things that were the same yesterday aren't today.  Could be an environment thing??



rm(list=ls())
load("output-processing/dates3.Rdata")
load("output-processing/filters.Rdata")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")

y <- dates3
i <- 1
dates3[25]


# extract essence of sub.egg1
test <- extractSubegg1(dates3, 25, ct=m1$ct, sa=m1$sa, sb=m1$sb)
test[[1]]@data
sub <- test[[1]]@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")]
str(sub)


#extract all sub.egg - remember that ice is unfiltered and will have lots more polygons
test1 <- lookAtIce(dates3, 25, ct=m1$ct, sa=m1$sa, sb=m1$sb)

test1$c@data # this has filtered out irrelevant polygons
orig <- test1$c@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD")]
str(orig)

# calculate total area from AREA for sub.egg1
AREA <- sum(test[[1]]@data$AREA)/1000000

# calculate total area from AREA of each ice type
A <- sum(test[[1]]@data$AREA_SA)/1000000
B <- sum(test[[1]]@data$AREA_SB)/1000000
C <- sum(test[[1]]@data$AREA_SC)/1000000
D <- sum(test[[1]]@data$AREA_SD)/1000000

TOT <- sum(A, B, C, D)

AREA/TOT # this compares favourably to AREA 

# look at AREA and AREA[SX]
compare <- test[[1]]@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "AREA_SA", "AREA_SB", "AREA_SC", "AREA_SD")]
# sum the columns of AREA[SX]
sumcols <- (test[[1]]@data$AREA_SA + test[[1]]@data$AREA_SB + test[[1]]@data$AREA_SC + test[[1]]@data$AREA_SD)
sum(sumcols)

# compare the values of AREA and sumcols - should be equal 
# note that this will change - sumcols should be a percentage of AREA corresponding to CT
as.data.frame(cbind(compare, sumcols))
as.data.frame(cbind(compare[,c("AREA")], sumcols))
# These match suggesting that the AREA_SX is summing correctly (internally consistent)


#why does sum of AREA not equal areas???? Externally inconsistent
#pick up where iceSubset leaves off
# this calculates area
sub.poly <- calcPolyA(test1$c) # calcPoly calculates the area of each polygon manually using either 1) gArea or 2) mannualy calculating the area from all the slot @area in each polygon

str(sub.poly)
sub.trend <- trendsCalc(test1$c, sub.poly) # converts the AREAS in the SPDF and sums all areas

# this further breaks down the code to run it line by line
z <- test1$c

  a <- try(gArea(z, byid = TRUE)) 
  str(a)
# from iceAREA  
  z@data$AREA <- a * 1e-6
  subarea <- iceArea(z@data)
  options(scipen=999)  

# sub.trend and subarea are equal    
  
  # calculation of areas does NOT match subarea (has been through all)
  calcAreaVolLat(dates3[25], ct=m1$ct, sa=m1$sa, sb=m1$sb)  
  sub.poly1 <- calcPolyA(test[[1]]) 
  str(sub.poly1)
  sub.trend1 <- trendsCalc(test[[1]], sub.poly1)

  # sum of AREA (has been through iceArea but not gArea)
  sum(test[[1]]@data$AREA/1000000) # this value is teh same as AREA - good!!!!! 14 values
  sum(test[[1]]$AREA_SA + test[[1]]$AREA_SB + test[[1]]$AREA_SC + test[[1]]$AREA_SD)/1000000
  
  # sum of AREAS[SX] - has been through all 
  sum(subarea[[1]]$AREA_SA + subarea[[1]]$AREA_SB + subarea[[1]]$AREA_SC + subarea[[1]]$AREA_SD)
  sum(subarea[[1]]$AREA)                  # 19 values
  sub.trend
  test[[2]]/1000000
  # this matches subarea perfectly - good
  
# subtrend and subarea and calcAreaVolLat and sum of AREA_S[X] give the same results (internal consistency)  

  

  # compare AREA and a ########## THESE ARE VERY, VERY different - why????
  cbind(AREA=test[[1]]@data$AREA/1000000, a=subarea[[1]]$AREA) # shows differences between approaches by polygon 
  
    
##################################################################

# compare data sets using the polygons
  # How is AREA calculated - its not using @area
str(test1$c, max.level = 5)
str(test[[1]], max.level = 5)
str(z)

# strucutre of polygons
str(test1[[1]]@polygons[4])
str(test1[[1]]@polygons[[5]]@Polygons)
test1$c@polygons[[10]]@Polygons

# dates[9]
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


#dates[25]
# not gArea
test1$c@polygons[[1]]@Polygons[[1]]@area/1000000
str(test1[[3]]@polygons)
# not gArea
test[[1]]@polygons[[1]]@Polygons[[1]]@area/1000000
str(test[[1]]@polygons)

# this perfectly replicates the area for subarea    
p1 <- z@polygons[[4]]@Polygons[[1]]@area/1000000
p2 <- z@polygons[[4]]@Polygons[[2]]@area/1000000
p1 + p2    # value for 
z@data$AREA[4]

str(z@polygons)
test1$c@data$AREA/1000000
test[[1]]@data$AREA/1000000
z@data$AREA

q1 <- z@polygons[[5]]@Polygons[[1]]@area/1000000
q2 <- z@polygons[[5]]@Polygons[[2]]@area/1000000
q1 - q2    # value for z@data$AREA[1]

z@data$AREA[5]
test[[1]]@data$AREA/1000000
z@polygons[[5]]@Polygons[[2]]@hole

#area of polygons is equal but the AREA is not equal bc in z, the AREA is computed based on the values of the polygons....so what is right?  Need to see it spatially

# this is a spatial display of the polygons with water, ice, filters, subset of ice and each polygon

# load water SPDF and convert to SP so that plot will work
load("sp_data/20160606.Rdata") # file doesnot exist - probably doesn't matter
water <- spTransform(ice[ice$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
water <- gUnaryUnion(water)

# plot water and minlongs/lats as points 
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


plot(water, xlim = c(-61, -54), ylim = c(50, 55), col = "lightblue", border = NA) 
plot(filters, border = "red", add = TRUE) # plot main map with filter
# letters do not match order in the polygons nor does egg data match what is on the map
plot(yy, border = "yellow", add = T)
plot(zz, border = "black", add = T)

plot(zzz[1], border = "grey", add = T)  # C filter?
plot(zzz[2], border = "red", add = T)  # ???
plot(zzz[3], border = "pink", add = T)  # E filter 1
plot(zzz[4], border = "orange", add = T) # F - northern
plot(zzz[5], border = "blue", add = T)  # I
plot(zzz[6], border = "dark green", add = T)  # H
plot(zzz[7], border = "orange", add = T)  # F - southern
plot(zzz[8], border = "orange", add = T)  # H
plot(zzz[9], border = "pink", add = T)  # K
plot(zzz[10], border = "orange", add = T)  # J
plot(zzz[11], border = "orange", add = T)  # L

# holes cannot be plotted - this indicates whether there is a hole in the polygon or not
plot(zzz[5]@Polygons[[2]]@hole, border = "blue", add = T) # nonsense
zzz[5]@Polygons[]

z@data
# plot polygon with specific ID
class(water)
class(filters)
class(zz)
methods(class = "Spatialpolygons")
methods("plot")
#################################### How do I get ID's to match up?????  
str(z@polygons[[1]]@ID)
z@polygons[[1]]@ID
z@polygons[[1]]
z@data
str(test1[[1]]@polygons[[5]]@Polygons)
test1[[1]]@polygons[[5]]@Polygons




#################################### Junk
test[[1]]@data$AREA/1000000


p12 <- test1[[1]]@polygons[[10]]@Polygons[[2]]@area

plot(test1[[10]])
(p11 + p12 )/1000000
test[[1]]@data$AREA[10]/1000000

#test - this is Subegg
test$c
str(test1$c, max.level = 5)
str(test1$c@polygons[1])
str(test1$c@polygons[[1]]@Polygons)

# this recreates the value for AREA[10]!!!!!!!!!!!
p11 <- test1$c@polygons[[10]]@Polygons[[1]]@area
p12 <- test1$c@polygons[[10]]@Polygons[[2]]@area


(p11 + p12)/1000000
test1$c@data$AREA[10]/1000000


sub.egg1@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD")]
z <- sub.egg1
y <- a
x <- sub.egg1
x@data[, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD")]
subarea[[1]][, c("AREA", "CT", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD")]



#####################junk
x <- test1$c
x@data
wtf <- withinPolyAreaC(x)
wtf@data
wtf@data[1, ]

x$AREA_SC[i] <- x$AREA[i] * as.numeric(x$CC[i])/as.numeric(x$CT[i])

x$AREA_SC[i] <- x$AREA[i] * as.numeric(x$CC[i])/as.numeric(x$CT[i])



Sr1 = Polygon(cbind(c(2,4,4,1,2),c(2,3,5,4,2)))
Sr2 = Polygon(cbind(c(5,4,2,5),c(2,3,2,2)))
Sr3 = Polygon(cbind(c(4,4,5,10,4),c(5,3,2,5,5)))
Sr4 = Polygon(cbind(c(5,6,6,5,5),c(4,4,3,3,4)), hole = TRUE)
Srs1 = Polygons(list(Sr1), "s1")
Srs2 = Polygons(list(Sr2), "s2")
Srs3 = Polygons(list(Sr3, Sr4), "s3/4")
SpP = SpatialPolygons(list(Srs1,Srs2,Srs3), 1:3)
str(SpP)
# plot single polygon
par(mfrow=c(3,1))
plot(SpP[1])
plot(SpP[3])

# or using IDs: retrieve list of all IDs
IDs = sapply(SpP@polygons, function(x) x@ID)

# plot polygon with specific ID
plot(SpP[which(IDs == 's2')])

