################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# This file was created simply to test the outputs of ice-chart-processing-function-v*.  It does not do anything perse but was used to explore how to extract objects from a Spatial Polygon Data frame, how to subset a SPDF, and how to use regular expressions to help with this process.

## Structure of SPDF-----
# based on the output of funciton calcAreaVolume where a single SPDF is extracted, the below is the code to explore the structure 
#sturcutre of data - various ways to drill into the data set
slotNames(sub.egg)
str(sub.egg)
str(sub.egg@data) # see data structure
str(sub.egg@data$EGG_ATTR) # see EGG_ATTR
str(sub.egg@data$E_CT)

#structure of polygons
str(sub.egg@polygons)
str(sub.egg@polygons[1])
str(sub.egg@polygons[[1]]@Polygons)
str(sub.egg@polygons[[1]]@Polygons[1])
str(sub.egg@polygons[[1]]@Polygons[[1]])
str(sub.egg@polygons[[1]]@Polygons[[1]]@area)

# structure at various levels - prevents huge printouts with loads of redundancy
plot(test)
test <- water@polygons  
str(water)
str(water, max.level = 2)
str(water, max.level = 3)
str(water, max.level = 4)
str(water, max.level = 5) # shows relevant structure, i.e, lenght of list
str(water, max.level = 6) 
str(water, max.level = 7)
str(water, max.level = 8) # same as str(water)

# see data
sub.egg@data
sub.egg@data$EGG_ATTR
sub.egg@data$E_CT
sub.egg@data$DATE_CARTE

# S3 or S4 check
isS4(sub.egg)
isS4(sub.egg@data)

# check of SBDF - compare EGG_ATTR and columns
test <- sub.egg@data[, c("AREA", "A_LEGEND", "EGG_ATTR", "E_CT", "E_CA", "E_CB", "E_SO", "E_SA", "E_SB", "E_SC")]
test
head(test)

test <- x[, c("AREA", "A_LEGEND", "EGG_ATTR", "E_CT", "E_CA", "CT", "E_SA", "SA")]
str(test)
head(test)

## Regular expressions and SPDF------
# grep, grepl

#create some data
x <- c("7", "7.", "4", "4.")
x

y <- c("6", "8.", "4", "4.")
y

xx <- c("7", "7", "4", "4")
xx
str(xx)

# make string
z <- (cbind(xx, y))
z
str(z) 

# make dataframe
zz <- as.data.frame(cbind(xx, y), stringsAsFactors = F)
zz
str(zz)

# extract parts of string
grep("7", x)
grep("7.", x)
grep(".", x)

# logical
grepl("7", x)
grepl("7.", x)
grepl(".", x)

# extract parts of string
grep("7", z)
grep("8.", z)
grep(".", z)

# logical
grepl("7", z)
grepl("8.", z)
grepl(".", z)
grepl(" .", z)
grepl("..", z)

# testing wildcards
test <- grepl("..", zz)
str(test)
test
zz
test1 <- as.data.frame(grepl("..", zz))
str(test1)
test1

# selecting across rows and columns
for(i in 1:nrow(zz)){
  for(j in 1:ncol(zz)){
    if(grepl("..", zz)){
      print(paste(i = "yeah", j="duh", sep=","))    
    }
  }  
}

for(i in 1:nrow(zz)){
  for(j in 1:ncol(zz)){
    if(grepl("7", zz)){
      print("yeah")    
    } else {
      print("duh")
    }
  }  
}

## subset SPDF------
#subset data - change as per discussions with others
# WARNING: THE BELOW WAS SIMPLY TO LEARN HOW TO SUBSET.  HOWEVER, THE USE OF THE "E_xx" DATA IS NOT VALID BC THE "FAST ICE" IS NOT INCLUDED - IT RESIDES ONLY IN THE EGG_ATTR TABLE

# create vectors for subsetting
r <- c("10", "9+", "9", "8", "7", "6", "5", "4")
r <- c("9+", "10")
r

# this subsets the dataframe 
test <- subset(sub.egg, E_CT %in% r) 
str(test@data)
test$E_CT
head(test[, c("AREA", "A_LEGEND", "EGG_ATTR", "E_CT", "E_CA", "E_CB", "E_SO", "E_SA", "E_SB", "E_SC")])
test@data[, c("AREA", "A_LEGEND", "EGG_ATTR", "E_CT", "E_CA", "E_CB", "E_SO", "E_SA", "E_SB", "E_SC")]

# subsets as above
test2 <- subset(sub.egg, E_CT == "9+"| E_CT == "10") 
str(test2@data)

# subsets as above
test3 <- subset(sub.egg, E_CT == r) 
str(test3@data)
test3$E_CT



# generalize to a range of columns-------------------
# problem is that we need to loop over z, w, and i.  Perhaps this can be done with mapply but beyond me for now.
x[[w[1]]]
withinPolyArea_GEN <- function(x = data, name = name, cols = z, conc = w, sa = sa) {
    for (i in 1:nrow(x)) {
      for(k in z){
        if(is.na(x[[z[k]]][i])) {
        x[,paste("AREA_", z[k], sep='')] <- 0     # if value of SA = NA, then set area to zero
      } else {
        for(j in w){
          if (x[[z[k]]][i] %in% sa) {
            x[,paste("AREA_", z[k], sep='')] <- x$AREA[i] * as.numeric(x[[w[j]]][i])/10 # calculate area for desired values of SA
          } else {
            x[,paste("AREA_", z[k], sep='')] <- 0 # if SA is not of desired value (but not NA) - then area is 0
          }
        }
      }
      }
    }
  return(x)
  }


# generalize to a range of columns - maaply-------------

withinPolyArea_GEN <- function(x = data, name = name, cols = z, conc = w, sa = sa) {
  for (i in 1:nrow(x)) {
    for(k in z){
      k <- paste
      if(is.na(x[[]][i])) {
        x[,paste("AREA_", z[k], sep='')] <- 0     # if value of SA = NA, then set area to zero
      } else {
        for(j in w){
          if (x[[z[k]]][i] %in% sa) {
            x[,paste("AREA_", z[k], sep='')] <- x$AREA[i] * as.numeric(x[[w[j]]][i])/10 # calculate area for desired values of SA
          } else {
            x[,paste("AREA_", z[k], sep='')] <- 0 # if SA is not of desired value (but not NA) - then area is 0
          }
        }
      }
    }
  }
  return(x)
}

test5 <- mapply(sum, df$a, df$b)
test5

test5 <- mapply(calc(df, "a", "b", "d"))
test5


###mapply#######--------------
withinPolyAreaALL <- function(x = data, sa = sa, conc = w) {
  for (i in 1:nrow(x)) {
    if (is.na(z[i])) {
      x[,paste("AREA_", z[k], sep='')] <- 0   # if value of SA = NA, then set area to zero
    } else {
      if (z[i] %in% sa) {
        x[,paste("AREA_", z[k], sep='')] <- x$AREA[i] * as.numeric(w[[j]][i])/10 # calculate area for desired values of SA
      } else {
        x[,paste("AREA_", z[k], sep='')] <- 0 # if SA is not of desired value (but not NA) - then area is 0
      }
    }
  }
  return(x)
}

withinPolyAreaALL <- function(x = data, sa = sa) {
    if (is.na(z[i])) {
      x[,paste("AREA_", z[k], sep='')] <- 0   # if value of SA = NA, then set area to zero
    } 
  return(x)
}

#############################################
test2 %>%
  mutate_at(.vars = z,
            .funs = withinPolyAreaALL)

test4 <- mapply(withinPolyAreaALL(test2, sa, "w"), test2$SA, test2$SB)

x[[w[1]]][1]
x[[z[1]]]
x[[z[i]]]
x[[z[i]]][1]
z[k]
z[i]
x$AREA[3]
x[[w[j]]][3]

x[[z[k]]][3]
str(x[[z[k]]])


test3 <- withinPolyAreaA(test2, sa)
test3
test3$AREA_NEW
test2$SA
test3$SA
test3[c("AREA", "CA", "SA", "SB", "AREA_NEW")]
AREA3 <- test2$AREA*(as.numeric(test2$CB)/10)

test4 <- withinPolyArea_GEN(test2, sa)
test4
test4[c("AREA", "SA", "SB", "AREA_NEW")]

x <- test$a
x@data
test
test$a
test$a@data
iceArea(test$a@data, m1$ct, m1$sa)

# for testing e00_to_SpatialPolygonDataframe()----------------------
#i = "19730528"
#i = "19780122"
#i = "19780423"
#i="19830102"
#i="19840313"
#i="19900226"
#i = "20170403"

#head(dates)

#setwd("D:/Keith/ice/2017-project/e00_data")
#test <- raster('20170626.e00')
#test <- readGDAL('20170626.e00')

# import the e00 data into R: all looks well ito of egg attribute - note HUGE number of polygons
shape1 <- readOGR(dsn = "D:/Keith/ice/2017-project/e00_data_conversion", layer = "test3")
str(shape1)
str(shape1@data)
shape1@data
class(shape1)
write.csv(shape1@data, "shapeflies.csv")
out <- capture.output(str(shape1))
write.csv(out, "str_test.csv")

# make shape1 a SPDF - shouldn't need this.
shape.spdf <- SpatialPolygonsDataFrame(shape1, data = shape1@data)
plot(shape.spdf)
save(shape1, file = "D:/Keith/ice/2017-project/e00_data_conversion/20170213.Rdata")

# calculate area using calAreaVolume - use data set that wouldn't load [1] and one that would [2]
dates4 <- as.POSIXlt(c("2017-02-13 NST", "2017-02-20 NST"))
y <- dates4
i <- 1
y[i]

load(format(y[i], "sp_data/%Y%m%d.Rdata"))
ice1 <- load(format(y[i], "sp_data/%Y%m%d.Rdata"))
plot(ice)
plot(ice1)
plot(shape1)

trends_update4 <- calcAreaVolume(dates4, ct=m1$ct, sa=m1$sa)


#################################################################
# area calculations
test <- lookAtIce(dates3, 626) # does not work for 1
str(test$a@data)
str(test$b)
head(test$b, 30)
str(test$a, max.level = 2)
class(test$a)
test$a@data

gl <- iceArea(test$a@data, ct=m3$ct, sa=m3$ct)
test$b


m1$ct
m1$sa
gl[[1]]
gl[[1]][, c("AREA", "EGG_ATTR", "CT", "CA", "CB", "SA", "SB", "AREA_SA", "AREA_SB")]
gl[[1]][7,][c("AREA", "EGG_ATTR", "CT", "CA", "CB", "SA", "SB", "AREA_SA", "AREA_SB")]
gl[[2]]


sum(gl[[1]]$AREA_SA)
area <- gl[[1]]$AREA_SA + gl[[1]]$AREA_SB + gl[[1]]$AREA_SC + gl[[1]]$AREA_SD
etab <- eggAttr(gl[[1]]$EGG_ATTR)
coverage <- as.numeric(gsub("\\+", ".5", etab$E_CT))/10


a <- gl[[1]][7,]$AREA
ca <- as.numeric(gl[[1]][7,]$CA)
cb <- as.numeric(gl[[1]][7,]$CB)
asa <- a*(ca/10)
asb <- a*(cb/10)
ar <- asa + asb
etab <- eggAttr(gl[[1]][7,]$EGG_ATTR)
coverage <- as.numeric(gsub("\\+", ".5", etab$E_CT))/10 

icesum <- sum(ar * coverage, na.rm = TRUE)

x$AREA_SA[i] <- x$AREA[i] * as.numeric(x$CA[i])/10


###Look at coordinates, boundary boxes, and extents---------------------
str(x)
x@polygons[[1]]
head(x@polygons[[1]]@Polygons[[1]]@coords)
iceTiming(test$a)

plot(x@bbox)
class(x)
x@proj4string
coordinates(x)

ie <- extent(ice)
ibb <- ice@bbox
plot(ie)
se <- extent(temp.ls[[1]])
plot(se, add=T)
intersect(ie, se)
gIntersects(ie, se)

extent(ice)
extent(sub.egg)
extent(sub.egg1)
extent(temp.ls[[1]])

y <- spTransform(x[x$A_LEGEND != "Land", ], CRS("+proj=longlat +datum=WGS84"))
ptest <- gUnaryUnion(test)
plot(test)
eastcoast <- gIntersection(water, test, byid = TRUE) 




###Work with environments---------------------
exists("i", globalenv())
exists("water", globalenv())
sys.frame()
ls()
ls.str()
library(pryr)
where("calcAreaVolLat")
where("calcLatLong")




