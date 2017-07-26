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

## test subsetting of the AreaVolume function------------------------------------------------------
## subset the ice data
# create "model sets"
ct <- c(3:10)
sa <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
m1 <- list(ct=ct, sa=sa)
m1

ct <- c(8, 9)
sa <- c("5")
m2 <- list(ct=ct, sa=sa)
m2

ct <- NULL
sa <- NULL
m3 <- list(ct=ct, sa=sa)
m3

#lookAtSubEgg function-------------------
test <- lookATSubEgg(dates3, 71)
str(test$a@data)
str(test$b)
head(test$b, 30)
str(test$a)

# test eggAttr.within()
test1 <- eggAttr.cols(test$b)
head(test1)
test1[c("AREA", "CA", "CB", "CC", "CD", "SA", "SB", "SC", "SD", "SE")]

test2 <- eggAttr.query(test1, m3$ct, m3$sa)
head(test2)
head(str(test2))
test2
#eggAttr.within(test2)

i=4
k=1
j=1
z <- c("SA", "SB")
w <- c("CA", "CB")
name <- "x"
x <- test2
x$AREA_SA
test2$AREA_SA
str(x)
# generalize to a range of columns
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

 
test4 <- withinPolyArea_GEN(test2, "z", "w", sa)
test4 <- mapply(withinPolyArea_GEN(test2, sa) SA, SB)


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


x <- test$a
x@data
test
test$a
test$a@data
iceArea(test$a@data, m1$ct, m1$sa)
