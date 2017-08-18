exists("i", globalenv())
exists("water", globalenv())
sys.frame()
ls()
ls.str()
library(pryr)
where("calcAreaVolLat")
where("calcLatLong")



calcAreaVolLat <- function(y, ct = NULL, sa = NULL){ 
  #browser()
  #my.env <- new.env()
  areas <- volumes <- minlats <- rep(NA, length(y)) # create an empty object
  for(i in seq_along(y)) {
    ## load map data
    print(y[i])                   #start here when making single object for testing
    load(format(y[i], "sp_data/%Y%m%d.Rdata"))
     ## set area and volume to 0 if no egg attributes
    #print("inside cAVL")
    #print(where("ice"))
    #crs <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +datum=NAD27")
    #ice <- ice
    if(all(is.na(ice$EGG_ATTR))) {  # START HERE
      area <- volume <- 0
      minlat <- 55 
    } else {
    egg <- subsetProject(ice)    # ice is not in the local environment
    sub.egg <- filterEgg(egg)
    temp.ls <- calcLatLong(sub.egg, ct=ct, sa=sa)
    temp.ls <- calcPolyA(temp.ls)
    calc <- trendsCalc(temp.ls)
    }
    areas[i] <- calc$area #???
    volumes[i] <- calc$volume
    minlats[i] <- calc$minlat
    minlongs[i] <- calc$minlong
  }
  list(areas = areas, volumes = volumes, minlats = minlats, minlongs = minlongs) 
}

test <- calcAreaVolLat(dates3[20], ct=m1$ct, sa=m1$sa)




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

temp.ls$sub.egg1@data
sub.egg1@data
sub.egg@data
egg@data

rm(ice)
rm(ice1)
rm(egg)
rm(sub.egg)
rm(sub.egg1)
rm(temp.ls)
rm(calc)
rm[i]
rm[y]
rm(ct)
rm(sa)
rm(coords_min)

#iceTiming(sub.egg)
y <- dates3[10:20]
i=11

y[i]
ct=m1$ct 
sa=m1$sa

# trying to get this to work without the if/else stuff going on
calcLatLong <- function(sub.egg, ct = NULL, sa = NULL){
  #browser()
  print("inside cLL")
  #print(environment())
  #print(parent.env())
  #search()
  #list(e = environment(), p = parent.env(environment()))
  #sys.frames()
  #print(ls())
  
#  e <- environment(calcLatLong)
 # p <- parent.env(calcAreaVolLat())
  #ice <- p$ice
  print(environment())

  
    sub.egg1 <- iceSubset(sub.egg, ct = ct, sa = sa) 
    coords_min <- iceTiming(sub.egg1)
  
      minlat <- coords_min$lat
      minlong <- coords_min$long
      #ice1 <- assign("ice", value = ice, pos=2) #doesn't do much
      #ice1 <- c(.GlobalEnv, ice)
      #ice1 <<- ice
      #ice1 <- assign("ice", ice, envir=parent.frame())
      #get0(ice, envir=
      #eval(quote(ice), parent.frame())
        sub.egg1 <- spTransform(sub.egg1, CRS(proj4string(ice)))
  
  return(list(sub.egg1=sub.egg1, coords_min=coords_min))
}

environment(calcLatLong)
str(calcLatLong())
temp.ls <- calcLatLong(sub.egg)
 

globalenv()

test <- iceTiming(sub.egg1)

x <- temp.ls
x$sub.egg1@data
seq_along(y)

test <- function(z){
  zz <- rep(0, length(z))
  for(i in seq_along(z)) {
    zz[i] <- z[i] + 1
    print(zz)
  }
}

z <- as.data.frame(1:10)

test(z)

loadMap <- function(y){
  ## load map data
  for(i in seq_along(y)){
    print(y[i])                   #start here when making single object for testing    
  }

 # load(format(y[i], "sp_data/%Y%m%d.Rdata"))
  #return(ice)
  #return(list(test=ice))
}
loadMap(y)


h <- function() {
  x <- 10
  function() {
    x
  }
}
i <- h()
x <- 20
i()
