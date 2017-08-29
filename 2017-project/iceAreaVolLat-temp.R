
calcAreaVolLat <- function(y, ct = NULL, sa = NULL){ 
  areas <- volumes <- minlats <- rep(NA, length(y)) # create an empty object
  for(i in seq_along(y)) {
    ice <- loadMap(y)  
    egg <- subsetProject(ice)
    l <- filterEgg(egg)
    a <- calcPolyA(l$sub.egg)
    calc <- trendsCalc(a)
    areas[i] <- calc$area #???
    volumes[i] <- calc$volume
    minlats[i] <- calc$minlat
     }
  list(areas = areas, volumes = volumes, minlats = minlats) 
}

outt <- calcAreaVolLat(dates3[1:10], ct=m1$ct, sa=m1$sa)



############################################################
rm(ice)
rm(y)
seq_along(y)
y <- dates3
i = 20
y[i]
rm(i)

rm(ice)

#######################################################
## load map data
loadMap <- function(y){
  print(y[i])                   #start here when making single object for testing
  load(format(y[i], "sp_data/%Y%m%d.Rdata"))
  return(ice)
}

#ice <- loadMap(y)

############################################################
subsetProject <- function (ice){
## subset and project to WGS84
  #browser()
    egg <- ice[!is.na(ice$EGG_ATTR), ]  # subset
    egg <- spTransform(egg, CRS("+proj=longlat +datum=WGS84")) # convert to WGS84 bc filters are in WGS84 - used on line 460
  return(egg)
}

#egg <- subsetProject(ice)

############################################################
## Apply filters and return to lcc projection
filterEgg <- function(egg) {
  #browser()
  attr.tab <- egg@data
  sub.egg <- try(gDifference(egg, filters, byid = TRUE))
  if(class(sub.egg) == "try-error") {
    buf.egg <- try(gBuffer(egg, byid = TRUE, width = 0)) # adding a small buffer often fixes geo issues
    sub.egg <- try(gDifference(buf.egg, filters, byid = TRUE)) # remove the egg areas with filters
  }
  if(class(sub.egg) != "try-error" & !is.null(sub.egg)) {
    row.names(attr.tab) <- paste(row.names(attr.tab), "- filters")
    attr.tab <- attr.tab[names(sub.egg), ]
    sub.egg <- SpatialPolygonsDataFrame(sub.egg, attr.tab) # recover data bc of gBuffer and gDifference which make a SP rather than an SPDF
  }
  if(class(sub.egg) == "try-error") { stop("There was a problem applying the filters") }
  return(sub.egg=sub.egg)
}


#sub.egg <- filterEgg(egg)      


############################################################
## proceed with area and volume calculations
calcLatLong <- function(sub.egg, ct = NULL, sa = NULL) {
  #browser()
  if(is.null(sub.egg)) {
    area <- volume <- 0
    minlat <- 55 # if NULL, then no ice was below 55 deg North
  } else {
      if(class(sub.egg) != "try-error") {
      sub.egg1 <- iceSubset(sub.egg, ct = ct, sa = sa) 
      
      #if(nrow(sub.egg1@data)==0){ # if you subset out all data
       # area <- volume <- 0
        #minlat <- 55
        #sub.egg1 <- sub.egg
      }
    #return(sub.egg1)
    }
      coords_min <- iceTiming(sub.egg1)
      minlat <- coords_min$lat
      minlong <- coords_min$long
      
      ## return to lcc projection and calculate area and volume bc e00 data is in lcc and meters - needed to get proper area/volume calculations
      #ls(envir = parent.frame())
      ice <- get("ice", parent.frame())
      print(proj4string(ice))
      print(proj4string(sub.egg1))
      #ice@bbox <- sub.egg1@bbox
      sub.egg1 <- sp::spTransform(sub.egg1, CRS(proj4string(ice)))
      #sub.egg1 <- sp::spTransform(sub.egg1, CRS("+proj=lcc"))
      
      #return(list(sub.egg1=sub.egg1, coords_min=coords_min))
      return(list(sub.egg1=sub.egg1, minlat=minlat, minlong=minlong))
}

#environment(temp.ls)
#temp.ls <- calcLatLong(sub.egg, ct = m1$ct, sa = m1$sa)
############################################################
## MAKE A FOLDER FOR THE FILTERED DATA???
############################################################
calcPolyA <- function(z){
  #browser()
  a <- try(gArea(z$sub.egg1, byid = TRUE)) # sometimes holes are not identified correctly, so try and extract max polygon area within each id
  if(class(a) == "try-error") { 
    message("gArea didn't work. Trying alternate appriach")
    a <- sapply(slot(z$sub.egg1, "polygons"), function(x) max(sapply(slot(x, "Polygons"), slot, "area"))) 
  }
  #z$sub.egg1$AREA <- a * 1e-6 # replace polygon area (use square km)
  #temp.ls$a <- a
  z$a <- a
  return(z=z)
}

#temp.ls <- calcPolyA(temp.ls)
#temp.ls[[2]]$sub.egg1@data

############################################################
trendsCalc <- function(x){
  #browser()
  #x$sub.egg1@data$AREA_SA <- x$sub.egg1@data$AREA_SB <- x$sub.egg1@data$AREA_SC <- x$sub.egg1@data$AREA_SD <- rep(NA, length(x$sub.egg1@data))
  if(class(x$a) != "try-error") {
    x$sub.egg1@data$AREA <- x$a * 1e-6 # replace polygon area (use square km)
    subarea <- iceArea(x$sub.egg1@data)
    area <- subarea[[2]]
    volume <- iceVolume(x$sub.egg1@data) 
    #minlat <- iceTiming(sub.egg)
  } else {
    area <- volume <-  minlat <- minlong <- NA 
  }
  return(list(area=area, volume=volume, minlat=x$minlat, minlong=x$minlong))
} 

#calc <- trendsCalc(temp.ls)


