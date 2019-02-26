## load libraries
require(GISTools)
require(maps)
library(raster)
library(rgdal)
library(spatial)
require(sfsmisc)
require(extrafont)
## import and project map data

setwd("map")

na <- readOGR("north_america/bound_p.shp", layer="bound_p")
ecan <- na[na$NAME=="Newfoundland and Labrador / Terre-Neuve-et-Labrador" |
          na$NAME=="Nova Scotia / Nouvelle-Écosse" |
          na$NAME=="Prince Edward Island / Île-du-Prince-Édouard" |
          na$NAME=="Saint-Pierre et Miquelon" |
          na$NAME=="New Brunswick / Nouveau-Brunswick" |
          na$NAME=="Quebec / QuÃ©bec" |
          na$NAME=="Nunavut" |
          na$NAME=="Maine" |
          na$NAME=="New Hampshire" |
          na$NAME=="Massachusetts" |
          na$NAME=="Kalaallit Nunaat",]
ecan.nad <- spTransform(ecan, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")) # projection with long lat
nafo<-  readOGR("NAFO_SHP/NAFO_Divns.shp", layer="NAFO_Divns")
nafo.nad <- spTransform(nafo, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")) # projection with long lat
bathy.proj.nad <- raster("bathy/nl_nad.grd", values=TRUE)
bathy.gull.nad <- crop(bathy.proj.nad, extent(-54, -50, 46, 49)) # projection with long lat



  latitudelabel<-expression(paste('Latitude (',degree,'N)'))
  longitudelabel<-expression(paste('Longitude (',degree,'W)'))

## map

 loadfonts(device="postscript")
cairo_ps(filename = "map_v1.ps", 
         width = 8.5, height = 8.5, pointsize = 12,
         onefile = FALSE, family = "Arial", bg = "white",
         antialias = c("default"))
par(mar=c(3,4,2,2))
plot(nafo.nad[c(37,31,3:10,13:14,16:19,20,22,27:30,33,34,35,36,39:49,51,53,57,58,60:69,71,72,74,75),],  xlim=c(-60,-47), ylim=c(46,56), bg="white")
  plot(ecan.nad, col="grey",  xlim=c(-60,-47), ylim=c(46,56), bg="white",add=T)
contour(bathy.proj.nad, levels=c(-100,-200,-400,-1000,-2000), add=TRUE, col="dark grey")
degAxis(1, seq(-60,-44,6), mgp=c(3,0.5,0))
  mtext(longitudelabel,1,line=2,cex=1)
degAxis(2, seq(46,56,5))  
  mtext(latitudelabel,2,line=2,cex=1)
#   lines(c(-55.65,-52.50),c(53.23,55.07),lwd=5,col='red')
#   lines(c(-52.97,-49),c(48.73,50),lwd=5,col='red')
#   text(-55,49.9,"Notre  Dame \nBay",srt=-19,cex=0.9)
#   text(-52.8,47.8,"CB",srt=55,cex=1.2)
#   text(-53.5,48,"Trinity Bay",srt=55,cex=0.9)
#   text(-52.3,47.6,"SJ",cex=1.2)
   text(-48,53.5,"2J",cex=1.5)
   text(-48,50.5,"3K",cex=1.5)
   text(-48,47,"3L",cex=1.5)
 #  text(-54.3,54.3,"Seal Island Line",srt=43,cex=0.9)
#   text(-51,49.7,"Bonavista Line",srt=28,cex=0.9)
   
   x <- c(-53.89560, -48.76948, -48.68590, -48.37945,
          -47.51581, -47.23721, -47.23721, -47.07006,
          -47.73868, -52.80909, -52.64193, -52.39120,
          -52.86481, -53.28270, -54.00704, -53.89560)     
   y <- c(49.93753, 49.93753, 49.34143, 48.76286, 48.23688,
          47.69338, 46.95701, 45.79987, 45.62455, 45.64208,
          47.13234, 47.79857, 48.58753, 49.28883, 49.63948, 49.93753)
   x <- c(-54.02727, -50.09795, -49.76211, -48.58668, -47.20974, -47.17615, -47.07540, -52.75108, -52.58316, -52.34807, -53.22125, -53.92652, -53.96010)
   y <- c(49.92000, 49.94114, 48.79984, 48.20806, 47.82763, 46.81315, 45.69299, 45.67186, 47.10904, 47.74309, 49.28595, 49.64524, 49.92000)
   #p1 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x, y))), ID = "LM")))
   p1 <- SpatialLines(list(Lines(list(Line(cbind(x, y))), ID = "LM")))
   proj4string(p1) <- proj4string(nafo.nad)
   stratum <- gIntersection(nafo.nad, p1, byid = TRUE) # intesection of water (ocean surronding NL) and Lake Melville
   stratum <- gBuffer(stratum, width = 0.1) # add a small buffer
   lines(p1, col = "red", lwd = 2) # plots L
   
   points(  -53.7818,	47.6358,pch=16,cex=1,col='black')
 #  points(-53.181375964,47.683192975,pch=16,cex=1,col='blue')
#   points(-52.712017,47.5612626,pch=16,col='black',cex=1)  
 #  text(-52.5,47.7,"St John's",srt=0,cex=0.9,adj=0)   
     points(-52.59, 47.55, pch=16, col='black', cex=1)
     text(-52.4, 47.5, "Station 27",srt=0,cex=0.9,adj=0)

     
#arrows(x0=-53.181375964,y0=47.683192975,x1=-52.55,y1=48, length = 0.1, angle = 20, code = 3,col='black',lwd=1.5)
#arrow.plot(a1=matrix(-53.181375964,47.683192975), a2=matrix(-52.55,48))
 #    points(-52.583611,47.534722,pch=18,col='black',cex=1.5)   # stn 27
#lines(x=c(-53,-52.55),y=c(47.75,48),lwd=2)
#plot(arrow(angle = 90, length = unit(0.25, "inches"),ends = "last", type = "closed"))
#p.arrows(x1=-53,x2=-52.6,y1=47.75,y2=48, fill = "black",size=0.5,lwd=1.85)
#text(-52.5, 48, "Bryant's Cove",cex=0.9,adj=0)
text(-56,46.2, "Bellevue Beach",cex=0.9,adj=0) 
p.arrows(x1=-53.7818,x2=-53.7818,y1=47.55,y2=46.4, fill = "black",size=0.5,lwd=1.85)
#map.scale(x=-52, y=56,relwidth=0.1)
#map.scale(x=-52, y=56, ratio=FALSE, relwidth=0.1)
#north.arrow(xb=-46.5, yb=56, len=0.1, lab="N",col="black")
box()
dev.off()


