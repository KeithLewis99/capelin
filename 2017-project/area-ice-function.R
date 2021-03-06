################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
#  See Ale's optimization functions at the end
################################################################

# The purpose of this file is to store the funcitons required for:
#1) eggAttr(): function for extracting egg attributes from EGG_ATTR string
#2) attrTab(): extracting attribute table from raw e00 file (not 100% successful)
#3) e00Download(): download e00 data from environment Canada - see function for notes
#4) eggAttr.cols(): Turn the eggAttr data into columns that can be querried
#5) eggAttr.query(): query the egg/ice data for concentration and stage among polygons 
#6) iceSubset(): sub.egg object (SPDF) generated by the calcAreaVolume function and modified by eggAttr.cols(),and egg.Attr.query()
#7) iceTiming (): calculate the minimum latitude of ice in a given year - this represents timing of the retreat
#8) withinPolyAreaA(): query the egg/ice data for concentration and stage among polygons (could not generalize this so 4 functions instead of one)
#9) iceArea(): using the attribute table to calculate total ice area
#10) iceVolume(): using attribute table to calculate total ice volume
#11) plotIce: function for plotting ice charts 
#12) e00_to_SpatialPolygonDataframe(): e00 to avc_data (coverages) and then to SpatialPolygonsDAtaframe in sp_data
#13) loadMap(): loads maps from Rdata file
#14) subsetProject(): subset and project to WGS84 (poorly named)
#15) filterEgg(): Apply filters and return to lcc projection
#16) calcLatLong(): subset the sub.egg data by ice concentration and stage of development and calculate the minimum latitude (timing of iceretreat) using sub functions - needed to proceed with area and volume calculations
#17) calcPolyA(): calculates the area of each polygon in SpatailPolygonDataframe
#18) trendsCalc(): use iceArea() and iceVolume() to calculate the area and volume of ice for specified area
#19) calcAreaVolLat(): Main function that brings together all the other functions
#20) multiplot
#21) several lookat functions
#22) Ale's optimization functions


#) LookAt function - functions for looking at various stages of the egg data and checking that modifications are correct.


# V1 - first attempt at taking Paul's code and making functions
# v2 - add functions to subset ice
# v3 - a substantial reworking of the file - the calcAreaVol() was reworked into calcAreaVolLat() and was broken up into multiple subfunctions
## Custom functions 


############################################################
## Below here are other functions for other work
####################################################################################
##' multiplot()--------- 
#' - like mfrow
#' # Multiple plot function http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



############################################################
##' iceMedian()---------------
#'
#' @param df the output of ice-chart-processing-data (e.g. trends.m1) and data1 is output from iceSummary()
#' @return a list of:
#' 1) the subsetted data
#' 2) the median values of minlats, tice, and area
#' 3) combines 2) with the values of iceSummary()
#' @export
#'
#' @examples sub1991 <- iceMedian(trends.m1, "year < 1992", "tice < 150", iceSum.m1)
iceMedian <- function(df, subset_yr, subset_ti, data1) {
  #browser()
  #print(yr)
  #print(ti)
  df1 <- subset(df, eval(parse(text = subset_yr)) & eval(parse(text = subset_ti)))
  #df1 <- subset(df, paste("year", operator1, yr) & tice < ti)
  #df1 <- subset(df, do.call(operator1, list(get(year), yr)))
  #print(df1)
  
  df2 <- df1 %>%
    group_by(year) %>%
    summarize(dminlats = median(minlats))
  
  df3 <- df1 %>%
    group_by(year) %>%
    summarize(dtice = median(tice))
  
  df4 <- df1 %>%
    group_by(year) %>%
    summarize(darea = median(area))

  df5 <- left_join(df2, df3, by = "year")
  df5 <- left_join(df5, df4, by = "year")
  
  df6 <- subset(data1, eval(parse(text = subset_yr)))
  df7 <- left_join(df6, df5, by = "year")
  
  return(list(data = df1, meds=df5, mall=df7))
}

##################################################################
##' iceSummary()-----
#'
#' @param df the output of ice-chart-processing-data (e.g. trends.m1)
#'
#' @return a dataframe with the max of iceArea and minLat for a given year
#' @export
#'
#' @examples iceSum.m1 <- iceSummary(trends.m1)
iceSummary <- function(df){
  #browser()
  ##modify trends
  # add year
  df$year <- year(df$date)
  #convert date.y to doy
  df$jday <- yday(df$date)
  
  # calculate the maximum area for ice by year
  df1 <- df[c("date", "area", "volume", "year", "jday")]
  temp1 <- df1 %>%
    filter(jday <= 150) %>%
    group_by(year) %>%
    slice(which.max(area)) 
  temp1 

  # calculate the minimum latitude for ice by year
  df2 <- df[c("date", "minlats", "minlongs", "year", "jday")]
  temp2 <- df2 %>%
    filter(jday <= 150) %>%
    group_by(year) %>%
    slice(which.min(minlats)) 
  
  # merge summarized dataframes
  df1 <- full_join(temp1, temp2, by = "year")
  
  #convert date.y to doy
  df1$tice <- yday(df1$date.y)
  return(df1)
}

##################################################################
##' iceMedianD5()-------
#'
#' @param df the output of ice-chart-processing-data (e.g. trends.m1)
#'
#' @return the median values of minlats, tice, and area for the top 5 values
#' @export
#'
#' @examples
iceMedianD5 <- function(df, subset_yr, subset_ti, data1){
  #browser()
  df1 <- subset(df, eval(parse(text = subset_yr)) & eval(parse(text = subset_ti)))
  
  t1 <- df1 %>%
    group_by(year) %>%
    arrange(desc(area)) %>%
    slice(1:5) %>%
    summarize(d5area=median(area))

  t2 <- df1 %>%
    group_by(year) %>%
    arrange(minlats) %>%
    slice(1:5) %>%
    summarize(d5tice=median(tice))

  t3 <- df1 %>%
    group_by(year) %>%
    arrange(minlats) %>%
    slice(1:5) %>%
    summarize(d5minlats=median(minlats))
  
  df2 <- left_join(t1, t2, by = "year")
  df2 <- left_join(df2, t3, by = "year")
  
  df3 <- subset(data1, eval(parse(text = subset_yr)))
  df4 <- left_join(df3, df2, by = "year")
  return(list(data = df1, meds=df2, mall=df4))
}


##################################################################
##' subsetTestPlot()------
#'
#' @param df1 results of model 1
#' @param df2 results of model 2
#' @param date the date which is used to merge the datasets
#'
#' @return a graph of the two datasets with an abline - no values should be above the line
#'
#' @examples subsetTestPlot(m1, m2, "date")
subsetTestPlot <- function(df1, df2, date){
  mtest <- merge(df1, df2, by = eval(parse(text = date)))
  p <- ggplot(mtest, aes(area.x, area.y)) +
          geom_point() + 
          geom_abline(aes(intercept = 0, slope = 1), colour = "red", size = 2)
    #windows()
    return(p)
}

##################################################################
##' iceSummarylm()------------
#'
#' @param dataframe 
#'
#' @return - list of lm summaries
#' @export
#'
#' @examples iceSummarylm (sub2017.m1)
iceSummarylm <- function(df, med = NULL, med1 = NULL) {
  #browser()
  a <- summary(lm(paste0(med, "minlats", " ~ ", med1, "area"), data=df$mall))
  b <- summary(lm(paste0(med, "minlats", " ~ ", med, "tice"), data=df$mall))
  c <- summary(lm(paste0(med1, "area", " ~ ", med, "tice"), data=df$mall))  
  return(list(minlats_v_area = a, minlats_v_tice = b, area_v_tice = c))
}

##################################################################
##' iceDateScatter()--------
#'
#' @param df 
#'
#' @return 1 graphs of date v minlats, and 2 faceted histograms of minlats by date
#' @export
#'
#' @examples
#' 
iceDateScatter <- function(df, d = NULL, vline = NULL){
  #browser()
p1 <- ggplot(data = df$data, aes(x = date, y = minlats)) + 
  geom_point() + 
  geom_smooth(method=lm)
  
p2 <- ggplot() + 
  geom_violin(data = sub1991.m1$data, aes(x = year, y = minlats, group = year)) +
  geom_point(data = sub1991.m1$mall, aes(x=year, y = dminlats), colour = "red")

  
p3 <- ggplot() + 
  geom_violin(data = sub1991.m1$data, aes(x = year, y = area, group = year)) +
  geom_point(data = sub1991.m1$mall, aes(x = year, y = darea), colour = "red")

return(list(p1=p1, p2=p2, p3=p3))
}

##################################################################
##' iceYearBox()----------
#'
#' @param df1 
#' @param df2 
#' @param df3 
#' @param title of the plots
#'
#' @return
#' @export
#'
#' @examples
iceYearBox <- function(df1, df2, df3, title) {
  p1 <- ggplot(data = df1, aes(x = year, y = area, group = year)) + 
    geom_boxplot()
  p2 <- ggplot(data = df1, aes(x = year, y = minlats, group = year)) + 
    geom_boxplot()
  p3 <- ggplot(data = df2, aes(x = year, y = area, group = year)) + 
    geom_boxplot()
  p4 <- ggplot(data = df2, aes(x = year, y = minlats, group = year)) + 
    geom_boxplot()
  p5 <- ggplot(data = df3, aes(x = year, y = area, group = year)) + 
    geom_boxplot()
  p6 <- ggplot(data = df3, aes(x = year, y = minlats, group = year)) + 
    geom_boxplot()
  #windows()
  #return(list(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6))
  cowplot::plot_grid(p1, p2, p3, p4, p5, p6, labels = c("m1", "", "m2", "", "m3", ""), ncol=2)
  #multiplot(p1, p3, p5, p2, p4, p6, cols=2)
}

##################################################################
##' iceScatterSummary()--------
#'
#' @param df 
#'
#' @return 1 graphs of date v minlats, and 2 faceted histograms of minlats by date
#' @export
#'
#' @examples
#' 
iceScatterSummary <- function(df, x = NULL, x1 = NULL, y = NULL, y1 = NULL) {
  #browser()
  p1  <- ggplot(data = df, aes_string(x = x, y = y)) + geom_point() + geom_smooth(method=lm)
  s1 <- summary(lm(paste0(y, "~", x), data=df))
  # plot tice against ice area
  p2 <- ggplot(data = df, aes_string(x = x1, y = y)) + geom_point() + geom_smooth(method=lm)
  s2 <- summary(lm(paste0(y, "~", x1), data=df))
  #plot ice area against minlats
  p3 <- ggplot(data = df, aes_string(x = x, y = y1)) + geom_point() + geom_smooth(method=lm)
  s3 <- summary(lm(paste0(y1, "~", x), data=df))
  return(list(minlats_v_area = s1, minlats_v_tice = s2, area_v_tice = s3, p1=p1, p2=p2, p3=p3))
  multiplot(p1, p3, p2, cols=2)
}

##################################################################
##' iceOuput()------
#'
#' @param df a list with areas, volumes, minlats, minlongs
#'
#' @return the dataframe from the list df
#' @export
#'
#' @examples iceOutput(df)
iceOuput <- function(df, date) {
  #browser()
  df1 <- rbind(data.frame(date = date, 
                                area = df$areas, 
                                volume = df$volumes,
                                minlats = df$minlats,
                                minlongs = df$minlongs)) 
  df1 <- df1[order(df1$date), ]  
  return(df1)
}

##################################################################
##' lmoutputSummary()------
#'
#' @param ls - a list of summary(lm(formula)) returned from iceScatterSummary() for the comparison of area v tice, minlats v tice, and area v minlats
#'
#' @return printed values of r-squared values for each formula
#' @export
#'
#' @examples lmoutputSummary(ls)
lmoutputSummary <- function(ls) {
#browser()
  a <- "tice_v_area"
  b <- "tice_v_minlats"
  c <- "minlats_v_area"
  comp <- rbind(a, b, c)
  rsq_val <- rbind(ls$minlats_v_area$r.squared, ls$minlats_v_tice$r.squared, ls$area_v_tice$r.squared)
  slope <- rbind(ls$minlats_v_area$coefficients[2,1], ls$minlats_v_tice$coefficients[2,1], ls$area_v_tice$coefficients[2,1])
  out <- cbind(comp = comp, rsq_val = round(rsq_val, 3), slope = round(rsq_val, 2))
  colnames(out) <- c("measure", "rsq", "slope")
  return(out)
}

##################################################################
##' area_ticeScatter()-----
#'
#' @param df1 max value
#' @param df2 median value
#' @param df3 median of top 5 values
#'
#' @return graphs of max/median/d5median v. tice
#' @export
#'
#' @examples atScat.m1 <- area_ticeScatter(iceSum.m1, iceMed.m1, iceMedD5.m1, "tice", "area", "m1")
#' ggsave("figs/41-atScat.m1.pdf")    
area_ticeScatter <- function(df1, df2, df3, x, y1, y2, y3, sub_val){
     p1 <- ggplot(data = df1, aes_string(x=x, y=y1)) + 
          geom_point() +         
          geom_smooth(method=lm)  +
          ggtitle(paste(sub_val))
     p2 <- ggplot(data = df2$mall, aes_string(x=x, y=y2)) +  
          geom_point() +         
          geom_smooth(method=lm)
     p3 <- ggplot(data = df3$mall, aes_string(x=x, y=y3)) +  
          geom_point() +         
          geom_smooth(method=lm)
    out_plot <- cowplot::plot_grid(p1, p2, p3, labels = c("max", "med", "d5med"), ncol=2)
}

##################################################################
##' area_ticeRsq()------
#'
#' @param df1 max value
#' @param df2 median value
#' @param df3 median of top 5 values
#' @param x values to be pasted in as the x value in the lm
#' @param y values to be pasted in as the y value in the lm
#'
#' @return
#' @export
#'
#' @examples area_ticeRsq(iceSum.m1, iceMed.m1, iceMedD5.m1, x = "tice", y = "area")
area_ticeRsq <- function(df1, df2, df3, x, y) {
     #browser()
     s1 <- summary(lm(paste0(y, "~", x), data=df1))
     # plot tice against ice area
     s2 <- summary(lm(paste0("d", y, "~", "d", x), data=df2$mall))
     #plot ice area against minlats
     s3 <- summary(lm(paste0("d5", y, "~", "d5", x), data=df3$mall))
     a <- "Max"
     b <- "Median"
     c <- "d5Median"
     comp <- rbind(a, b, c)
     rsq_val <- rbind(s1$r.squared, s2$r.squared, s3$r.squared)
     slope <- rbind(s1$minlats_v_area$coefficients[2,1], s2$minlats_v_tice$coefficients[2,1], s3$area_v_tice$coefficients[2,1])
     out <- cbind(comp = comp, rsq_val = round(rsq_val, 3), slope = round(rsq_val, 2))
     colnames(out) <- c("measure", "rsq", "slope")
     return(out)
}
