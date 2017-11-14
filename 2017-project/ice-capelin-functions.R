################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
#  See Ale's optimization functions at the end
################################################################

##################################################################################################### OPtimization funcitons
##  SSQCapelinDome()--------------
### Feb 3 2015 - ADB
### Code to fit the ice capelin model presented in Buren et al 2014
### this is not elegant code, but at the time I ran the nalyses was good enough
### I define 2 functions: one is the objective function SSQCapelinDome   and the second function
### obtains the Expected Log Capelin BIomass (CapelinDome)

### The data must be stored in a dataframe called capelin with at least 3 columns: 'year','tice','logcapelin'
### logcapelin is the observed capelin biomass (in log scale). The data frame capelin cancontain more columns
### NA values of logcapelin (i.e. years when there were no capelin surveys) are not used during optimization

## Objective function - part of optimization function

#' SSQCapelinDome()------
#' most params passed from calcFit_all():
#' @param params - from par: a series of values to help the optim function 
#' @param dataf - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param var - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#'
#' @return - the expected value of log capelin biomass
#' @export
#'
#' @examples CapelinDomeFit <- optim(par = par,dataf = df[which(df$logcapelin!='NA'), c('year', paste(var),'logcapelin')], form1=form1, form2=form2, var=var, fn = SSQCapelinDome, method=c("BFGS"))
#' 
SSQCapelinDome <- function(params, dataf, form1, form2, var){
     #browser()
     dataf <- as.data.frame(dataf) # needed because optim doesn't work with tibbles.....grrrrrrrr!!!!!
     Alpha <- params[1]
     Beta <- params[2]
     Gamma <- params[3]
     year <- dataf[,1]
     #tice <- dataf[,2]
     tmp <- dataf[,2] #variable of interest in form1 and form2
     assign(var, tmp)
     logcap <- dataf[,3]
     # this is based on MSY
     ELogCapBiom <- ifelse(year<1991, eval(parse(text = form1)), eval(parse(text = form2)))
     sum((logcap-ELogCapBiom)^2)
}

# Alpha*tice*(1-(tice/Beta)), Alpha*tice*(1-(tice/Beta))*Gamma

## Function to obtain Expected Log Capelin Biomass    
#' most params passed from calcFit_all():
#' @param params - from par: a series of values to help the optim function 
#' @param dataf - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param var - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' 
#' @return - #### check this
#' @export
#'
#' @examples df$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par), dataf = df[,c('year', paste(var))], form1, form2, var)
#' 
CapelinDome <- function(params,dataf, form1, form2, var){
     #browser()
     dataf <- as.data.frame(dataf)
  Alpha <- params[1]
  Beta <- params[2]
  Gamma <- params[3]
  year <- dataf[,1]
  tmp <- dataf[,2]
  assign(var, tmp)
  ELogCapBiom <- ifelse(year<1991, eval(parse(text = form1)), eval(parse(text = form2)))
  ELogCapBiom
}

#####################################################################################################
##' modelGraphs()-------------
#' figure for plotting predicted values of models over a range of years
#' @param x = years
#'
#' @return graphs from 1970-2016 (or based on vector)
#' @export
#'
#' @examples modelGraphs(years)
modelGraphs <- function(years){
  for(i in seq_along(years)) {
    
    ## subset data
    mod.dat <- trends[year == years[i]]
    mod.dat$t <- mod.dat$doy
    mod.dat$a <- mod.dat$area/1000
    mod.p <- ggplot(mod.dat) + 
      geom_point(aes(x = t, y = a)) + 
      theme_bw() + 
      xlab('Day of Year (NOV - JUL)') +
      ylab('Area') +
      theme(panel.grid.major = element_line(linetype = "dotted")) # make base layer and plot data
    
    
    ## extract basic parameters
    raw.amax[i] <- max(mod.dat$a)
    raw.tmax[i] <- mod.dat$t[which.max(mod.dat$a)]
    
    ## approx start vals
    j <- which(mod.dat$a > quantile(mod.dat$a, prob = 0.50)) # top values
    h <- max(mod.dat$a[j]) / sqrt(2 * pi)
    s <- max(mod.dat$t[j]) - min(mod.dat$t[j])
    tm <- mean(mod.dat$t[j])
    # b <- h/tm
    
    ## run models
    mod1 <- try(nls(a ~ (h / sqrt(2 * pi)) * (exp(-((t - tm) ^ 2) / (2 * s ^ 2))), 
                    data = mod.dat, 
                    start = list(h = h, s = s, tm = tm),
                    lower = c(0, 0, -60), # except for tm (let search back to Nov 1 [doy -60]), constrain to positive parameter space 
                    control = nls.control(maxiter = 100), 
                    algorithm = "port"))
    mod2 <- try(nls(a ~ (h / sqrt(2 * pi)) * (exp(-((t - tm) ^ 2) / (2 * ifelse(t < tm, s1, s2) ^ 2))),  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2993707/   
                    data = mod.dat, 
                    start = list(h = h, s1 = s/2, s2 = s/2, tm = tm),
                    lower = c(0, 0, 0, -60), # except for tm (let search back to Nov 1 [doy -60]), constrain to positive parameter space
                    control = nls.control(maxiter = 100),
                    algorithm = "port"))
    
    if(class(mod1) != "try-error" | class(mod2) != "try-error") {
      
      ## predict
      preds <- data.frame(t = min(mod.dat$t):max(mod.dat$t))
      if(class(mod1) != "try-error") preds$mod1 <- predict(mod1, newdata = preds)
      if(class(mod2) != "try-error") preds$mod2 <- predict(mod2, newdata = preds)
      preds <- data.table(preds)
      preds <- melt(preds, id.var = "t", variable.name = "mod", value.name = "a")
      
      if(class(mod2) != "try-error") {
        amax[i] <- coefficients(mod2)["h"] / sqrt(2 * pi)
        tmax[i] <- coefficients(mod2)["tm"]
        tstart[i] <- coefficients(mod2)["tm"] - sqrt(-2 * log(0.20)) * coefficients(mod2)["s1"] # start of advance (20% of amax)
        tend[i] <- coefficients(mod2)["tm"] + sqrt(-2 * log(0.20)) * coefficients(mod2)["s2"] # end of retreat (20% of amax)
        adv.mag[i] <- coefficients(mod2)["h"] * coefficients(mod2)["s1"] / 2
        ret.mag[i] <- coefficients(mod2)["h"] * coefficients(mod2)["s2"] / 2
        #round(mod.amax[i] * 0.20, 3) == round(predict(mod2, newdata = data.frame(t = ret.end[i])), 3) # check math
        ## calculate advance and retreat maginitude
        
        x <- c(tmax[i], tstart[i], tend[i])
        y <- predict(mod2, newdata = data.frame(t = x))
        seg.df <- data.frame(x1 = x, x2 = x, y1 = rep(0, length(y)), y2 = y)
        mod.p <- mod.p + 
          #geom_ribbon(aes(x = t, ymin = 0, ymax = a), 
          #                           data = preds[mod == "mod2"],
          #                           fill = "#00BFC4", alpha = 0.15) +
          geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                       data = seg.df, linetype = "dashed", col = "#00BFC4") +
          geom_hline(yintercept = amax[i], linetype = "dashed", col = "#00BFC4") # plots intersection - not sure of what yet
      }
      
      ## add prediction to plot 
      mod.p <- mod.p + 
        geom_line(aes(x = t, y = a, col = mod), preds, size = 0.7) # plots the predicted values 
      #labs(mod = "Model") #- tried to change the legend
      
      #scale_fill_continuous(guide = guide_legend(title = "Model"))  
    }
    
    ## print plot and clean-up
    print(mod.p + ggtitle(years[i])) # add the years
    rm(mod1, mod2)
    
  }
}


############################################################
`+.uneval` <- function(a,b) {
     `class<-`(modifyList(a,b), "uneval")
}
##################################################################
##' optimGraphs()------
#'
#' @param df 
#' @param reg1 
#' @param reg2 
#' @param yearInt 
#' @param lnbiomassInt 
#' @param title
#' @param var
#'
#' @return
#' @export
#'
#' @examples
optimGraphs <- function(df, reg1, reg2, yearInt, lnbiomassInt, title, var){
  #browser()
  p1 <- ggplot(df, aes(x = year, y = logcapelin)) + 
    geom_errorbar(width = 0.3, colour = "black", aes(ymin=logcapelinlb, ymax=logcapelinub)) + 
    geom_point(shape=16, size=3)  +
    geom_line(aes(y=ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
    geom_line(aes(y=ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
    scale_y_continuous(limits = c(0,10), breaks = lnbiomassInt) +
    scale_x_continuous(limits = c(1982,2018), breaks = yearInt) +
    xlab('Year') +
    ylab('ln (Capelin biomass (ktons))') + 
    theme_bw()
  
  p2 <- ggplot(df, aes(x = year, y = capelin)) +
    geom_errorbar(width = 0.3, colour = "black", aes(ymin=capelinlb, ymax=capelinub)) + 
    geom_point(shape=16, size=3)  +
    geom_line(aes(y=exp(ExpectedLogBiomass)), colour="red", linetype=1, size=1.25) +
    geom_line(aes(y=exp(ExpectedLogBiomassOld)), colour="blue", linetype=1, size=1.25) +
    #scale_y_continuous(limits = c(0,8500), breaks = biomassInt) +
    scale_x_continuous(limits = c(1982,2018), breaks = yearInt) +
    xlab('Year') +
    ylab('Capelin biomass (ktons)') + 
    theme_bw() +
    annotate("text", x = 2008, y = 8400, label = "Model estimates to 2014") + # all following for the legend
    annotate("text", x = 2008, y = 8000, label = "Model estimates to 2010") +
    annotate("segment", x = 2001, xend = 2003, y = 8400, yend = 8400, colour = "red") +
    annotate("segment", x = 2001, xend = 2003, y = 8000, yend = 8000, colour = "blue")
  
  
  p3 <- ggplot() +
    geom_line(data = reg1, aes_string(x = var) + aes(y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
       geom_line(data = reg2, aes_string(x = var) + aes(y = ExpectedLogBiomass), colour="green", linetype=1, size=1.25) +
       geom_line(data = reg1, aes_string(x = var) + aes(y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
       geom_line(data = reg2, aes_string(x = var) + aes(y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
       geom_point(data = subset(df, year < 1991), aes_string(x = var) + aes(y = logcapelin), shape=2, size=3) +
       geom_point(data = subset(df, year > 1991), aes_string(x = var) + aes(y = logcapelin), shape=15, size=3) + 
       geom_errorbar(data = subset(df, year < 1991),  aes_string(x = var) + aes(ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
       geom_errorbar(data = subset(df, year > 1991), aes_string(x = var) + aes(ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
       xlab(paste(var)) +
       ylab("ln (Capelin biomass (ktons))") + 
       #ylim(0,9) +
       theme_bw()
  
  # make multiplot
  #pdf('ice-capelin-update-2014-new.pdf',height=9,width=9,pointsize =8)
  cowplot::plot_grid(p1, p2, p3, labels = c(paste(title)), ncol=2)
  #multiplot(p1, p3, p2, cols=2)
  #dev.off()
  
}
# make optimization graphs by year and in comparison to ice

################
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


###############################################################
##' loadSubsetDatasets()-----
#'
#' @param df - the capelin dataset
#' @param pat - a pattern of file names
#' @param N - the number of files
#'
#' @return a list of dataframes with the summaries of year, area, minlats, tice, and capelin info for each subset (m1-m6)
#' @export
#'
#' @examples cape <- loadSubsetDatasets(capelin_join, pattern, 6)
loadSubsetDatasets <- function(df, pat, N){
  #browser()
  name <- paste0("capelin_m", 1:N) # create names of dataframes
  x <- vector("list", N) # create a list to fill
  for (i in seq_along(pat)) {
    print(pat[i])   # print to see what has been loaded
    cap <- read_csv(paste0("output-processing/", pat[i])) 
    cap <- left_join(cap, df, by = "year")
    cap <- rename(cap, maxarea = area)
    cap <- rename(cap, minlat = minlats)
    x[[i]] <- cap # put cap in list
    #x[[name]] <- x
  }
  names(x) <- paste(name) #name elements of list
  return(x)
}     


##########################################################################
##' calcFit()----
#' calcFit performs the optim function for a single dataframe or list
#' most param are passed from calcFit_all()
#' @param df - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param var - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param par - a series of values to help the optim function
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param x_range - range of possible valus on the x-axis
#'
#' @return a list with the new expected values of the caplein, the CapelinDomeFit, regime 1 and regime 2
#' @export
#'
#' @examples AleMaxArea <- calcFit(capelin_ale, var = "maxarea", par = c(1, 500, 0.6), form1 = "Alpha*tmp*(1-(tmp/Beta))", form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma", x_range = c(0:650))

#' note that this returns a warning "Unknown or uninitialised column: 'par'." which apparently is a tibble problem!!
#' calcFit_all <- function(ls, titlenames, var, par, form1 = NULL, form2 = NULL, x_range)
calcFit <- function(df, var, par, form1 = NULL, form2 = NULL, x_range) {
  #browser()
     #print(environment())
     CapelinDomeFit <- optim(par = par,
                          dataf = df[which(df$logcapelin!='NA'), c('year', paste(var),'logcapelin')], form1=form1, form2=form2, var=var,
                           fn = SSQCapelinDome, 
                          method=c("BFGS"))
  
  # before 2010
  CapelinDomeFitOld <- optim(par=par,
                         dataf = df[which(df$year < 2011 & df$logcapelin!='NA'),
                         c('year', paste(var), 'logcapelin')],
                         form1=form1, form2=form2, var=var,
                         fn = SSQCapelinDome, method=c("BFGS"))
  
  ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
  df$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par), dataf = df[,c('year', paste(var))], form1, form2, var)
  
  df$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par), dataf = df[,c('year', paste(var))], form1, form2, var)
  
  # attach the optimization curves of capelin abundance to ice data
  xtice <- expand.grid(year = c(1990,2000), var = as.numeric(paste(x_range)))
  colnames(xtice)[2] <- c(paste(var))
  xtice <- xtice[order(xtice[2]),]
  xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par),dataf = xtice, form1, form2, var)
  xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par),dataf = xtice, form1, form2, var)
  #not sure what these are for but used in plots below but creates a data set where all values of year are the same????
  regime1 <- xtice[which(xtice$year == 1990),]
  regime2 <- xtice[which(xtice$year == 2000),]
return(list(df = df, cdf = CapelinDomeFit, regime1 = regime1, regime2 = regime2))  
}


###############################################################
##' loadSubsetDatasets()-----
#'
#' @param df - the capelin dataset
#' @param name - name of output list
#' @param pat - a pattern of file names
#' @param N - the number of files
#' @param var1 - area, darea, d5area
#' @param var2 - minlats, dminlats, d5minlats
#' @param nvar1 - new name of var1 
#' @param nvar2 - new name of var2
#'
#' @return a list of dataframes with the summaries of year, area, minlats, tice, and capelin info for each subset (m1-m6)
#' @export
#'
#' @examples cape <- loadSubsetDatasets(capelin_join, pattern, 6, area, maxarea, minlats, minlat)
loadSubsetDatasets1 <- function(df, name, pat, N, var1 = NULL, var2 = NULL, nvar1 = NULL, nvar2 = NULL){
     #browser()
     #var1 <- enquo(var1)
     #nvar1 <- enquo(nvar1)
     name <- paste0(name, 1:N) # create names of dataframes
     x <- vector("list", N) # create a list to fill
     for (i in seq_along(pat)) {
          print(pat[i])   # print to see what has been loaded
          cap <- read_csv(paste0("output-processing/", pat[i])) 
          cap <- left_join(cap, df, by = "year")
          #cap <- rename(.dots = setNames(var1, nvar1))
          cap <- rename(cap, !!nvar1 := !!var1)  # !! unquotes the input but the code is not valid.  The ":=" makes it valid.....somehow
          cap <- rename(cap, !!nvar2 := !!var2) 
          x[[i]] <- cap # put cap in list
          #x[[name]] <- x
     }
     names(x) <- paste(name) #name elements of list
     return(x)
}     

###############################################################

# prints plots from a list in sequence
testPlot <- function(ls1, titlenames){
     #browser()
     for(i in 1:length(ls1)){
          df1 <- as.data.frame(ls1[[i]])     
          p1 <- ggplot(data = df1,
                  aes(x = max_area, y = capelin)) +
                    geom_point() + 
                    ggtitle(paste(titlenames[i]))
     print(p1)
     }
}

###############################################
##' capelinAreaPlot()-----
#'
#' @param ls1 list of max area
#' @param ls2 list of med area
#' @param ls3 list of median of top 5 values
#' @param i the subset m1-m6
#' @param titlenames m1-m6
#'
#' @return cowplot of capelin v area
#' @export
#'
#' @examples capelinAreaPlot(cape, med_cape_all, d5med_cape_all, 1, titlenames)
capelinAreaPlot <- function(ls1, ls2, ls3, i, titlenames) {
     df1 <- as.data.frame(ls1[[i]])
     df2 <- as.data.frame(ls2[[i]])
     df3 <- as.data.frame(ls3[[i]])
     p1 <- ggplot(data = df1,
                  aes(x = max_area, y = capelin)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     p2 <- ggplot(data = df2,
                  aes(x = med_area, y = capelin)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     p3 <- ggplot(data = df3,
                  aes(x = d5med_area, y = capelin)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     
     cowplot::plot_grid(p1, p2, p3, labels = c("Max Area", "Median Area", "D5Median Area"), nrow=1)    
}

###############################################
##' lnCapelinArea()----
#'
#' @param ls1 - list (max_area, med, d5med)
#' @param ls2 - as above but <= 1991
#' @param ls3 - as above but > 1991
#' @param xaxis - column for xaxis
#' @param yaxis - column for yaxis
#' @param i - subset
#' @param titlenames - names of subset 
#'
#' @return - graphs
#' @export
#'
#' @examples testAnotherPlot(cape, cape_1991, cape_2017, "max_area", "logcapelin", 1, titlenames)

lnCapelinArea <- function(ls1, ls2, ls3, xaxis, yaxis, i, titlenames){
     #browser()
     df1 <- as.data.frame(ls1[[i]])
     df2 <- as.data.frame(ls2[[i]])
     df3 <- as.data.frame(ls3[[i]])
     p1 <- ggplot(data = df1,
                  aes_string(x = xaxis, y = yaxis)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     p2 <- ggplot(data = df2,
                  aes_string(x = xaxis, y = yaxis)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     p3 <- ggplot(data = df3,
                  aes_string(x = xaxis, y = yaxis)) +
          geom_point() + 
          ggtitle(paste(titlenames[i])) +
          theme(plot.title = element_text(hjust = 0.5))
     cowplot::plot_grid(p1, p2, p3, labels = c("All years", "<1992", ">=1992"), nrow=1)
}
###############################################

#' calcFit_all()----
#' #calcFit performs the optim function for a single dataframe.  CalcFit_all is a wrapper that performs the optim function on a list of dataframes for the different ice subsets
#'
#' @param ls a list of dataframes for the different ice subsets (e.g. m1-m6) with associated capelin data
#' @param titlenames for the different ice subsets (e.g. m1-m6)
#' @param var - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param par - a series of values to help the optim function
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param x_range - range of possible valus on the x-axis
#'
#' @return - a list of lists.  A list for each of the ice subsets (i.e., list of 6). Each list is a list of 4 with a dataframe (data), cdf (list of 5 which is output of the optim funtion), and 2 dataframes (regime1 and regime2) that are used to produce the associated graphs  
#' @export
#'
#' @examples MaxTice <- calcFit_all(cape, titlenames, par = c(1, 200, 0.6), var = "tice", form1 = "Alpha*tmp*(1-(tmp/Beta))", form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma", x_range = c(0:190,173.515,187.768))
#' 
calcFit_all <- function(ls, titlenames, var, par, form1 = NULL, form2 = NULL, x_range) {
     #browser()
     optim_ls <- rep(list(list()), length(titlenames))
     names(optim_ls) <- titlenames
     for(i in 1:length(ls)){
          df1 <- as.data.frame(ls[[i]])
          optim <- calcFit(df1, var=var, par=par, form1 = form1, form2 = form2, x_range = x_range)
          optim_ls[[i]] <- optim
     }
     return(list(optim_ls = optim_ls))
}     
