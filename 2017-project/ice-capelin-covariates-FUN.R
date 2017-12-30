#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-08, R version 3.3.3 (2017-03-06)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-11.08 #

# The purpose of this file is to:
# 1) store functions in support of the file ice-capelin-covariates.  Have updated some to account for multiple covariates
     # a) load the ice subsets
     # b) perform the optimization 
     # c) graph the outputs of the optimization.


#' loadSubsetDatasets(deprecated)
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

#' loadSubsetDatasets1
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


#' optimGraphs1_all(deprecated)
#'
#' @param ls a list of ice subsets with capelin and pseudocalanus data
#' @param var the first variable of interest (usually tice)
#' @param var2 the second variable of interest
#' @param var2val the value of the second variable - median value
#' @param file_name an evocative name signifying the type of analysis
#' @param saveGraph "yes" or "no" - "no" allows you to view graphs
#'
#' @return a list of output from function optim (parameter estimates, likelihood values, and convergence)
#' @export
#'
#' @examples optimGraphs2_all(MaxTice7, "tice", "Ssurface_tows_lag2", var2val = 280, "opt_div3c", "yes")
optimGraphs1_all <- function(ls, var, var2=NULL, file_name){
     for(i in 1:length(ls$optim_ls)){
          df1 <- as.data.frame(ls$optim_ls[[i]]$df)
          df2 <- as.data.frame(ls$optim_ls[[i]]$regime1)
          df3 <- as.data.frame(ls$optim_ls[[i]]$regime2)
          mm <- optimGraphs2(df1, df2, df3, yearLim, yearInt, lnbiomassInt,  titlenames[i], var, var2)
          ggsave(mm, filename = paste0("figs/covariates/", file_name, titlenames[i], ".pdf"), width=10, height=8, units="in")
     }
}


#' optimGraphs1 (deprecated)
#'
#' @param df the dataframe provided by the list
#' @param reg1 a set of expected values for pre1990 - deprecated
#' @param reg2 a set of expected values for post 2002
#' @param yearLim the year limits
#' @param yearInt the year intervals
#' @param lnbiomassInt the ln(Biomass) intervats
#' @param title - names of the output graphs
#' @param var - the first variable of interest (usually tice)
#' @param var2 - the second variable of interest
#'
#' @return a set of three graphs for each subset
#' @export
#'
#' @examples nested in optimGraphs1_all(MaxTice7, "tice", "Ssurface_tows_lag2", "opt_div3b")
optimGraphs1 <- function(df, reg1, reg2, yearLim, yearInt, lnbiomassInt, title, var, var2){
     #browser()
p1 <- ggplot(df, aes(x = year, y = logcapelin)) + 
     geom_errorbar(width = 0.3, colour = "black", aes(ymin=logcapelinlb, ymax=logcapelinub)) + 
     geom_point(shape=16, size=3)  +
     geom_line(aes(y=ExpectedLogBiomass), colour="blue", linetype=1, size=1.25) +
     #geom_line(aes(y=ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
     scale_y_continuous(limits = c(-4,10), breaks = lnbiomassInt) +
     scale_x_continuous(limits = yearLim, breaks = yearInt) +
     xlab('Year') +
     ylab('ln (Capelin biomass (ktons))') + 
     theme_bw()
     
p2 <- ggplot(df, aes(x = year, y = capelin)) +
     geom_errorbar(width = 0.3, colour = "black", aes(ymin=capelinlb, ymax=capelinub)) + 
     geom_point(shape=16, size=3)  +
     geom_line(aes(y=exp(ExpectedLogBiomass)), colour="blue", linetype=1, size=1.25) +
     #geom_line(aes(y=exp(ExpectedLogBiomassOld)), colour="blue", linetype=1, size=1.25) +
          #scale_y_continuous(limits = c(0,8500), breaks = biomassInt) +
     scale_x_continuous(limits = yearLim, breaks = yearInt) +
     xlab('Year') +
     ylab('Capelin biomass (ktons)') + 
     theme_bw() 
#annotate("text", x = 2008, y = 8400, label = "Model estimates to 2014") + # all following for the legend
     #annotate("text", x = 2008, y = 8000, label = "Model estimates to 2010") +
     #annotate("segment", x = 2001, xend = 2003, y = 8400, yend = 8400, colour = "red") +
     #annotate("segment", x = 2001, xend = 2003, y = 8000, yend = 8000, colour = "blue")
     
#browser()
p3 <- ggplot() + 
     geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomass"), colour="red", linetype=1, size=1.25) + 
     #geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomassOld"), colour="blue", linetype=1, size=1.25) +
     geom_line(data = subset(reg2, eval(parse(text=var2))==230), aes_string(x = var, y = "ExpectedLogBiomass"), colour="green", linetype=1, size=1.25) +
     geom_point(data = subset(df, year > 1991), aes_string(x = var, y = "logcapelin"), shape=15, size=3) + 
     geom_errorbar(data = subset(df, year > 1991), aes_string(x = var, ymin="logcapelinlb", ymax="logcapelinub"), width = 0.3, colour = "black") +
     xlab(paste(var)) +
     ylab("ln (Capelin biomass (ktons))") + 
     #ylim(0,9) +
     theme_bw()
cowplot::plot_grid(p1, p2, p3, labels = c(paste(title)), ncol=2)
}


#' optimGraphs2_all
#'
#' @param ls a list of ice subsets with capelin and pseudocalanus data
#' @param var the first variable of interest (usually tice)
#' @param var2 the second variable of interest
#' @param var2val the value of the second variable - median value
#' @param file_name an evocative name signifying the type of analysis
#' @param saveGraph "yes" or "no" - "no" allows you to view graphs
#'
#' @return a list of output from function optim (parameter estimates, likelihood values, and convergence)
#' @export
#'
#' @examples optimGraphs2_all(MaxTice7, "tice", "Ssurface_tows_lag2", var2val = 280, "opt_div3c", "yes")
optimGraphs2_all <- function(ls, var, var2=NULL, var2val, file_name, saveGraph){
     for(i in 1:length(ls$optim_ls)){
          #browser()
          df1 <- as.data.frame(ls$optim_ls[[i]]$df)
          df2 <- as.data.frame(ls$optim_ls[[i]]$regime1)
          df3 <- as.data.frame(ls$optim_ls[[i]]$regime2)
          mm <- optimGraphs2(df1, df2, df3, yearLim, yearInt, lnbiomassInt,  titlenames[i], var, var2, var2val)
          if(saveGraph == "yes"){
               ggsave(mm, filename = paste0("figs/covariates/", file_name, titlenames[i], ".pdf"), width=10, height=8, units="in")
          }
          if(saveGraph == "no"){
               #x <- optimGraphs2(df1, df2, df3, yearLim, yearInt, lnbiomassInt,  titlenames[i], var, var2, var2val)
               #print(x)
               print(mm)
          }
     }
}


#' optimGraphs2
#'
#' @param df the dataframe provided by the list
#' @param reg1 a set of expected values for pre1990 - deprecated
#' @param reg2 a set of expected values for post 2002
#' @param yearLim the year limits
#' @param yearInt the year intervals
#' @param lnbiomassInt the ln(Biomass) intervats
#' @param title - names of the output graphs
#' @param var - the first variable of interest (usually tice)
#' @param var2 - the second variable of interest
#' @param var2val - the median value of the second variable - determined by x_range in the calcFit function and used in graph
#'
#' @return a set of three graphs for each subset
#' @export
#'
#' @examples nested in optimGraphs2_all(MaxTice7, "tice", "Ssurface_tows_lag2", var2val = 280, "opt_div3c", "yes")
optimGraphs2 <- function(df, reg1, reg2, yearLim, yearInt, lnbiomassInt, title, var, var2, var2val){
     #browser()
     
     p1 <- ggplot(df, aes(x = year, y = logcapelin)) + 
          geom_errorbar(width = 0.3, colour = "black", aes(ymin=logcapelinlb, ymax=logcapelinub)) + 
          geom_point(shape=16, size=3)  +
          geom_line(aes(y=ExpectedLogBiomass), colour="blue", linetype=1, size=1.25) +
          #geom_line(aes(y=ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
          scale_y_continuous(limits = c(-4,10), breaks = lnbiomassInt) +
          scale_x_continuous(limits = yearLim, breaks = yearInt) +
          xlab('Year') +
          ylab('ln (Capelin biomass (ktons))') + 
          theme_bw()
     
     p2 <- ggplot(df, aes(x = year, y = capelin)) +
          geom_errorbar(width = 0.3, colour = "black", aes(ymin=capelinlb, ymax=capelinub)) + 
          geom_point(shape=16, size=3)  +
          geom_line(aes(y=exp(ExpectedLogBiomass)), colour="blue", linetype=1, size=1.25) +
          #geom_line(aes(y=exp(ExpectedLogBiomassOld)), colour="blue", linetype=1, size=1.25) +
          #scale_y_continuous(limits = c(0,8500), breaks = biomassInt) +
          scale_x_continuous(limits = yearLim, breaks = yearInt) +
          xlab('Year') +
          ylab('Capelin biomass (ktons)') + 
          theme_bw() 
#browser()
#p3
     if(!is.null(var2)){
          p3 <- ggplot() + 
               geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomass"), colour="red", linetype=1, size=1.25) + 
               #geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomassOld"), colour="blue", linetype=1, size=1.25) +
#               geom_line(data = reg2, eval(parse(text=var2))==var2val), aes_string(x = var, y = "ExpectedLogBiomass"), colour="green", linetype=1, size=1.25) +
               
               geom_line(data = subset(reg2, eval(parse(text=var2))==var2val), aes_string(x = var, y = "ExpectedLogBiomass"), colour="green", linetype=1, size=1.25) +
               geom_point(data = subset(df, year > 1991), aes_string(x = var, y = "logcapelin"), shape=15, size=3) + 
               geom_errorbar(data = subset(df, year > 1991), aes_string(x = var, ymin="logcapelinlb", ymax="logcapelinub"), width = 0.3, colour = "black") +
               xlab(paste(var)) +
               ylab("ln (Capelin biomass (ktons))") + 
               #ylim(0,9) +
               theme_bw() 
     } else if(is.null(var2)){
     p3 <- ggplot() + 
          geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomass"), colour="red", linetype=1, size=1.25) + 
          geom_point(data = subset(df, year > 1991), aes_string(x = var, y = "logcapelin"), shape=15, size=3) +      geom_errorbar(data = subset(df, year > 1991), aes_string(x = var, ymin="logcapelinlb", ymax="logcapelinub"), width = 0.3, colour = "black") +
          xlab(paste(var)) +
          ylab("ln (Capelin biomass (ktons))") + 
          #ylim(0,9) +
          theme_bw()
     }
#p4
#browser()
     if(!is.null(var2)){
          p4 <- ggplot() + 
               geom_point(data = df, aes_string(x = var2, y = "logcapelin"), shape=15, size=3) +
               geom_errorbar(data = df, aes_string(x = var2, ymin= "logcapelinlb", ymax="logcapelinub"), width = 0.3, colour = "black") +
               xlab(paste(var2)) +
               ylab("ln (Capelin biomass (ktons))") + 
               #ylim(0,9) +
               theme_bw()
     }
     #browser()
     mm <- if(is.null(var2)){
          cowplot::plot_grid(p1, p2, p3, labels = c(paste(title)), ncol=2)
     } else if(!is.null(var2)){
          cowplot::plot_grid(p1, p2, p3, p4, labels = c(paste(title)), ncol=2)          }
return(mm)
}



#' #calcFit (deprecated)
#' performs the optim function for a single dataframe.  CalcFit_all is a wrapper that performs the optim function on a list of dataframes for the different ice subsets
#'
#' @param ls a list of dataframes for the different ice subsets (e.g. m1-m6) with associated capelin data
#' @param titlenames for the different ice subsets (e.g. m1-m6)
#' @param var1 - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
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
calcFit_all1 <- function(ls, titlenames, var1, var2, par, form1 = NULL, form2 = NULL, x1_range, x2_range) {
     optim_ls <- rep(list(list()), length(titlenames))
     names(optim_ls) <- titlenames
     for(i in 1:length(ls)){
          df1 <- as.data.frame(ls[[i]])
          #browser()
          optim <- calcFit1(df1, var1=var1, var2=var2, par=par, form1 = form1, form2 = form2, x1_range = x1_range, x2_range = x2_range)
          optim_ls[[i]] <- optim
     }
     return(list(optim_ls = optim_ls))
}     


#' calcFit1 (deprecated)
#' performs the optim function for a single dataframe or list
#' most param are passed from calcFit_all(); like calcFit but for two variables
#' @param df - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param var1 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
#' @param par - a series of values to help the optim function
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param x_range - range of possible valus on the x-axis
#'
#' @return a list with the new expected values of the caplein, the CapelinDomeFit, regime 1 and regime 2
#' @export
#'
#' @examples 
#' note that this returns a warning "Unknown or uninitialised column: 'par'." which apparently is a tibble problem!!

calcFit1 <- function(df, var1, var2=NULL, par, form1 = NULL, form2 = NULL, x1_range, x2_range) {
     #browser()
     #print(environment())
     CapelinDomeFit <- optim(par = par,
                             dataf = df[which(df$logcapelin!='NA'), c(paste(var1), paste(var2), 'logcapelin')], 
                             form1=form1, form2=form2, var1=var1, var2=var2,
                             fn = SSQCapelinDome1, method=c("L-BFGS-B"), lower = c(0.00001, 0))
     
     # before 2010
#     CapelinDomeFitOld <- optim(par=par,
#                                dataf = df[which(df$year < 2011 & df$logcapelin!='NA'), c('year',  paste(var1), paste(var2), 'logcapelin')], 
#                                form1=form1, form2=form2, var1=var1, var2=var2,
#                                fn = SSQCapelinDome1, method=c("Nelder-Mead"))
     
     ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
     #browser()
     df$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit$par), dataf = df[,c(paste(var1), paste(var2))], form1, form2, var1, var2)
     
 #    df$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld$par), dataf = df[,c('year', paste(var1), paste(var2))], form1, form2, var1, var2)
     
     # attach the optimization curves of capelin abundance to ice data
     
     #browser()
     xtice <- expand.grid(col1 = as.numeric(paste(x1_range)), col2=as.numeric(x2_range))
     
     #    xtice <- expand.grid(year = c(1990,2000), col2 = as.numeric(paste(x1_range)), col3=as.numeric(x2_range))
     colnames(xtice)[1] <- c(paste(var1))
     colnames(xtice)[2] <- c(paste(var2))
     xtice <- xtice[order(xtice[1]),]
#     xtice$ExpectedLogBiomassOld <- CapelinDome1(params = c(CapelinDomeFitOld$par), dataf = xtice, form1, form2, var1, var2)
     #browser()
     xtice$ExpectedLogBiomass <- CapelinDome1(params = c(CapelinDomeFit$par),dataf = xtice, form1, form2, var1, var2)
     #not sure what these are for but used in plots below but creates a data set where all values of year are the same????
     #regime1 <- xtice[which(xtice$year == 1990),]
     #regime2 <- xtice[which(xtice$year == 2000),]
     regime2 <- xtice
     #return(list(df = df, cdf = CapelinDomeFit, regime1 = regime1, regime2 = regime2))  
     return(list(df = df, cdf = CapelinDomeFit, regime2 = regime2))  
}


#' SSQCapelinDome1(deprecated)
## Objective function - part of optimization function
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
SSQCapelinDome1 <- function(params, dataf, form1, form2, var1, var2){
     #browser()
     dataf <- as.data.frame(dataf) # needed because optim doesn't work with tibbles.....grrrrrrrr!!!!!
     Alpha <- params[1]
     Beta <- params[2]
     Gamma <- params[3]
     #Delta <- params[4]
     #year <- dataf[,1]
     #tice <- dataf[,2]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     tmp2 <- dataf[,2] #variable of interest in form1 and form2
     assign(var2, tmp2)
     logcap <- dataf[,3]
     # this is based on MSY
     # ELogCapBiom must be positive
     # It is - always
     ELogCapBiom <- eval(parse(text = form1))
#     ELogCapBiom <- ifelse(year < 2017, eval(parse(text = form1)), eval(parse(text = form2)))
     sum((logcap-ELogCapBiom)^2)
}


#' CapelinDome1(depredated)
#' Function to obtain Expected Log Capelin Biomass    
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
CapelinDome1 <- function(params, dataf, form1, form2, var1, var2){
     #browser()
     dataf <- as.data.frame(dataf)
     Alpha <- params[1]
     Beta <- params[2]
     Gamma <- params[3]
     #Delta <- params[4]
     #year <- dataf[,1]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     tmp2 <- dataf[,2] #variable of interest in form1 and form2
     assign(var2, tmp2)
     ELogCapBiom <- eval(parse(text = form1))
#     ELogCapBiom <- ifelse(year<2017, eval(parse(text = form1)), eval(parse(text = form2)))
     ELogCapBiom
}


#' #calcFit performs the optim function for a single dataframe.  CalcFit_all is a wrapper that performs the optim function on a list of dataframes for the different ice subsets
#'
#' @param ls a list of dataframes for the different ice subsets (e.g. m1-m6) with associated capelin data
#' @param titlenames for the different ice subsets (e.g. m1-m6)
#' @param var1 - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
#' #' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "pseudo calanus density"
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
calcFit_all2 <- function(ls, titlenames, var1, var2, var3=NULL, par, 
                         form1 = NULL, form2 = NULL, 
                         x1_range, x2_range=NULL, x3_range=NULL, 
                         method = c("L-BFGS-B"), lowerLim = "no") {
     
     optim_ls <- rep(list(list()), length(titlenames))
     names(optim_ls) <- titlenames
     for(i in 1:length(ls)){
          df1 <- as.data.frame(ls[[i]])
          #browser()
          optim <- calcFit2(df1, var1=var1, var2=var2, var3=var3, par=par, form1 = form1, form2 = form2, x1_range = x1_range, x2_range = x2_range, x3_range = x3_range)
          optim_ls[[i]] <- optim
     }
     return(list(optim_ls = optim_ls))
}     


#' calcFit1 performs the optim function for a single dataframe or list
#' most param are passed from calcFit_all(); like calcFit but for two variables
#' @param df - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param var1 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
#' @param par - a series of values to help the optim function
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param x_range - range of possible valus on the x-axis
#'
#' @return a list with the new expected values of the caplein, the CapelinDomeFit, regime 1 and regime 2
#' @export
#'
#' @examples 
#' note that this returns a warning "Unknown or uninitialised column: 'par'." which apparently is a tibble problem!!

calcFit2 <- function(df, var1, var2=NULL, var3=NULL, par, 
                     form1 = NULL, form2 = NULL, 
                     x1_range, x2_range=NULL, x3_range=NULL, 
                     method = c("BFGS"),
                     lowerLim = "no") {
     
     #browser()
     #print(environment())
     if(lowerLim == "yes"){
          CapelinDomeFit <- optim(par = par,
                                  dataf = df[which(df$logcapelin!='NA'), c(paste(var1), paste(var2), paste(var3), 'logcapelin')], #age2_log
                                  form1=form1, form2=form2, var1=var1, var2=var2, var3=var3,
                                  fn = SSQCapelinDome2, method=method, lower = c(0.00001, 0))
          
     }
     if(lowerLim == "no"){
          CapelinDomeFit <- optim(par = par,
                                  dataf = df[which(df$logcapelin!='NA'), c(paste(var1), paste(var2), paste(var3), 'logcapelin')], #logcapelin
                                  form1=form1, form2=form2, var1=var1, var2=var2, var3=var3, 
                                  fn = SSQCapelinDome2, method=method)
          
     }
     ## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
     #browser()
     if(is.null(var3)){
          df$ExpectedLogBiomass <- CapelinDome2(params = c(CapelinDomeFit$par), dataf = df[,c(paste(var1), paste(var2))], form1, form2, var1, var2, var3)
          
     } else {
          df$ExpectedLogBiomass <- CapelinDome2(params = c(CapelinDomeFit$par), dataf = df[,c(paste(var1), paste(var2), paste(var3))], form1, form2, var1, var2, var3)
     }
     
     #browser()
     xtice <- expand.grid(col1 = as.numeric(paste(x1_range)), col2=as.numeric(x2_range), col3=as.numeric(x3_range))
     colnames(xtice)[1] <- c(paste(var1))
     colnames(xtice)[2] <- c(paste(var2))
     colnames(xtice)[3] <- c(paste(var3))
     xtice <- xtice[order(xtice[1]),]
     xtice$ExpectedLogBiomass <- CapelinDome2(params = c(CapelinDomeFit$par),dataf = xtice, form1, form2, var1, var2, var3)
     regime2 <- xtice
     return(list(df = df, cdf = CapelinDomeFit, regime2 = regime2))  
}


#' SSQCapelinDome2
#'
#' @param params 
#' @param dataf 
#' @param form1 
#' @param form2 
#' @param var1 
#' @param var2 
#' @param var3 
#'
#' @return
#' @export
#'
#' @examples
SSQCapelinDome2 <- function(params, dataf, form1, form2, var1, var2, var3){
     #browser()
     dataf <- as.data.frame(dataf) # needed because optim doesn't work with tibbles.....grrrrrrrr!!!!!
     Alpha <- params[1]
     Beta <- params[2]
     Gamma <- params[3]
     Delta <- params[4]
     #year <- dataf[,1]
     #tice <- dataf[,2]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     tmp2 <- dataf[,2] #variable of interest in form1 and form2
     assign(var2, tmp2)
     tmp3 <- dataf[,3] #variable of interest in form1 and form2      
     assign(var3, tmp3)
     logcap <- dataf[,4]
     # this is based on MSY
     # ELogCapBiom must be positive
     # It is - always
     ELogCapBiom <- eval(parse(text = form1))
     #     ELogCapBiom <- ifelse(year < 2017, eval(parse(text = form1)), eval(parse(text = form2)))
     sum((logcap-ELogCapBiom)^2)
}


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
CapelinDome2 <- function(params, dataf, form1, form2, var1, var2, var3){
     #browser()
     dataf <- as.data.frame(dataf)
     Alpha <- params[1]
     Beta <- params[2]
     Gamma <- params[3]
     Delta <- params[4]
     #year <- dataf[,1]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     tmp2 <- dataf[,2] #variable of interest in form1 and form2
     assign(var2, tmp2)
     tmp3 <- dataf[,3] #variable of interest in form1 and form2
     assign(var3, tmp3)
     ELogCapBiom <- eval(parse(text = form1))
     #     ELogCapBiom <- ifelse(year<2017, eval(parse(text = form1)), eval(parse(text = form2)))
     ELogCapBiom
}

optimSummary <- function(ls, titlenames){
     output <- NULL
     for(i in seq_along(ls$optim_ls)){
          par <- ls$optim_ls[[i]]$cdf$par
          par <- as.data.frame(t(par))
          if(length(par)==2) {par <- rename(par, par1 = V1, par2 = V2)}
          par
          if(length(par)==3) {par <- rename(par, par1 = V1, par2 = V2, par3 = V3)}
          if(length(par)==4) {par <- rename(par, par1 = V1, par2 = V2, par3 = V3, par4 = V4)}
          val <- ls$optim_ls[[i]]$cdf$value
          con <- ls$optim_ls[[i]]$cdf$convergence
          mname <- titlenames[i]
          dat <- as.data.frame(cbind(mname, par, val, con))
          output <- rbind(output, dat)
     }     
     return(output)
}


#' #calcFit performs the optim function for a single dataframe.  CalcFit_all is a wrapper that performs the optim function on a list of dataframes for the different ice subsets
#'
#' @param ls a list of dataframes for the different ice subsets (e.g. m1-m6) with associated capelin data
#' @param titlenames for the different ice subsets (e.g. m1-m6)
#' @param var1 - the variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
#' #' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "pseudo calanus density"
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
calcFit_all3 <- function(ls, titlenames, var1, var2=NULL, var3=NULL, rv, par, 
                         form1 = NULL, form2 = NULL, 
                         x1_range, x2_range=NULL, x3_range=NULL, 
                         method = c("L-BFGS-B"), lowerLim = "no") {
     
     optim_ls <- rep(list(list()), length(titlenames))
     names(optim_ls) <- titlenames
     for(i in 1:length(ls)){
          df1 <- as.data.frame(ls[[i]])
          #browser()
          optim <- calcFit3(df1, var1=var1, var2=var2, var3=var3, rv=rv, par=par, form1 = form1, form2 = form2, x1_range = x1_range, x2_range = x2_range, x3_range = x3_range)
          optim_ls[[i]] <- optim
     }
     return(list(optim_ls = optim_ls))
}     


#' calcFit1 performs the optim function for a single dataframe or list
#' most param are passed from calcFit_all(); like calcFit but for two variables
#' @param df - a dataframe for a single ice subset (e.g. Ale's data) with associated capelin data or part of a list used in calcFit_all
#' @param var1 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "tice"
#' @param var2 - a variable of interest to enter into the model, i.e. form1 and form2, e.g. "surface_tows_lag2"
#' @param par - a series of values to help the optim function
#' @param form1 - a formula to describe the relationship of some ice parameter to capelin abundance post 1991
#' @param form2 - a formula to describe the relationship of some ice parameter to capelin abundance  pre-1991; has the gamma param
#' @param x_range - range of possible valus on the x-axis
#'
#' @return a list with the new expected values of the caplein, the CapelinDomeFit, regime 1 and regime 2
#' @export
#'
#' @examples 
#' note that this returns a warning "Unknown or uninitialised column: 'par'." which apparently is a tibble problem!!

calcFit3 <- function(df, var1, var2=NULL, var3=NULL, rv, par, 
                     form1 = NULL, form2 = NULL, 
                     x1_range, x2_range=NULL, x3_range=NULL, 
                     method = c("BFGS"),
                     lowerLim = "no") {
     #browser()
     #print(environment())
     if(lowerLim == "yes"){
          CapelinDomeFit <- optim(par = par,
               dataf = df[which(df$logcapelin!='NA'), 
               c(paste(var1), paste(var2), paste(var3), 
                 paste(rv))], 
               form1=form1, form2=form2, var1=var1, 
               var2=var2, var3=var3, rv = rv,
               fn = SSQCapelinDome3, method=method, 
               lower = c(0.00001, 0))
          
     }
     if(lowerLim == "no"){
          CapelinDomeFit <- optim(par = par,
                                  dataf = df[which(df$logcapelin!='NA'), c(paste(var1), paste(var2), paste(var3), paste(rv))], 
                                  form1=form1, form2=form2, var1=var1, var2=var2, var3=var3, rv = rv,
                                  fn = SSQCapelinDome3, method=method)
     }

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
     #browser()
     if(is.null(var3)){
          df$ExpectedLogBiomass <- CapelinDome3(params = c(CapelinDomeFit$par), dataf = df[,c(paste(var1), paste(var2))], form1, form2, var1, var2, var3)
     } else {
          df$ExpectedLogBiomass <- CapelinDome3(params = c(CapelinDomeFit$par), dataf = df[,c(paste(var1), paste(var2), paste(var3))], form1, form2, var1, var2, var3)
     }
     
     #browser()
     if(is.null(var3) && is.null(var2)){
          xtice <- expand.grid(col1 = as.numeric(x1_range))
     } else if(is.null(var3)) {
     xtice <- expand.grid(col1 = as.numeric(x1_range), col2=as.numeric(x2_range))
     } else {
     xtice <- expand.grid(col1 = as.numeric(x1_range), col2=as.numeric(x2_range), col3=as.numeric(x3_range))
     }
     colnames(xtice)[1] <- c(paste(var1))
     if(!is.null(var2)){colnames(xtice)[2] <- c(paste(var2))}
     if(!is.null(var3)){colnames(xtice)[3] <- c(paste(var3))}
     #xtice <- xtice[order(xtice[1]),]
     #browser()
     xtice$ExpectedLogBiomass <- CapelinDome3(params = c(CapelinDomeFit$par),dataf = xtice, form1, form2, var1, var2, var3)
     regime2 <- xtice
     return(list(df = df, cdf = CapelinDomeFit, regime2 = regime2))  
}


#' SSQCapelinDome2
#'
#' @param params 
#' @param dataf 
#' @param form1 
#' @param form2 
#' @param var1 
#' @param var2 
#' @param var3 
#'
#' @return
#' @export
#'
#' @examples
SSQCapelinDome3 <- function(params, dataf, form1, 
                            form2, var1, var2=NULL, var3=NULL, rv){
     #browser()
     dataf <- as.data.frame(dataf) # needed because optim doesn't work with tibbles.....grrrrrrrr!!!!!
     Alpha <- params[1]
     Beta <- params[2]
     if(!is.null(var2)){Gamma <- params[3]}
     if(!is.null(var3)){Delta <- params[4]}
     #year <- dataf[,1]
     #tice <- dataf[,2]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     if(!is.null(var2)){
          tmp2 <- dataf[,2] #variable of interest in form1 and form2
          assign(var2, tmp2)
     }
     if(!is.null(var3)){
          tmp3 <- dataf[,3] #variable of interest in form1 and form2      
          assign(var3, tmp3)
     }
     logcap <- dataf[,paste(rv)]
     # this is based on MSY
     # ELogCapBiom must be positive
     # It is - always
     ELogCapBiom <- eval(parse(text = form1))
     #     ELogCapBiom <- ifelse(year < 2017, eval(parse(text = form1)), eval(parse(text = form2)))
     sum((logcap-ELogCapBiom)^2)
}


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
CapelinDome3 <- function(params, dataf, form1, form2, var1, var2=NULL, var3=NULL, rv){
     #browser()
     dataf <- as.data.frame(dataf)
     Alpha <- params[1]
     Beta <- params[2]
     if(!is.null(var2)){Gamma <- params[3]}
     if(!is.null(var3)){Delta <- params[4]}
     #year <- dataf[,1]
     tmp1 <- dataf[,1] #variable of interest in form1 and form2
     assign(var1, tmp1)
     if(!is.null(var2)){
          tmp2 <- dataf[,2] #variable of interest in form1 and form2
          assign(var2, tmp2)
     }
     if(!is.null(var3)){
          tmp3 <- dataf[,3] #variable of interest in form1 and form2
          assign(var3, tmp3)
     }
          
     ELogCapBiom <- eval(parse(text = form1))
     #     ELogCapBiom <- ifelse(year<2017, eval(parse(text = form1)), eval(parse(text = form2)))
     ELogCapBiom
}


#' var2val - helps to give value for the var2val arguement in the optimGraphs2_all
#'
#' @param df 
#'
#' @return - median values for 2nd variable + range of values
#' @export
#'
#' @examples
var2val <- function(df){
     med <- median(df)
     lev <- levels(as.factor(df))
     return(list(med=med, lev=lev))    
}