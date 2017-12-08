#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-08, R version 3.3.3 (2017-03-06)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-11.08 #


# The purpose of this file is to:
# 1) Use function optim to test additional covariates to improve the Tice model fit
# 2) focus on capelin recruitment and pseudocalanus abundance.

rm(list=ls())

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(plotly)

## read in source code-----
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
#source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")


## load Ale original data----
capelin <- read.csv('capelin-ice-2014.csv',header=T)
# get rid of extra columns in csv file
capelin <- capelin[c("year", "maxarea", "minlat", "tice", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]
capelin_ale <- capelin
# from "acoustic_estimates_2017.xlsx" in capelin folder
# add capelin, years, and lncapelin
capelin[c(47, 48, 49), 5] <- matrix(c(662, NA, 158), ncol = 1) 
capelin[c(47, 48, 49), 1] <- matrix(c(2015, 2016, 2017), ncol = 1) 
capelin[c(47, 48, 49), 8] <- matrix(c(log(662), NA, log(158)), ncol = 1)

# data set for joining to m1-m6 datasets
capelin_join <- capelin[c("year", "capelin", "capelinlb", "capelinub", "logcapelin", "logcapelinlb", "logcapelinub")]

# read in larval data
larvae <- read.csv('capelin_larval_indices.csv',header=T)
larvae$surface_tows_lag2 <- lag(larvae$surface_tows, 2)
#normalize the surface_tows 
larvae$Nsurface_tows_lag2 <- (larvae$surface_tows_lag2 - mean(larvae$surface_tows_lag2, na.rm=T)) / sd(larvae$surface_tows_lag2, na.rm = T)


# read in Pseudocalanus data
## load original data----
pscal <- read_csv("data/pseudocal_1999_2016.csv")
str(pscal)

# filter out year and use just the "total" stage.  Calcuate mean and sd by year
ps_tot <- pscal %>%
     gather(key = stage, value = density, 4:10) %>%
     filter(doy <= 274 & doy >= 152) %>%
     filter(stage == "total") %>%
     group_by(year) %>%
     summarise(ps_meanTot = mean(density), ps_sdTot = sd(density))
# lage mean and sd and /1000 to scale to tice
ps_tot$ps_meanTot_lag1 <- lag(ps_tot$ps_meanTot, 1)/1000
ps_tot$ps_sdTot_lag1 <- lag(ps_tot$ps_sdTot, 1)/1000
range(ps_tot$ps_meanTot_lag1, na.rm = TRUE)

ps_tot$ps_meanTot_lag2 <- lag(ps_tot$ps_meanTot, 2)/1000
ps_tot$ps_sdTot_lag2 <- lag(ps_tot$ps_sdTot, 2)/1000

# get rid of ex
#join the acoustic and larvae data
capelin_join <- left_join(capelin_join, larvae, by = "year")
capelin_join <- left_join(capelin_join, ps_tot, by = "year")
str(capelin_join)

### load files generated in area-ice----
## max area
# look at files in "output-processing" and create a file pattern
data_path <- "output-processing"
files <- list.files(data_path)
pattern <- grep("capelin", files, value = T)
# read the datasets and join them to capelin_join
cape <- loadSubsetDatasets1(df = capelin_join, name = "capelin_m", pat = pattern, 6, var1 = "area", var2 = "minlats", nvar1 = "max_area", nvar2 = "minlats")

# Divide area by 1000 and normalize Ntice
for(i in 1:length(cape)){
     cape[[i]]$max_area1000 <- cape[[i]]$max_area/1000
     cape[[i]]$Ntice <- ((cape[[i]]$tice - mean(cape[[i]]$tice))/sd(cape[[i]]$tice)) + 5
     cape[[i]]$logsurface_tows_lag2 <- log(cape[[i]]$surface_tows_lag2)
     cape[[i]]$logtice <- log(cape[[i]]$tice)
     cape[[i]]$Htice <- cape[[i]]$tice*100
     cape[[i]]$Ssurface_tows_lag2 <- cape[[i]]$surface_tows_lag2/10
}
range(cape$capelin_m1$Ssurface_tows_lag2, na.rm = T)

# seperate cape into time components
cape_1991 <- map(cape, ~filter(.x, year <= 1991))
cape_2017 <- map(cape, ~filter(.x, year > 1991))
cape_2001 <- map(cape, ~filter(.x, year > 2002 & year < 2017))


## EXPLORATORY ANALYSIS----
# make a simple data set and plot for m1
sdf <- cape_2001$capelin_m1[c("year", "tice", "logcapelin", "surface_tows", "surface_tows_lag2", "ps_meanTot", "ps_meanTot_lag2")]

str(cape_2001$capelin_m1)
plot(sdf$tice, sdf$Ssurface_tows_lag2)
plot(sdf$tice, sdf$ps_meanTot_lag1)
plot(sdf$Ssurface_tows_lag2, sdf$ps_meanTot_lag1)

library(plot3D)
View(sdf)
sdf$logTows <- log10(sdf$surface_tows_lag2)
max(sdf$logTows)/min(sdf$logTows)
scatter3D(sdf$Ntice, sdf$logcapelin, sdf$Nsurface_tows_lag2)

windows()
ggplot(data=cape_2001$capelin_m6) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=tice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)

---------------------------------------------------------------
## values for optim functions----
yearInt <- seq(2000, 2017, by=4)
yearLim <- c(2000, 2017)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures


# replicate graphs with calcFit_all - gives extra lines
MaxTice <- calcFit_all(cape_2001, titlenames, par = c(1, 200, 0.6), var = "tice", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))",
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:190,173.515,187.768))
str(MaxTice, max.level = 3)

optimGraphs1_all(MaxTice, "tice", "opt_")


## new values for optim functions----
yearInt <- seq(2002, 2017, by=3)
yearLim <- c(2002, 2017)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")

##MaxTice1 - raw data, 2 parm, 1 var, dome----
# replicate graphs with calcFit_all1 and only 2000 data
MaxTice1a <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "tice", var2 = "surface_tows_lag2",
                       form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                       x1_range = c(0:150),
                       x2_range = seq(500, 4000, 25))
optimSummary(MaxTice1a, titlenames = titlenames)
rawTice <- optimSummary(MaxTice1a, titlenames = titlenames)
optimGraphs2_all(MaxTice1a, "tice", var2 = NULL, var2val = NULL, "rawTice", "no")

# Surface tows
MaxTice1b <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "Ssurface_tows_lag2", var2 = "tice",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         x1_range = c(0:500),
                         x2_range = c(0:150))
scalST <- optimSummary(MaxTice1b, titlenames = titlenames)
optimGraphs2_all(MaxTice1b, "Ssurface_tows_lag2", var2 = NULL, var2val = NULL, "scalST", "no")

# Pseudo calanus
MaxTice1c <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "ps_meanTot_lag1", var2 = "tice",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                          x1_range = seq(10, 200, 20),
                          x2_range = c(0:500))
scalPS <- optimSummary(MaxTice1c, titlenames = titlenames)
optimGraphs2_all(MaxTice1c, "ps_meanTot_lag1", var2 = NULL, var2val = NULL, "scalPS", "no")

# Surface tows_linear model
MaxTice1d <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "Ssurface_tows_lag2", var2 = "tice",
                          form1 = "Alpha*tmp1",
                          x1_range = c(0:500),
                          x2_range = c(0:150))
scalST_linear <- optimSummary(MaxTice1d, titlenames = titlenames)
optimGraphs2_all(MaxTice1d, "Ssurface_tows_lag2", var2 = NULL, var2val = NULL, "scalST", "no")

# Pseudo calanus_linear model
MaxTice1e <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "ps_meanTot_lag1", var2 = "tice",
                          form1 = "Alpha*tmp1",
                          x1_range = seq(10, 200, 20),
                          x2_range = c(0:500))
scalPS_linear <- optimSummary(MaxTice1e, titlenames = titlenames)
optimGraphs2_all(MaxTice1e, "ps_meanTot_lag1", var2 = NULL, var2val = NULL, "scalPS", "no")

##MaxTice2a - normalized data, 2 parm, 1 var, dome----
MaxTice2a <- calcFit_all1(cape_2001, titlenames, par = c(1, 7), var1 = "Ntice", var2 = "Nsurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         x1_range = seq(0, 10, 0.067), 
                         x2_range = c(-3:3))

normTice <- optimSummary(MaxTice2a, titlenames = titlenames)
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
optimGraphs2_all(MaxTice2a, "Ntice", var2 = NULL, var2val = 0, "normTice", "no")

optimGraphs2_all(MaxTice2a, "Ntice", var2 = "Nsurface_tows_lag2", var2val = 0, "normTice", "no")
MaxTice2a$optim_ls$`MaxTice-m1`$df


### MaxTice4 log data, 2 parm, 1 var, dome----
MaxTice4 <- calcFit_all1(cape_2001, titlenames, par = c(1, 1), var1 = "logtice", var2 = "logsurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         x1_range = seq(3, 5, 0.1),
                         x2_range = seq(3, 5, 0.1))

logTice <- optimSummary(MaxTice4, titlenames = titlenames)
optimGraphs2_all(MaxTice4, "logtice", var2 = NULL, var2val = 4.0, file_name = "logTice", "no")


### MaxTice6 - scaled data, 3 parm, 2 var, dome + line----
MaxTice6 <- calcFit_all2(cape_2001, titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2", var3 ="ps_meanTot_lag1", 
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50), 
                         x3_range = c(1), 
                         lowerLim = "no")

scalTiceSTow <- optimSummary(MaxTice6, titlenames = titlenames)
optimGraphs2_all(MaxTice6, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTow", saveGraph = "no")


### MaxTice6a----
# psuedocalanus
#scaled data, 3 parm, 2 var, dome + line
MaxTice6a <- calcFit_all2(cape_2001, titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", var2 = "ps_meanTot_lag1", var3 = "Ssurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(10, 200, 20),
                         x3_range = c(1))

scalTicePS <- optimSummary(MaxTice6a, titlenames = titlenames)
optimGraphs2_all(MaxTice6a, "tice", var2 = "ps_meanTot_lag1", var2val = 90, file_name = "scalTicePS", saveGraph = "no")

### MaxTice6b scaled data, 3 parm, 2 var, dome + line----
MaxTice6b <- calcFit_all2(cape_2001, titlenames, 
                          par = c(1, 200, 1), 
                          var1 = "Ssurface_tows_lag2", var2 ="ps_meanTot_lag1", var3 = "tice",
                          form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                          #form2 = "Alpha*tmp2*Beta*Gamma",
                          x1_range = seq(30, 500, 50),
                          x2_range = seq(10, 200, 20),
                          x3_range = c(1))

scalSTowPS <- optimSummary(MaxTice6b, titlenames = titlenames)
optimGraphs2_all(MaxTice6b, "Ssurface_tows_lag2", var2 = "ps_meanTot_lag1", var2val = 90, file_name = "scalSTowPS", saveGraph = "no")

### MaxTice6c - Holling IV: scaled data, 3 parm, 2 var ----
#DON'T CHANGE 
MaxTice6c <- calcFit_all1(cape_2001, titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2",
                         form1 = "Alpha*tmp1^2/(Beta + Gamma*tmp1+tmp1^2)",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50))

scalTiceST_HollingIV <- optimSummary(MaxTice6c, titlenames = titlenames)
optimGraphs2_all(MaxTice6c, "tice", var2 = NULL, var2val = 230, file_name = "scalTiceST_HollingIV", saveGraph = "no")


### MaxTice6d - Holling IV: scaled data, 3 parm, 2 var ----
MaxTice6d <- calcFit_all2(cape_2001, titlenames, 
                          par = c(1, 200, 1, 1), 
                          var1 = "tice", var2 = "Ssurface_tows_lag2", var3 ="ps_meanTot_lag1",
                          form1 = "Alpha*tmp1^2/(Beta + Gamma*tmp1+tmp1^2)+ Delta*tmp2",
                          #form2 = "Alpha*tmp2*Beta*Gamma",
                          x1_range = seq(0, 150, 10),
                          x2_range = seq(30, 500, 50),
                          x3_range = c(1))

    
scalTiceST_HollingIV_line <- optimSummary(MaxTice6d, titlenames = titlenames)
optimGraphs2_all(MaxTice6d, "tice", var2 = NULL, var2val = 230, file_name = "scalTiceST_HollingIV", saveGraph = "no")


### MaxTice8 scaled data, 4 parm, 3 var, dome + 2line----
MaxTice8 <- calcFit_all2(cape_2001, titlenames, 
                         par = c(1, 200, 1, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2", var3 ="ps_meanTot_lag1",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50), 
                         x3_range = seq(10, 200, 20))


scalTiceSTowPS <- optimSummary(MaxTice8, titlenames = titlenames)
optimGraphs2_all(MaxTice8, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowPS", saveGraph = "no")


MaxTice8a <- calcFit_all2(cape_2001, titlenames, 
                         par = c(1, 200, 1, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2", var3 ="ps_meanTot_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2 + Delta*tmp3",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50), 
                         x3_range = seq(10, 200, 20))
#c(250))
scalTiceSTowPS_lag2 <- optimSummary(MaxTice8a, titlenames = titlenames)
optimGraphs2_all(MaxTice8a, "tice", var2 = "Ssurface_tows_lag2", var2val = 230, file_name = "scalTiceSTowPS_l2", saveGraph = "no")

# Table of optimSummary results
optimSummary_ls <- list(rawTice=rawTice, 
                        scalST = scalST, 
                        scalPS = scalPS, 
                        normTice = normTice, 
                        logTice = logTice, 
                        scalTiceST_HollingIV = scalTiceST_HollingIV,
                        scalTiceSTow = scalTiceSTow, 
                        scalTicePS = scalTicePS, 
                        scalSTowPS = scalSTowPS,
                        scalTiceSTowPS = scalTiceSTowPS)
optimSummary_df <- do.call("bind_rows", lapply(optimSummary_ls, as.data.frame))
write_csv(optimSummary_df, "optimSummary_df.csv")

test <- data.frame(ID = rep(names(optimSummary_ls), sapply(optimSummary_ls, length)), Obs=unlist(optimSummary_ls))

spread(test, "ID")
data.table::rbindlist(optimSummary_ls, fill = TRUE, use.names = TRUE)

library(data.table)
DT1 = data.table(A=1:3,B=letters[1:3])
DT2 = data.table(B=letters[4:5],A=4:5)
l = list(DT1,DT2)
rbindlist(l, use.names=TRUE)

DT1 = data.table(A=1:3,B=letters[1:3])
DT2 = data.table(B=letters[4:5],C=factor(1:2))
l = list(DT1,DT2)
rbindlist(l, use.names=TRUE, fill=TRUE)

# generate index column, auto generates indices
rbindlist(l, use.names=TRUE, fill=TRUE, idcol=TRUE)
# let's name the list
setattr(l, 'names', c("a", "b"))
rbindlist(l, use.names=TRUE, fill=TRUE, idcol="ID")

################################END#################################################
## med area -----

# read the datasets and join them to capelin_join
pattern <- grep("sub2017-m+[^rsq].csv", files, value = T)
med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "med_m", pat = pattern, N = 6, var1 = "darea", var2 = "dminlats", nvar1 = "med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(med_cape_2017)){
     med_cape_2017[[i]]$med_area1000 <- med_cape_2017[[i]]$med_area/1000
}

# seperate cape into time components
med_cape_2001 <- map(med_cape_2017, ~filter(.x, year > 2001))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- med_cape_2001$med_m1[c("year", "dtice", "logcapelin", "surface_tows_lag2")]

windows()
ggplot(data=med_cape_2001$med_m1) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=dtice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)

---------------------------------------------------------------
## D5 med area----
pattern <- grep("iceMedD5p2017-m.a", files, value = T)

d5med_cape_2017 <- loadSubsetDatasets1(df = capelin_join, name = "d5med_m", pat = pattern, N = 6, var1 = "d5area", var2 = "d5minlats", nvar1 = "d5med_area", nvar2 = "minlats")

# Divide area by 1000
for(i in 1:length(d5med_cape_2017)){
     d5med_cape_2017[[i]]$d5med_area1000 <- d5med_cape_2017[[i]]$d5med_area/1000
}

# seperate cape into time components
d5med_cape_2001 <- map(d5med_cape_2017, ~filter(.x, year > 2001))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- d5med_cape_2001$d5med_m1[c("year", "d5tice", "logcapelin", "surface_tows_lag2")]

windows()
ggplot(data=d5med_cape_2001$d5med_m1) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
     geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=d5tice), size=3) + 
     geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)

titlenames <- c("D5MedArea-m1", "D5MedArea-m2", "D5MedArea-m3", "D5MedArea-m4", "D5MedArea-m5", "D5MedArea-m6")

