#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-11-08, R version 3.3.3 (2017-03-06)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-11.08 #


# The purpose of this file is to:
# 1) Use function optim to test additional covariates to improve the Tice model fit

rm(list=ls())

## libraries------
library(plotrix)
library(ggplot2)
library(readr)
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

#join the acoustic and larvae data
capelin_join <- left_join(capelin_join, larvae, by = "year")

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

# seperate cape into time components
cape_1991 <- map(cape, ~filter(.x, year <= 1991))
cape_2017 <- map(cape, ~filter(.x, year > 1991))
cape_2001 <- map(cape, ~filter(.x, year > 2002 & year < 2017))

# EXPLORATORY ANALYSIS
# make a simple data set and plot for m1
sdf <- cape_2001$capelin_m1[c("year", "tice", "logcapelin", "surface_tows_lag2", "Ntice", "Nsurface_tows_lag2")]
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
## optimization----
# values for optim functions
yearInt <- seq(2000, 2017, by=4)
yearLim <- c(2000, 2017)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

labtice <- expression(paste(italic(t[ice]), '(day of year)')) # label for figures
titlenames <- c("MaxTice-m1", "MaxTice-m2", "MaxTice-m3", "MaxTice-m4", "MaxTice-m5", "MaxTice-m6")

# replicate graphs with calcFit_all - gives extra lines
MaxTice <- calcFit_all(cape_2001, titlenames, par = c(1, 200, 0.6), var = "tice", 
                       form1 = "Alpha*tmp*(1-(tmp/Beta))",
                       form2 = "Alpha*tmp*(1-(tmp/Beta))*Gamma",
                       x_range = c(0:190,173.515,187.768))
str(MaxTice, max.level = 3)

optimGraphs1_all(MaxTice, "tice", "opt_")

##MaxTice1----
# replicate graphs with calcFit_all1 and only 2000 data
MaxTice1 <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "tice", var2 = "surface_tows_lag2",
                       form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
#                       form2 = "Alpha*tmp1*(1-(tmp1/Beta))*Gamma",
                       x1_range = c(0:150),
                       x2_range = seq(500, 4000, 25))
str(MaxTice1, max.level = 4)
head(MaxTice1$optim_ls$`MaxTice-m1`$df)
MaxTice1$optim_ls$`MaxTice-m1`$df
MaxTice1$optim_ls$`MaxTice-m1`$cdf

optimGraphs1_all(MaxTice1, "tice", "opt1")


## works above here----
MaxTice2 <- calcFit_all1(cape_2001, titlenames, par = c(1, 7), var1 = "Ntice", var2 = "Nsurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         #form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + tmp2*Gamma",
#                         form2 = "Alpha*tmp1*(1-(tmp1/Beta))*Gamma",
                         x1_range = seq(0, 10, 0.067), #c(0:10),
                         x2_range = c(-3:3))
str(MaxTice, max.level = 3)
head(MaxTice2$optim_ls$`MaxTice-m1`$regime2)

MaxTice2$optim_ls$`MaxTice-m1`$cdf
MaxTice1$optim_ls$`MaxTice-m1`$cdf

# this suggests to me that the normalizatoin is squishing the variable bc ratios are different
x <- MaxTice2$optim_ls$`MaxTice-m1`$df
rm(x)
min(x$tice) - mean(x$tice)
range(x$tice)
max(x$tice)/min(x$tice)
max(x$tice)
min(x$Ntice) - mean(x$Ntice)
range(x$Ntice)
max(x$Ntice)/min(x$Ntice)
max(x$Ntice)

#create and save graphs
#source("D:/Keith/capelin/2017-project/ice-capelin-functions.R")
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")

optimGraphs1_all(MaxTice4, "Ntice", "opt_norm")

### MaxTice4 log----
yearInt <- seq(2002, 2017, by=3)
yearLim <- c(2002, 2017)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")

MaxTice4 <- calcFit_all1(cape_2001, titlenames, par = c(1, 1), var1 = "logtice", var2 = "logsurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         #                       form2 = "Alpha*tmp1*(1-(tmp1/Beta))*Gamma",
                         x1_range = seq(3, 5, 0.1),
                         x2_range = seq(3, 5, 0.1))

MaxTice4$optim_ls$`MaxTice-m1`$cdf
optimGraphs1_all(MaxTice4, "logtice", file_name = "optLog2")

str(MaxTice4$optim_ls$`MaxTice-m5`$regime2)
head(MaxTice4$optim_ls$`MaxTice-m4`$df)
x <- MaxTice4$optim_ls$`MaxTice-m1`$regime2
plot(x$logtice, x$ExpectedLogBiomass)

head(x)
tail(x)

#attempt at back transforming
for(i in 1:length(MaxTice4$optim_ls)){
     for(col in names(MaxTice4$optim_ls[[i]]$regime2)){
          new_names <- c(1:3, NA)
          new_names <- paste0("E", names(MaxTice4$optim_ls[[i]]$regime2[col]))
          MaxTice4$optim_ls[[i]]$regime2[[new_names]] <- exp(MaxTice4$optim_ls[[i]]$regime2[[col]])
     }
}

for(col in names(x)){
     new_names <- c(1:3, NA)
     new_names <- paste0("E", names(x[col]))
     x[[new_names]] <- exp(x[[col]])
}

head(x)

#won't work unless you add another var
optimGraphs1_all(MaxTice4, "logtice", "optLog2")
for(i in 1:length(MaxTice4$optim_ls)){
     df1 <- as.data.frame(MaxTice4$optim_ls[[i]]$df)
     df2 <- as.data.frame(MaxTice4$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(MaxTice4$optim_ls[[i]]$regime2)
     mm <- optimGraphs1(df1, df2, df3, yearLim, yearInt, lnbiomassInt,  titlenames[i], "Elogtice", "tice")
     ggsave(mm, filename = paste0("figs/covariates/optLog", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

     df1 <- as.data.frame(MaxTice4$optim_ls$`MaxTice-m1`$df)
     df2 <- as.data.frame(MaxTice4$optim_ls$`MaxTice-m2`$regime2)
     df3 <- as.data.frame(MaxTice4$optim_ls$`MaxTice-m1`$regime2)
     mm <- optimGraphs1(df1, df2, df3, yearLim, yearInt, lnbiomassInt,  "MaxTice-m1", "Elogtice", "tice")
     ggsave(mm, filename = paste0("figs/covariates/optLog", titlenames[i], ".pdf"), width=10, height=8, units="in")


graphics.off()
var1 <- "Elogtice"
var2 <- "tice"
p3 <- ggplot() + 
     geom_line(data = df3, aes_string(x = var1, y = "ExpectedLogBiomass"), colour="red", linetype=1, size=1.25) + 
     #geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomassOld"), colour="blue", linetype=1, size=1.25) +
     geom_point(data = subset(df1, year > 1991), aes_string(x = var2, y = "logcapelin"), shape=15, size=3) + 
     geom_errorbar(data = subset(df1, year > 1991), aes_string(x = var2, ymin="logcapelinlb", ymax="logcapelinub"), width = 0.3, colour = "black") +
     xlab(paste(var2)) +
     ylab("ln (Capelin biomass (ktons))") + 
     #ylim(0,9) +
     theme_bw()
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")


###
MaxTice5 <- calcFit_all1(cape_2001, titlenames, par = c(1, 200), var1 = "tice", var2 = "Ssurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = c(0:150),
                         x2_range = c(0:500))

MaxTice5$optim_ls$`MaxTice-m1`$cdf

optimGraphs1_all(MaxTice5, "tice", "opt_div")


### MaxTice6 log----
#DON'T CHANGE 
MaxTice6 <- calcFit_all1(cape_2001, titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50))
#c(250))
                              #
optimGraphs1_all(MaxTice6, "tice", "opt_div2")

x <- MaxTice6$optim_ls$`MaxTice-m1`
MaxTice6$optim_ls$`MaxTice-m1`$cdf
str(MaxTice6$optim_ls$`MaxTice-m1`$regime2)
head(MaxTice6$optim_ls$`MaxTice-m1`$regime2, 30)
range(x$df$Ssurface_tows_lag2)

MaxTice6$optim_ls$`MaxTice-m1`$df
a <- MaxTice6$optim_ls$`MaxTice-m1`$cdf$par[1]
b <- MaxTice6$optim_ls$`MaxTice-m1`$cdf$par[2]
c <- MaxTice6$optim_ls$`MaxTice-m1`$cdf$par[3]

a*80*(1-(80/b)) + c*118.1


###
MaxTice7 <- calcFit_all1(cape_2001, titlenames, 
                         par = c(1, 200, 1), 
                         var1 = "tice", var2 = "Ssurface_tows_lag2",
                         form1 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                         #form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = seq(0, 150, 10),
                         x2_range = seq(30, 500, 50))
#c(250))
#
optimGraphs1_all(MaxTice7, "tice", "opt_div3")
MaxTice7$optim_ls$`MaxTice-m1`$cdf

###
MaxTice3 <- calcFit_all1(cape_2001, titlenames, par = c(1, 20, 0.6), var1 = "Ntice", var2 = "Nsurface_tows_lag2",
                         form1 = "Alpha*tmp2*Beta",
                         form2 = "Alpha*tmp2*Beta*Gamma",
                         x1_range = c(-3:3),
                         x2_range = c(-3:3))
optimGraphs1_all(MaxTice3, "Nsurface_tows_lag2", "opt_norm_sur")

---------------------------------------------------------------
     ggplot() + 
     geom_line(data = MaxTice6$optim_ls$`MaxTice-m1`$regime2, aes(x = tice, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +  
     #geom_line(data = reg2, aes_string(x = var, y = "ExpectedLogBiomassOld"), colour="blue", linetype=1, size=1.25) +
     geom_line(data = subset(MaxTice6$optim_ls$`MaxTice-m1`$regime2, Ssurface_tows_lag2==230), aes(x = tice, y = ExpectedLogBiomass), colour="green", linetype=1, size=1.25) +
     geom_point(data = MaxTice6$optim_ls$`MaxTice-m1`$df,  aes(x = tice, y = logcapelin), shape=15, size=3) + 
     geom_errorbar(data = MaxTice6$optim_ls$`MaxTice-m1`$df, aes(x = tice, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
     #     xlab(paste(tic)) +
     ylab("ln (Capelin biomass (ktons))") + 
     #ylim(0,9) +
     theme_bw()
## IN DEVELOPMENT-----
#
# test of a generalization of the formula to multiple variables
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
titlenames <- c("test-m1", "test-m2", "test-m3", "test-m4", "test-m5", "test-m6")

test <- calcFit_all1(cape_2001, titlenames, par = c(5, 300, 5, 0.6), var1 = "Ntice", var2 = "Nsurface_tows_lag2",
                       form1 = "Alpha*tmp1*(1-(tmp1/Beta))",
                       form2 = "Alpha*tmp1*(1-(tmp1/Beta)) + Gamma*tmp2",
                       x1_range = c(-5:5), x2_range = c(-5:5))

range(cape_2001$capelin_m1$Ntice)
range(cape_2001$capelin_m1$Nsurface_tows_lag2)

test <- calcFit_all1(cape_2001, titlenames, par = c(5, 300, 5, 0.6), var1 = "Ntice", var2 = "Nsurface_tows_lag2",
                     form1 = "Delta*tmp2",
                     form2 = "Alpha*(tmp1*(1-(tmp1/Beta)) + Delta*tmp2)",
                     x1_range = c(-2:2), x2_range = c(-2:2))

Alpha <- 1
Beta <- 300
Delta <- 1
tmp1 <- cape_2001$capelin_m1$tice[2:15]
tmp2 <- cape_2001$capelin_m1$surface_tows_lag2[2:15]
Alpha*tmp1*(1-(tmp1/Beta)) + Delta*tmp2
View(cape_2001$capelin_m1)

yearInt <- seq(2000, 2017, by=2)
lnbiomassInt <- seq(-6, 10, by=2)
biomassInt <- seq(0, 8500)

head(test$optim_ls$`test-m1`$regime1)

for(i in 1:length(test$optim_ls)){
     df1 <- as.data.frame(test$optim_ls[[i]]$df)
     df2 <- as.data.frame(test$optim_ls[[i]]$regime1)
     df3 <- as.data.frame(test$optim_ls[[i]]$regime2)
     mm <- optimGraphs(df1, df2, df3, yearInt, lnbiomassInt,  titlenames[i], "Nsurface_tows_lag2")
     ggsave(mm, filename = paste0("figs/optimization/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}

head(test$optim_ls$`test-m1`$df)

x <- c(-20:20)
z <- c(-20:20)
Alpha <- 0.167
Beta <- 141
Gamma <- 0.5
y <- x + 3
y <- x^2
y <- -x^2
y <- x^-2
y <- Alpha*x - Gamma*z
y <- 1-x/Beta
y <- x*(1-x/Beta)
y <- Alpha*x*(1-x/Beta)
y <- Alpha*(x*(1-x/Beta) + Gamma*z)
y <- Gamma*z
y <- Gamma*z^2
y <- Gamma*z^3 + 3000
y <- Gamma*z^3 + Gamma*z^2 - Gamma
y <- Gamma*z^3 - Gamma*z^2 - Gamma
y <- Gamma*z^3 - Alpha*z^2 - Gamma + 30
plot(x, y)
plot(z, y)
plot(x, y, ylim=c(-20, 20))

cbind(z, y)


---------------------------------------------------------------
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

for(i in 1:seq_along(d5med_cape_2001)){
mm <- ggplot(data=d5med_cape_2001[[i]]) + geom_point(aes(x = surface_tows_lag2, y = logcapelin)) +
          geom_text(aes(x = surface_tows_lag2, y = logcapelin + 0.05, label=d5tice), size=3) + 
          geom_text(aes(x = surface_tows_lag2, y = logcapelin - 0.05, label=year), nudge_x = 100, size=3)
     ggsave(mm, filename = paste0("figs/covariates/scatter/", titlenames[i], ".pdf"), width=10, height=8, units="in")
}
