################################################################
#  Script written by Alejandro???? (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# The purpose of this file is to:
#1) Plot the optimization curves in comparing capelin abundance with timing of sea ice retreat.

# Notes:
# This does an outstanding job except for 2013 and 2014 - why?  What covariates might explain the difference?

setwd("D:/Keith/capelin/2017-project")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v3.R")
library(plotrix)
library(ggplot2)


capelin <- read.csv('capelin-ice-2014.csv',header=T)

capelin$myear <- capelin$year
capelin$myear[45:46] <- c(1960,1961)

## Optimization - produces lists of the optimization curve for use in figures below
# up to present
CapelinDomeFit <- optim(par = c(1,200,0.6),
                        dataf = capelin[which(capelin$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome, method=c("BFGS"))
# before 2010
CapelinDomeFitOld <- optim(par=c(0.2,180,0.6),
                        dataf = capelin[which(capelin$year < 2011 & capelin$logcapelin!='NA'),
                        c('year','tice','logcapelin')],
                        fn = SSQCapelinDome, method=c("BFGS"))

## Obtain Expected Log Capelin Biomass using parameters estimated in lines above
capelin$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par), dataf = capelin[,c('year','tice')])
capelin$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par), dataf = capelin[,c('year','tice')])


labtice <- expression(paste(italic(t[ice]), ' (day of year)')) # label for figures

# attach the optimization curves of capelin abundance to ice data
xtice <- expand.grid(year = c(1990,2000),tice = c(0:190,173.515,187.768))
xtice <- xtice[order(xtice$tice),]
xtice$ExpectedLogBiomassOld <- CapelinDome(params = c(CapelinDomeFitOld$par),dataf = xtice)
xtice$ExpectedLogBiomass <- CapelinDome(params = c(CapelinDomeFit$par),dataf = xtice)
#not sure what these are for but used in plots below but creates a data set where all values of year are the same????




# make optimization graphs by year and in comparison to ice
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for

# set values
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

yearInt <- seq(1982, 2014, by=4)
lnbiomassInt <- seq(0, 10, by=2)
biomassInt <- seq(0, 8500)

p1 <- ggplot(capelin, aes(x = year, y = logcapelin)) + 
  geom_errorbar(width = 0.3, colour = "black", aes(ymin=logcapelinlb, ymax=logcapelinub)) + 
  geom_point(shape=16, size=3)  +
  geom_line(aes(y=ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
  geom_line(aes(y=ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
  scale_y_continuous(limits = c(0,10), breaks = lnbiomassInt) +
  scale_x_continuous(limits = c(1982,2014), breaks = yearInt) +
  xlab('Year') +
  ylab('ln (Capelin biomass (ktons))') + 
  theme_bw()
p1

p2 <- ggplot(capelin, aes(x = year, y = capelin)) +
  geom_errorbar(width = 0.3, colour = "black", aes(ymin=capelinlb, ymax=capelinub)) + 
  geom_point(shape=16, size=3)  +
  geom_line(aes(y=exp(ExpectedLogBiomass)), colour="red", linetype=1, size=1.25) +
  geom_line(aes(y=exp(ExpectedLogBiomassOld)), colour="blue", linetype=1, size=1.25) +
  #scale_y_continuous(limits = c(0,8500), breaks = biomassInt) +
  scale_x_continuous(limits = c(1982,2014), breaks = yearInt) +
  xlab('Year') +
  ylab('Capelin biomass (ktons)') + 
  theme_bw() +
  annotate("text", x = 2008, y = 8400, label = "Model estimates to 2014") + # all following for the legend
  annotate("text", x = 2008, y = 8000, label = "Model estimates to 2010") +
  annotate("segment", x = 2001, xend = 2003, y = 8400, yend = 8400, colour = "red") +
  annotate("segment", x = 2001, xend = 2003, y = 8000, yend = 8000, colour = "blue")
p2


p3 <- ggplot() +
  geom_line(data = regime1, aes(x = tice, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
  geom_line(data = regime2, aes(x = tice, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
  geom_line(data = regime1, aes(x = tice, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
  geom_line(data = regime2, aes(x = tice, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
  geom_point(data = subset(capelin, year < 1991), aes(x = tice, y = logcapelin), shape=2, size=3) +
  geom_point(data = subset(capelin, year > 1991), aes(x = tice, y = logcapelin), shape=15, size=3) + 
  geom_errorbar(data = subset(capelin, year < 1991), aes(x = tice, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
  geom_errorbar(data = subset(capelin, year > 1991), aes(x = tice, ymin=logcapelinlb, ymax=logcapelinub), width = 0.3, colour = "black") +
  xlab("labtice") +
  ylab("ln (Capelin biomass (ktons))") + 
  ylim(0,9) +
  theme_bw()
p3

ggplot() +
  geom_line(data = regime1, aes(x = tice, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) + 
  geom_line(data = regime2, aes(x = tice, y = ExpectedLogBiomass), colour="red", linetype=1, size=1.25) +
  geom_line(data = regime1, aes(x = tice, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
  geom_line(data = regime2, aes(x = tice, y = ExpectedLogBiomassOld), colour="blue", linetype=1, size=1.25) +
  geom_point(data = subset(capelin, year < 1991), aes(x = tice, y = capelin), shape=2, size=3) +
  geom_point(data = subset(capelin, year > 1991), aes(x = tice, y = capelin), shape=15, size=3) + 
  geom_errorbar(data = subset(capelin, year < 1991), aes(x = tice, ymin=capelinlb, ymax=capelinub), width = 0.3, colour = "black") +
  geom_errorbar(data = subset(capelin, year > 1991), aes(x = tice, ymin=capelinlb, ymax=capelinub), width = 0.3, colour = "black") +
  xlab("labtice") +
  ylab("Capelin biomass (ktons)") + 
  ylim(0,9) +
  theme_bw()


# make multiplot
#pdf('ice-capelin-update-2014-new.pdf',height=9,width=9,pointsize =8)
multiplot(p1, p3, p2, cols=2)
#dev.off()


#####################################################################################################
#graphics.off()
#####################################################################################################
par(mfrow=c(1,1))

capelin$logdiff <- NA
for(i in 2:nrow(capelin)){
  capelin$logdiff[i] <- (capelin$logcapelin[i]-capelin$logcapelin[i-1])*100
}

# new version of figure below
#png(filename = "ice-capelin-update2.png",width = 1000, height = 1000, units = "px", pointsize = 20, bg = "white", res = NA, family = "", restoreConsole = TRUE, type = c("cairo-png")) 

ggplot(capelin, aes(x = year, y = logdiff)) + 
  geom_point(pch = 16, size = 3) + 
  geom_hline(yintercept = 100, lty = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -100, lty = 3) +
  xlab("Year") +
  ylab("L%") + 
  #ylim(0,9) +
  theme_bw()

#plot of logdiff in capelin from previous year    
   with(capelin, plot(year,logdiff,xlab='year',ylab='L%',pch=16))
   abline(h=100,lty=3)
   abline(h=0)
   abline(h=-100,lty=3)
   
#dev.off()   