################################################################
#  Script written by Alejandro???? (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# Plot of the optimization curves in comparing capelin abundance with timing of sea ice retreat.
# This does an outstanding job except for 2013 and 2014 - why?  What covariates might explain the difference?

setwd("D:/Keith/ice/2017-project")
source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v2.R")
library(plotrix)



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
regime1 <- xtice[which(xtice$year == 1990),]
regime2 <- xtice[which(xtice$year == 2000),]

# make optimization graphs by year and in comparison to ice
pdf('ice-capelin-update-2014.pdf',height=9,width=9,pointsize =8)
par(mfrow=c(2,2))
# plot of capelin biomass v. year with ice models  
# plot data, CI, optimzation curves (red = 2014, blue = 2010), not sure what last line if for
       with(capelin,plot(x=year,y=logcapelin,lty=1,pch=16,xlab='Year',ylab='ln (Capelin biomass (ktons))',xlim=c(1982,2014),xaxp=c(1982,2014,4),ylim=c(2,9),yaxp=c(2,9,2),col='black')) 
        with(capelin,plotCI(x=year,y=logcapelin,uiw=logcapelinub-logcapelin,liw=logcapelin-logcapelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=16,xlab='',xaxt='n',ylab='',xlim=c(1982,2014),ylim=c(2,9)))
           with(subset(capelin,year>1981),lines(year,ExpectedLogBiomass,col='red', lwd=2,type='l'))
           with(subset(capelin,year>1981),lines(year,ExpectedLogBiomassOld,col='blue', lwd=2,type='l'))
           with(capelin,points(x=year,y=logcapelin,lty=1,pch=16,xlab='Year',ylab='ln (Capelin biomass (ktons))',xlim=c(1982,2014),xaxp=c(1982,2014,4),ylim=c(2,9),yaxp=c(2,9,2),col='black')) # think that this is redundant

           
           # plot of capelin biomass v. year with ice models  red = 2014, blue = 2010
       with(capelin,plot(x=year,y=capelin,lty=1,pch=16,xlab='Year',ylab='Capelin biomass (ktons)',xlim=c(1982,2014),xaxp=c(1982,2014,4),ylim=c(0,8500),yaxp=c(0,8500,2),col='black')) 
        with(capelin,plotCI(x=year,y=capelin,uiw=capelinub-capelin,liw=capelin-capelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=16,xlab='',xaxt='n',ylab='',xlim=c(1982,2014),ylim=c(0,8500)))       
           with(subset(capelin,year>1981),lines(year,exp(ExpectedLogBiomass),col='red', lwd=2,type='l'))
           with(subset(capelin,year>1981),lines(year,exp(ExpectedLogBiomassOld),col='blue', lwd=2,type='l'))
           legend(1998,8500,col=c('red','blue'),pch=NA,lwd=2,bty='n',cex=0.85,legend=c('Model estimates including data up to 2014', 'Model estimates including data up to 2010'))
           with(capelin,points(x=year,y=capelin,lty=1,pch=16,xlab='Year',ylab='Capelin biomass (ktons)',xlim=c(1982,2014),xaxp=c(1982,2014,4),ylim=c(0,8500),yaxp=c(0,8500,2),col='black')) 

# optimization curves           
  plot(regime1$tice[which(regime1$ExpectedLogBiomass>0)],regime1$ExpectedLogBiomass[which(regime1$ExpectedLogBiomass>0)],col='red', type='l', lwd=2,ylim=c(0,9),yaxp=c(0,9,2),xlab=labtice, ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
  
  lines(regime2$tice[which(regime2$ExpectedLogBiomass>0)], regime2$ExpectedLogBiomass[which(regime2$ExpectedLogBiomass>0)],col='red', type='l', lwd=2)
  
  lines(regime1$tice[which(regime1$ExpectedLogBiomassOld>0)], regime1$ExpectedLogBiomassOld[which(regime1$ExpectedLogBiomassOld>0)], col='blue', type='l', lwd=2,ylim=c(0,9),yaxp=c(0,9,2), xlab=labtice,ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
  
  lines(regime2$tice[which(regime2$ExpectedLogBiomassOld>0)],regime2$ExpectedLogBiomassOld[which(regime2$ExpectedLogBiomassOld>0)],col='blue', type='l', lwd=2,ylim=c(0,9),yaxp=c(0,9,2),xlab=labtice,ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
  axis(1, at = c(0,50,70,120,140), labels = TRUE, tick = TRUE, line = NA,pos = NA, outer = FALSE, font = NA, lty = "solid", lwd = 1, col = NULL)

     points(capelin$tice[which(capelin$year<1991)],capelin$logcapelin[which(capelin$year<1991)],type='p',ylim=c(0,9),xlim=c(0,200),pch= 2)
    points(capelin$tice[which(capelin$year>1990)],capelin$logcapelin[which(capelin$year>1990)],type='p',ylim=c(0,9),xlim=c(0,200),pch= 15)
     with(capelin[which(capelin$year<1991),],plotCI(x=tice,y=logcapelin,uiw=logcapelinub-logcapelin,liw=logcapelin-logcapelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=2,xlab='',xaxt='n',ylab='',xlim=c(0,200),ylim=c(0,9)))           
     with(capelin[which(capelin$year>1990),],plotCI(x=tice,y=logcapelin,uiw=logcapelinub-logcapelin,liw=logcapelin-logcapelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=15,xlab='',xaxt='n',ylab='',xlim=c(0,200),ylim=c(0,9)))           
#dev.off()

# other optimization curves
 # png(filename = "ice-capelin-update.png",width = 1000, height = 1000, units = "px", pointsize = 20, bg = "white", res = NA, family = "", restoreConsole = TRUE, type = c("cairo-png"))
  plot(regime1$tice[which(regime1$ExpectedLogBiomass>0)],regime1$ExpectedLogBiomass[which(regime1$ExpectedLogBiomass>0)],col='red', type='l', lwd=3,ylim=c(0,9),yaxp=c(0,9,2),xlab=labtice,ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
   lines(regime2$tice[which(regime2$ExpectedLogBiomass>0)],regime2$ExpectedLogBiomass[which(regime2$ExpectedLogBiomass>0)],col='red', type='l', lwd=3)
   lines(regime1$tice[which(regime1$ExpectedLogBiomassOld>0)],regime1$ExpectedLogBiomassOld[which(regime1$ExpectedLogBiomassOld>0)],col='blue', type='l', lwd=3,ylim=c(0,9),yaxp=c(0,9,2),xlab=labtice,ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
  lines(regime2$tice[which(regime2$ExpectedLogBiomassOld>0)],regime2$ExpectedLogBiomassOld[which(regime2$ExpectedLogBiomassOld>0)],col='blue', type='l', lwd=3,ylim=c(0,9),yaxp=c(0,9,2),xlab=labtice,ylab='ln (Capelin biomass (ktons))', xlim=c(0,190),xaxp=c(0,190,2))
  axis(1, at = c(0,50,70,120,140), labels = TRUE, tick = TRUE, line = NA,pos = NA, outer = FALSE, font = NA, lty = "solid", lwd = 1, col = NULL)
   points(capelin$tice[which(capelin$year<1991)],capelin$logcapelin[which(capelin$year<1991)],type='p',ylim=c(0,9),xlim=c(0,200),pch= 2)
    points(capelin$tice[which(capelin$year>1990)],capelin$logcapelin[which(capelin$year>1990)],type='p',ylim=c(0,9),xlim=c(0,200),pch= 15)
     with(capelin[which(capelin$year<1991),],plotCI(x=tice,y=logcapelin,uiw=logcapelinub-logcapelin,liw=logcapelin-logcapelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=2,xlab='',xaxt='n',ylab='',xlim=c(0,200),ylim=c(0,9)))           
     with(capelin[which(capelin$year>1990),],plotCI(x=tice,y=logcapelin,uiw=logcapelinub-logcapelin,liw=logcapelin-logcapelinlb,sfrac=0,add=T,lty=1, yaxt='n',gap=0,pch=15,xlab='',xaxt='n',ylab='',xlim=c(0,200),ylim=c(0,9)))           
dev.off()      
#####################################################################################################
#graphics.off()
#####################################################################################################
par(mfrow=c(1,1))

capelin$logdiff <- NA
for(i in 2:nrow(capelin)){
  capelin$logdiff[i] <- (capelin$logcapelin[i]-capelin$logcapelin[i-1])*100
}
 
#plot of logdiff in capelin from previous year    
png(filename = "ice-capelin-update2.png",width = 1000, height = 1000, units = "px", pointsize = 20, bg = "white", res = NA, 
family = "", restoreConsole = TRUE, type = c("cairo-png")) 
   with(capelin, plot(year,logdiff,xlab='year',ylab='L%',pch=16))
   abline(h=100,lty=3)
   abline(h=0)
   abline(h=-100,lty=3)
   
dev.off()   