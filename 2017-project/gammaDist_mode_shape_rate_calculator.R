# http://doingbayesiandataanalysis.blogspot.ca/2012/01/parameterizing-gamma-distribution-by.html

# Specify desired mode and sd of gamma distribution:
mode = 1.86 # this is Ale's value from 2014 paper for Beta (width of dome)
sd = .6 # this is choosen arbitrarily but signifies two months (60days/100) for the width which means a month for the peak (Beta/2) which allows for reasonable values of the peak.

# Here are the corresponding rate and shape parameter values:
ra = ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra 

show(sh)
show(ra)
show(1/ra)

# Graph it:
x = seq(0,mode+5*sd,len=1001)
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" , 
      main=paste("dgamma, mode=",mode,", sd=",sd,sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )

plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" , 
      main=paste("dgamma, mode=", round(mode, 2),", \n, shape = " , round(sh, 2),", rate=", round(ra, 2), sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )

plot( x , dgamma( x , shape=2.8 , rate=1 ) , type="l" , 
      main=paste("dgamma, mode=", round(mode, 2),", \n, shape = " , round(sh, 2),", rate=", round(ra, 2), sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )

x <- seq(1, 100, 1)
y <- log(x)
plot(x, y)


x <- seq(100, 5000, 10)
y <- log(x)
plot(x, y)


estimate_mode <- function(x) {
     d <- density(x)
     d$x[which.max(d$y)]
}


mode = estimate_mode(ice1$tice*2/100)
mean(ice1$tice*2/100)
sd = sd(ice1$tice*2/100)
