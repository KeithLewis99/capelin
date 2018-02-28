# http://doingbayesiandataanalysis.blogspot.ca/2012/01/parameterizing-gamma-distribution-by.html

# Specify desired mode and sd of gamma distribution:
mode = 1.86
sd = 1

# Here are the corresponding rate and shape parameter values:
ra = ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra 

show(sh)
show(ra)

# Graph it:
x = seq(0,mode+5*sd,len=1001)
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" , 
      main=paste("dgamma, mode=",mode,", sd=",sd,sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )



x <- seq(1, 100, 1)
y <- log(x)
plot(x, y)


x <- seq(100, 5000, 10)
y <- log(x)
plot(x, y)