# http://doingbayesiandataanalysis.blogspot.ca/2012/01/parameterizing-gamma-distribution-by.html

# Specify desired mode and sd of gamma distribution:
mode = 1.74 # this is Ale's value from 2014 paper for Beta (width of dome)
sd = .6 # this is choosen arbitrarily but signifies two months (60days/100) for the width which means a month for the peak (Beta/2) which allows for reasonable values of the peak.

# Here are the corresponding rate and shape parameter values:
ra = ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra 

show(sh)
show(ra)
show(1/ra)

# Graph it:
x = seq(0,mode+5*sd,len=1001)
# mode and sd
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" , 
      main=paste("dgamma, mode=",mode,", sd=",sd,sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )

#mode, shape, and rate
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" , 
      main=paste("dgamma, mode=", round(mode, 2),", \n, shape = " , round(sh, 2),", rate=", round(ra, 2), sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )


# an alternate way to estimate mode
estimate_mode <- function(x) {
     d <- density(x)
     d$x[which.max(d$y)]
}


mode = estimate_mode(ice1$tice*2/100)
mean(ice1$tice*2/100)
sd = sd(ice1$tice*2/100)


# the following is what Eric Pedersen did for us 2018-06-11
# We first considered the gamma and then the beta and lognormal distributions before coming back to the gamma.
# Eric created the below simulation to see how the priors actually perform which is an important step that I had overlooked.


alpha = rgamma(100, shape = 5,rate= 1/3) # randomly produce vector of length 100
beta = rgamma(100, shape = sh, rate=ra)
out = matrix(0, nrow = 500,ncol = 100) # produces 500x100 matrix
# produce values for tice
x = seq(0, 2, length.out = 500)
#loop through and make 100 values for each of the 500 x values
for(i in 1:100){
     out[,i] = alpha[i]*x*(1-x/beta[i])
}

# plot the columns of one matrix against the colums of another
matplot(x, out, type="l", lty=1,col="black")
matplot(x, out, type="l", lty=1)

# 
mean= 1.2
eta = (mean/3.65)*(1-mean/3.65)/(sd/3.65)^2 - 1

#alt_beta = rbeta(100,shape1 = mean/3.65*eta, 
#                 shape2 = (1-mean/3.65)*eta)*3.65

alt_beta = rnorm(100, mode, sd); 
alt_beta = ifelse(alt_beta<0, 0, alt_beta)
alt_beta = ifelse(alt_beta>3.65, 3.65, alt_beta)
