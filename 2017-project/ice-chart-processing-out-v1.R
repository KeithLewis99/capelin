################################################################
#  Script written by Paul Regular (Paul.Regular@dfo-mpo.gc.ca)  #
#  Created 201X-XX-XX, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Paul Regular, Alejandro Buren, and Keith Lewis 2017-07.04 #
################################################################

# The purpose of this file is to:
#1) 

## Area analysis ---------------------------------------------------------------

source("D:/Keith/capelin/2017-project/ice-chart-processing-function-v2.R")
library(data.table)

load("ice_trends_2017.Rdata") # this is date, area, and volume
trends <- data.table(trends) # shows head and tail of "trends"
trends$doy <- as.POSIXlt(trends$date)$yday
trends$doy <- ifelse(trends$doy > 305, trends$doy - 365, trends$doy) # force doy to negative if between 305 to 365 (nov-dec)
trends$day <- as.integer(format(trends$date, "%d")) 
trends$month <- as.integer(format(trends$date, "%m"))
trends$year <- as.integer(format(trends$date, "%Y"))
trends$year <- ifelse(trends$month %in% 11:12, trends$year + 1, trends$year) # if nov or dec, plus 1 to year
trends$psudo.date <- as.Date(ISOdate(ifelse(trends$month %in% 11:12, 2015, 2016), 
                                     trends$month, trends$day)) # for plotting
lookAt(trends)


############create a bunch of figures
# area v. date
ggplot() + geom_line(aes(x = date, y = area), trends, size = 1)

# graph in three year chunks
p1 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 1990:1992], size = 1)
p1
p2 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 1995:1997], size = 1)

p3 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 2006:2008], size = 1)

p4 <- ggplot() + geom_line(aes(x = psudo.date, y = area, col = factor(year)), 
                     trends[year %in% 2007:2010], size = 1)

multiplot(p1, p3, p2, p4, cols=2)

# facet graph over 11 years
ggplot() + geom_line(aes(x = psudo.date, y = area), 
                     trends[year %in% 2005:2016], size = 1) + facet_wrap(~ year)


## Not sure if volume is appropriate...too big of an approximation maybe?


years <- unique(trends$year)
adv.mag <- ret.mag <- amax <- tstart <- tend <- tmax <-
  raw.amax <- raw.tmax <- rep(NA, length(years))

# this prints out ~ 50 graphs plotting the predicted values of two differnt models as well as the data
# I think it is relating ice area and day of year - presumably ice retreat.

for(i in seq_along(years)) {
  
  ## subset data
  mod.dat <- trends[year == years[i]]
  mod.dat$t <- mod.dat$doy
  mod.dat$a <- mod.dat$area/1000
  mod.p <- ggplot(mod.dat) + geom_point(aes(x = t, y = a)) + 
    theme_bw() + theme(panel.grid.major = element_line(linetype = "dotted")) # make base layer and plot data
  
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
    mod.p <- mod.p + geom_line(aes(x = t, y = a, col = mod), preds, size = 0.7) # plots the predicted values 
    
  }
  
  ## print plot and clean-up
  print(mod.p + ggtitle(years[i])) # add he years
  rm(mod1, mod2)
  
}



pars <- data.table(year = years, adv.mag, ret.mag, amax, tstart, 
                   tend, tmax, raw.amax, raw.tmax)
pars$tdur <- pars$tend - pars$tstart
pars$ret.dur <- pars$tend - pars$tmax
pars$adv.dur <- pars$tmax - pars$tstart
pars$mag <- pars$adv.mag + pars$ret.mag

save(pars, file = "ice_pars.Rdata")


## Area
ggplot(pars, aes (x = year, y = amax)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Timing
tice <- melt(pars, id.vars = "year", measure.vars = c("tstart", "tmax", "tend"),
             variable.name = "timing", value.name = "doy")
ggplot(tice, aes(x = year, y = doy, col = timing)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Duration
dice <- melt(pars, id.vars = "year", measure.vars = c("tdur", "adv.dur", "ret.dur"),
             variable.name = "type", value.name = "duration")
ggplot(dice, aes(x = year, y = duration, col = type)) + 
  geom_line(size = 0.7) + geom_point(size = 2)

## Maginitude
mice <- melt(pars, id.vars = "year", measure.vars = c("mag", "adv.mag", "ret.mag"),
             variable.name = "type", value.name = "magnitude")
ggplot(mice, aes(x = year, y = magnitude, col = type)) + 
  geom_line(size = 0.7) + geom_point(size = 2)



# how to filter the ice by type?

