## ADB - June 28, 2018
## produce 3d plot for capelin SPERA paper


# load libraries ----------------------------------------------------------
library(lattice)

folder <- "figs"
# read data ---------------------------------------------------------------
cap <- read.csv('figs/vectors.csv', header = T, as.is = T)
range(cap$surface_tows_lag2, na.rm = T)


# produce plot ------------------------------------------------------------
## parameters
alpha <- -1.05336389345007
beta <- 0.634902231628861
gamma <- 17.1053505086351
delta <- 1.67715261064784
epsilon <- 0.273220043571445

## Explanatory variables
TI <- c(0, seq(0.1, 1.6, by = 0.1), delta)
ST <- seq(min(cap$surface_tows_lag2, na.rm = T), max(cap$surface_tows_lag2, na.rm = T), by = 0.12)
CO <- c(min(cap$meanCond_lag, na.rm = T), 
        mean(cap$meanCond_lag, na.rm = T), 
        max(cap$meanCond_lag, na.rm = T) )
plotdata <- expand.grid(TI, ST, CO)
names(plotdata) <- c('TI', 'ST', 'CO')

## predicted values
plotdata$fit <- alpha + beta * plotdata$ST + gamma * plotdata$TI * (1 - plotdata$TI/delta) + epsilon * plotdata$CO


## labels
tice_label <- expression(paste(t[italic(ice)]))
larval_label <- 'Larval abundance'
zlabel <- 'ln(capelin biomass (ktonnes))'

## standardize variables to make it easier to try different visualizations
colnames(plotdata)[colnames(plotdata) == "ST"] <- "y"
colnames(plotdata)[colnames(plotdata) == "TI"] <- "x"

## axis rotations
zrot <- 25
xrot <-  -85

## stanrdadize predicted values between 0 and1
range01 <- function(x){(x - min(x))/(max(x) - min(x))}
plotdata$fit01 <- range01(plotdata$fit)

## colors for the 3 planes
greycs <- rev(grey.colors(4))[1:3]

 predictedcap <-    wireframe(fit/max(fit01) ~  x * y , 
                              group = CO,
                              col.groups = greycs,
                              data = plotdata, 
                              xlab = list(tice_label, rot = 90 + xrot, cex = 1), 
                              ylab = list(larval_label,rot =  270 + zrot + 40, cex = 1),
                              zlab = list(zlabel, rot = 90),
                              scales = list(col = "black") ,
                              par.settings = list(axis.line = list(col = "transparent"))
                             ,zoom = 0.85
                             ,screen = list(z = zrot,  x = xrot)
                             ,distance = 0.1
)
 predictedcap


# save plot ---------------------------------------------------------------
jpeg(filename = 'predicted_CSM3a.jpg', res = 400, width = 177, height = 177, units = 'mm', quality = 100)
 print(predictedcap)
dev.off()

 