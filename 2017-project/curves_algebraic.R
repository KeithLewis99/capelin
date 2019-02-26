x <- c(-100:100)
alpha <- 2 # increaes alpha, decrease |y|, negative turns parabola upside down
beta <- 100 # increase beta, decrease  |y|
gamma <- 2 # decrease gamma, decrease y
delta <- 1
z <- rnorm(201, mean=50, sd=10)
w <- rnorm(201, mean=30, sd=20)
hist(z)


y1 <- alpha*x
y2 <- 1-x/beta
y3 <- alpha*x/(1-x/beta)
y4 <- alpha*x*(1-x/beta)
y5 <- y4 + gamma*z
y6 <- y5 + delta*w

df <- cbind(y5, y4, y3, x, z, w)
head(df)

plot(x, y1)
plot(x, y2)
plot(x, y3)
plot(x, y4)
lines(x, y4)

plot(x, y5)

lines(x, y5)
plot(x, y6)
lines(x, y6)


alp <- MaxTice4$optim_ls$`MaxTice-m1`$cdf$par[1]
bet <- MaxTice4$optim_ls$`MaxTice-m1`$cdf$par[2]

x <- seq(0, 10, 0.1)
y <- alp*x*(1-x/bet)
plot(x, y)


lx <- 1:1000

yln <- log(lx)
ylg <- log10(lx)
ylg2 <- log2(lx)

plot(lx, yln)



# GOTO ice-capelin-covariates and run up to the MaxTice2 option - then, reproduce the below
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


#attempt at back transforming
# GOTO ice-capelin-covariates and run up to the MaxTice4 option - then, reproduce the below
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

# from Bolker
x = 1:20
a = 2
b = 1

y_det = a + b *x
y = rnorm(20, mean = y_det, sd = 2)
plot(x, y)
lines(x, y)



# Bolker pg 173
data(ReedfromgPred)
k <- 30
binomNLL1 <- function(p, k, N) {
     -sum(dbinom(k, prob = p, size = N, log = TRUE))
}

opt1 <- optim(fn = binomNLL1, par = c(p=0.5), N = 10, k = k, method = "BFGS")
