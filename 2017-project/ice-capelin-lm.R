#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-01-18, R version 3.3.3 (2017-03-06)  
# the purpose of this file is to explore relationships between capelin abundance/biomass and some explanatory variables using simple linear models.

# testing branching
# test again

library(psych)



rm(list=ls())
## load data----
df <- read.csv('figs/covariates/capelin_covariates_2001.csv',header=T)
str(test)

# some subsets of the larger dataset - useful for viewing
sdf1 <- df[c("year", "logcapelin", "age2_log10", "tice", "log10surface_tows_lag2", "ps_meanTot_lag2", "resids_adj")]
View(sdf1)
paul <- cape_2001$capelin_m1[c("year", "logcapelin", "surface_tows_lag1", "age1")]
plot(paul$surface_tows_lag1, paul$age1)

# relationships and correlations among RV and EV
pairs.panels(df[c("age2_log10", "surface_tows_lag2", "ps_meanTot_lag2", "resids_adj",  "tice")], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses,
             cex.labels = 2,
             cex.cor = 1
)

summary(lm(age2_log10 ~ log10surface_tows_lag2, data = sdf1))
m1 <- lm(age2_log10 ~ log10surface_tows_lag2, data = sdf1)
str(m1)
str(m1$resid)
op <- par(mfrow = c(2,2))
plot(m1, add.smooth=F)
par(op)
par(mfrow = c(1,1))

resids_t <- data.frame(resids = rep(NA, 14))
resids_t[1:3, ] <- m1$resid[1:3]
resids_t[5:13, ] <- m1$resid[4:12]

x <- cbind(sdf1, resids_t)
plot(x$resids, x$tice, na.rm=T)
plot(x$resids, x$logcapelin, na.rm=T)

##
summary(lm(age2_log10 ~ ps_meanLog_lag2, data = sdf))
layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(m1)
layout(matrix(1,1))
View(sdf)

sdf_outlier <- sdf[-12, ]
summary(lm(age2_log10 ~ ps_meanLog_lag2, data = sdf_outlier))
m1 <- lm(age2_log10 ~ ps_meanLog_lag2, data = sdf_outlier)
plot(sdf_outlier$ps_meanLog_lag2, sdf_outlier$age2_log10)
layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(m1)
layout(matrix(1,1))

summary(lm(age2_log10 ~ resids_adj, data = sdf))

summary(lm(age2_log10 ~ log10surface_tows_lag2 + ps_meanLog_lag2, data = sdf))

# age2_log10 ~ log10surface_tows_lag2 + resids_adj
m2 <- lm(age2_log10 ~ log10surface_tows_lag2 + resids_adj, data = df)
summary(m2)
op <- par(mfrow = c(2,2))
plot(m2, add.smooth=F)
par(op)
par(mfrow = c(1,1))

resids_t <- data.frame(resids = rep(NA, 14))
resids_t[1:3, ] <- m2$resid[1:3]
resids_t[5:13, ] <- m2$resid[4:12]

x <- cbind(sdf1, resids_t*100)
plot(x$tice, x$resids, na.rm=T)
plot(x$log10surface_tows_lag2, x$resids)
plot(x$resids_adj, x$resids)

library(nlme)
m2.lm <- gls(age2_log10 ~ log10surface_tows_lag2 + resids_adj, data = df, na.action = na.exclude)
vf1Fix <- varFixed(~log10surface_tows_lag2)
m2.gls <- gls(age2_log10 ~ log10surface_tows_lag2 + resids_adj, data = df, na.action = na.exclude, weights = vf1Fix)

anova(m2.lm, m2.gls)

summary(lm(age2_log10 ~ log10surface_tows_lag2*ps_meanLog_lag2, data = sdf))
summary(lm(age2_log10 ~ log10surface_tows_lag2 + ps_meanLog_lag2 + resids_adj, data = sdf))

# age2_log10 ~ log10surface_tows_lag2 + ps_meanLog_lag2 + resids_adj
m3 <- lm(age2_log10 ~ log10surface_tows_lag2 + ps_meanLog_lag2 + resids_adj, data = df)
summary(m3)


resids_t <- data.frame(resids = rep(NA, 14))
resids_t[1:3, ] <- m3$resid[1:3]
resids_t[5:13, ] <- m3$resid[4:12]

x <- cbind(sdf1, resids_t*100)
plot(x$tice, x$resids, na.rm=T)
View(x)

m2 <- glm(age2_log10 ~ log10surface_tows_lag2, data = sdf, family = poisson)
summary(m2)
layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(m2)
layout(matrix(1,1))

summary(lm(age2_log10 ~ resids_adj, data = sdf ))
