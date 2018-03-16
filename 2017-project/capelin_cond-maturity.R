rm(list=ls())

source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

####----
capelin_data_set <- "age" #
cap <- capelin_data(capelin_data_set)
cap <- filter(cap, year > 1998)

capelin_data_set <- "biomass" #
cap_b <- capelin_data(capelin_data_set)
cap_b <- filter(cap_b, year > 1998)

cap_all <- left_join(cap, cap_b, by = "year")
plot(cap_all$age2_log, cap_all$ln_biomass_med)
summary(lm(ln_biomass_med ~ age2_log, data = cap_all))
plot(cap$age2_log, cap$perAge2)

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(cap$age2PerMat, cap$perAge2)
plot(cap$age2PerMat, cap$age2_log10)
points(cap$age2PerMat[12], cap$age2_log10[12], col = "red", bg = "red", pch = 21) #show 2004 
points(cap$age2PerMat[16], cap$age2_log10[16], col = "blue", bg = "blue", pch = 21) 
plot(cap$perAge2, cap$age2_log10)
par(mfrow = c(1,1))

summary(lm(perAge2 ~ age2PerMat, data = cap))
summary(lm(cap$age2_log10 ~ age2PerMat, data = cap))
summary(lm(cap$age2_log10 ~ perAge2, data = cap))

####----
cond <- "cond" # 
cond1 <- condition_data(cond, 'data/condition_ag1_out.csv')
cond1 <- filter(cond1, year > 1998)

df <- left_join(cap_all, cond1, by = "year")
glimpse(df)
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(df$meanCond_lag, df$age2PerMat)
plot(df$meanCond_lag, df$perAge2)
plot(df$meanCond_lag, df$age2_log10)
par(mfrow = c(1,1))


ice <- read_csv('output-processing/capelin-m1.csv')
glimpse(ice)


df <- left_join(df, ice, by = "year")
View(df)
plot(df$tice, df$age2PerMat)
plot(df$tice, df$perAge2)
plot(df$tice, df$meandCond_lag)

df[12, c("year", "tice", "perAge2", "age2PerMat", "age2_log10", "ln_biomass_med", "meanCond_lag")]

####----
cond2 <- condition_data(cond, 'data/condition_ag2_out.csv')
cond2 <- filter(cond2, year > 1998)

dfc <- left_join(cond1, cond2, by = "year")
plot(dfc$meanCond.x, dfc$meanCond.y)

df2 <- left_join(cap, cond2, by = "year")
plot(df2$meanCond_lag, df2$age3)
