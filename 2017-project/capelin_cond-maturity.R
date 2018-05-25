#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-??-??, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to explore the influece of the percetnage of mature age2 Capelin and the percentage of age 2 capelin in the population to other variables. :
# 1) compare overall biomass and age-2 abundance - they are highly correlated (r-sq = 0.86)
# 2) see if there is a relationship between age-2 capelin abundance and the percentage in the pouplation that is Age-2 and percent mature Age2
# 3) see if there is a relationship between condition and perAge2 and age2PerMat
# 4) see if there is a relationship between tice and perAge2 and age2PerMat
rm(list=ls())

library(readr)
library(tidyr)
library(dplyr)

source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

####----
# load age and biomass datasets and merge
capelin_data_set <- "age" #
cap <- capelin_data(capelin_data_set)
cap <- filter(cap, year > 1998)

capelin_data_set <- "biomass" #
cap_b <- capelin_data(capelin_data_set)
cap_b <- filter(cap_b, year > 1998)

# join the data sets
cap_all <- left_join(cap, cap_b, by = "year")

# plot biomass v abundance
plot(cap_all$age2_log, cap_all$ln_biomass_med)
summary(lm(ln_biomass_med ~ age2_log, data = cap_all))

# plot abundance age-2 v percent Age2 in population
plot(cap$age2_log, cap$perAge2)

# relationship between percent Age2 in population and percent Age2 that are mature (rsq = 0.24)
# relationship between age2 percent mature and abundance of age2 - negative correlations (rsq = 0.74)
# relationship between percent Age2 that are mature and abundance of age2 (rsq = 0.14)
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(cap$age2PerMat, cap$perAge2)
plot(cap$age2PerMat, cap$age2_log10)
points(cap$age2PerMat[12], cap$age2_log10[12], col = "red", bg = "red", pch = 21) #show 2004 
points(cap$age2PerMat[16], cap$age2_log10[16], col = "blue", bg = "blue", pch = 21) #show 2014
plot(cap$perAge2, cap$age2_log10)
plot(cap$age2PerMat, cap_all$ln_biomass_med)
par(mfrow = c(1,1))

summary(lm(perAge2 ~ age2PerMat, data = cap))
summary(lm(age2_log10 ~ age2PerMat, data = cap))
summary(lm(age2_log10 ~ perAge2, data = cap))
summary(lm(ln_biomass_med ~ age2PerMat, data = cap_all))

####----
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")
cond_dat_a1 <- "cond" # 
cond1 <- condition_data(cond_dat_a1, 'data/condition_ag1_out.csv')
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
cond_dat_a1 <- "cond"
cond2 <- condition_data(cond_dat_a1_2, 'data/condition_ag2_out.csv')
cond2 <- filter(cond2, year > 1998)

dfc <- left_join(cond1, cond2, by = "year")
plot(dfc$meanCond.x, dfc$meanCond.y)

df2 <- left_join(cap, cond2, by = "year")
plot(df2$meanCond_lag, df2$age3)

