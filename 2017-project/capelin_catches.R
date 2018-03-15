
library(readr)
library(tidyr)
library(dplyr)
library(psych)
library(ggplot2)

rm(list=ls())
df <- read_csv('data/internationalcatches_1960-2017.csv')
glimpse(df)
df$total_2j3kl_l2 <- lag(df$total_2j3kl, 2)

plot(df$year, df$total_2j3kl)

df1 <- filter(df, year > 1984)
glimpse(df1)

cap <- read_csv('data/capelin-2017.csv')
glimpse(cap)
#View(cap)
cap$ln_abun_med <- log(cap$abundance_med)
cap$ln_ab_lci <- log(cap$ab_lci)
cap$ln_ab_uci <- log(cap$ab_uci)
cap$ln_biomass_med <- log(cap$biomass_med)
cap$ln_bm_lci <- log(cap$bm_lci)
cap$ln_bm_uci <- log(cap$bm_uci)

df2 <- left_join(df1, cap, by = "year")

plot(df2$year, df2$total_2j3kl)

plot(df2$total_2j3kl, df2$ln_biomass_med)
summary(lm(ln_biomass_med ~ total_2j3kl, data = df2))

plot(df2$total_2j3kl, df2$abundance_med)
summary(lm(abundance_med ~ total_2j3kl, data = df2))

View(df2)
plot(df2$total_2j3kl_l2, df2$abundance_med)
summary(lm(abundance_med ~ total_2j3kl_l2, data = df2))

plot(df2$total_2j3kl_l2, df2$ln_biomass_med)
summary(lm(ln_biomass_med ~ total_2j3kl_l2, data = df2))

df3 <- df2[c("year", "ln_biomass_med", "abundance_med", "total_2j3kl", "total_2j3kl_l2")]
View(df3)

df4 <- filter(df3, year > 1992)
plot(df4$total_2j3kl_l2, df4$ln_biomass_med)
summary(lm(ln_biomass_med ~ total_2j3kl_l2, data = df4))
