library(readr)
library(tidyr)
library(dplyr)
library(psych)
library(ggplot2)

#just makes a graph of ln capelin biomass by year

rm(list=ls())

## read in source code-----
source('D:/Keith/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('D:/Keith/R/zuur_rcode/HighstatLibV7.R')
source("D:/Keith/capelin/2017-project/ice-capelin-covariates-FUN.R")
source("D:/Keith/capelin/2017-project/ice-capelin-jags_sequential-FUN.R")

capelin_data_set <- "biomass" 
cap <- capelin_data(capelin_data_set)

df1 <- filter(cap, year >1998)
p <- ggplot()  
p <- p + geom_point(data = cap, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 2)
p <- p + geom_errorbar(data = cap, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))

p <- p + xlab("Year") + ylab("ln(capelin biomass(ktons))")
p <- p + theme_bw(base_size = 30) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
p
ggsave("capelin_biomass.pdf", width=10, height=8, units="in")



p <- ggplot()  
p <- p + geom_point(data = df1, 
                    aes(y = ln_biomass_med, x = year),
                    shape = 16, 
                    size = 1.5)
p <- p + geom_errorbar(data = df1, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))

p <- p + xlab("Year") + ylab("ln(capelin biomass(ktons))")
p <- p + theme_bw(base_size = 30) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
p
ggsave("capelin_biomass_1998_2017.pdf", width=10, height=8, units="in")
