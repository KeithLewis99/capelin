
# load data sets through ice-capelin-jags_sequential
# a brief exploration of the Pseudocalanus data.  The big problem here is that there is one MASSIVE outlier.
# however, when this is eliminated, the results is not significant.  HM said other wise because they used robust regression....i'm not clear on how she arrived at this.  I have tested it against teh data back to 2001......

df3

temp <- df3[c("year", "ln_biomass_med", "ps_meanTot_lag2")]

# plot from 1999 to 2017
plot(temp$year, temp$ps_meanTot_lag2)
plot(temp$ps_meanTot_lag2, temp$ln_biomass_med)
abline(5.3576, 0.4831) # not sure where I got this

# remove 2010
plot(temp$ps_meanTot_lag2[-12], temp$ln_biomass_med[-12])
abline(5.4888, 0.9215) # not sure where I got this
summary(lm(ln_biomass_med ~ ps_meanTot_lag2, data = temp[-12, ]))

# back to 2001
temp <- df1[c("year", "ln_biomass_med", "ps_meanTot_lag2")]
plot(temp$year, temp$ps_meanTot_lag2)

summary(lm(ln_biomass_med ~ ps_meanTot_lag2, data = temp))
plot(temp$ps_meanTot_lag2, temp$ln_biomass_med)
abline(5.2806, 0.4945) # got this from summary two lines above

#remove the massive outlier
summary(lm(ln_biomass_med ~ ps_meanTot_lag2, data = temp[-16, ]))
plot(temp$ps_meanTot_lag2[-16], temp$ln_biomass_med[-16])
abline(5.3796, 0.8271) # got this from summary two lines above

