#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-01-25, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to:
# 1) Calculate a condition index for capelin: relative condition = observed/predicted for the capelin project
# 2) turn this into a markdown document for Fran
# From "Fall acoustic and Campelen strat samples 1981-2015 from 2J3KL_v1.xlsx" from F. Mowbray
# Uses "capelin_condition_maturation_v1.csv" - copied and pasted from above
# generates "condition_out.csv - used in ice-capelin-jagps_sequential.R


## libraries------
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(magrittr)

rm(list=ls())

## load data----
# this from the ice-capelin-covariates file
df <- read_csv('data/capelin_condition_maturation_v1.csv')
glimpse(df)
head(df)

#create one filter for ice-capelin project and one for Fran and the markdown doc
df1 <- df %>%
     filter(year > 1992 &  sex==1 & age == 2 & maturity != 6 & project != 10 & as.factor(month) %in% c("10", "11", "12")) %>% #one-year males after 1992, just project 23 sex == 1 &
     filter(!is.na(weight)) %>%
     filter(!is.na(length)) 

cols <- c("project", "nafo_div", "project_code", "sex", "maturity")
df1 %<>%
     mutate_each_(funs(factor(.)),cols)
glimpse(df1)

levels(df1$nafo_div)[levels(df1$nafo_div) == "23"] <- "2J"
levels(df1$nafo_div)[levels(df1$nafo_div) == "31"] <- "3K"
levels(df1$nafo_div)[levels(df1$nafo_div) == "32"] <- "3L"
levels(df1$sex)[levels(df1$sex) == "1"] <- "Male"
levels(df1$sex)[levels(df1$sex) == "2"] <- "Female"

# look at levels of data
levels(df1$project)
levels(as.factor(df1$sample_number))
levels(as.factor(df1$year))
levels(as.factor(df1$month))
levels(df1$nafo_div)
levels(df1$maturity)
range(df1$weight)
range(df1$length)
levels(as.factor(df1$age))

# filter out 99% by year - how to do this????
summary(df1$weight)
mean_w <- mean(df1$weight)
sd_w <- sd(df1$weight)
nrow(df1)
uci <- mean_w + qnorm(0.995)*sd_w/sqrt(nrow(df1))
quantile(df1$weight, c(0.5, 0.95, 0.99, 0.995, 0.999))

count_n <- df1 %>%
     group_by(year) %>%
     summarize(n = n())

quant <- df1 %>%
     group_by(year) %>%
     summarize(uq = quantile(weight, c(0.99)), lq = quantile(length, c(0.01)))

# for age 2
df1 <- df1 %>%
     group_by(year) %>%
     right_join(x=df1, y=quant, by ="year") %>%
     filter(weight < uq & length > lq)
# for age 1
#df1 <- df1 %>%
 #    group_by(year) %>%
  #   right_join(x=df1, y=quant, by ="year") %>%
   #  filter(weight < uq)

## DATA EXPLORATION - ZUUR 2010----
## Step 1 Are there outliers in X and Y?

#dotplot
ggplot(df1, aes(x = length)) + geom_dotplot() + facet_wrap(~year)

ggplot(df1, aes(x = weight)) + geom_dotplot() + facet_wrap(~year)

# Fran - should we got with just project 23 for ice-capelin? yes, based on next analysis?

# project - L/W larger in 23
# outliers generally small
ggplot(data=df1) +
     geom_boxplot(aes(x = project, y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = project, y = weight))

# year
# shows some significant outliers in a few years.  See Fran
ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(year), y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(year), y = weight))

#month L/W increase with month bc growing 1cm/month
ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(month), y = length))
ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(month), y = weight))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(month), y = weight)) +      facet_wrap(~nafo_div)

# nafo_div (3L lower length and weight)- probably older bc of latter in season
ggplot(data=df1) +
     geom_boxplot(aes(x = nafo_div, y = length))
ggplot(data=df1) +
     geom_boxplot(aes(x = nafo_div, y = weight))

# maturity L/W higher for 2
ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(maturity), y = length))
ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(maturity), y = weight))

# Conclusion: 
# Differences between projects
# Length/weight low in early 1990s
# only a few values for maturity 6 - filtered earlier

##Step 3 Are the data normally distributed?
##Step 4 Are there lots of zeros in the data
p <- ggplot(df1, aes(x=weight))
p + geom_histogram() 

#or density
p <- ggplot(df1, aes(x=weight))
p + geom_density() 

filter(df1, weight < 5) # no lower limit

# Conc: these are skewed to the left: still don't get the normality thing.  There are 0 values < 1, 1090 < 5: see Fran for what to filter

##Step 5 Is their collinearity among covariates
##Step 6 What are the relationships between Y and X variables

ggplot(data=df1) + geom_point(aes(x=length, y = weight))

ggplot(data=df1) + geom_point(aes(x=log10(length), y = log10(weight), colour=nafo_div)) + 
     facet_wrap(~nafo_div)

#Conclusion: collinearity not an issue unless we divide by some factor.  Exponential relationship between weight and length but log-log takes care of that.

##Step 7 Are there interactions?  None unless we add terms


## Step 8 Are the observations of the response variable independent?
### No reason to expect unless there is observer error

## LINEAR MODELLING----
m1 <- lm(log10(weight) ~ log10(length), data= df1)
#m1 <- lm(log10(weight) ~ log10(length), weights = nafo_div, data= df1)
summary(m1)
str(m1)
str(summary(m1))

sum <- summary(m1)
rsq <- sum$r.squared
slope <- sum$coefficients[2,1]
ste <- sum$coefficients[2,2]

op <- par(mfrow = c(2,2))
plot(m1, add.smooth=F)
par(op)
par(mfrow = c(1,1))

# problems with residuals but this may be due to outliers
# 3 outliers; homogeneity looks good other than than.  Normality stinks but probably not a problem - very few values.
# 2405, 3962, 2227 - outliers before filtering for 99%
data.frame(df1[2405,])# short and heavy
data.frame(df1[3962,]) # long and light
data.frame(df1[2227,])# long and light

# new outliers
data.frame(df1[2006,])# short and heavy
data.frame(df1[1828,]) # long and light
data.frame(df1[3542,])# long and light
#leverage
data.frame(df1[332,])# short and light

# for age 2
data.frame(df1[2363,])

##output----
df1$fits <- fitted(m1)
m1$fitted.values
str(m1)
df1$rel.cond <- df1$weight/10^df1$fits

df1$resids <- df1$weight-10^df1$fits

plot(df1$rel.cond, df1$resids)
# produce table
df1 %>%
     group_by(year, nafo_div) %>%
     summarize(meanCond = round(mean(rel.cond),2), stdCond= round(sd(rel.cond),2)) %>% 
     unite(mean, meanCond:stdCond, sep = " +/- ") %>%
     spread(key = nafo_div, value = mean)

# produce output for Bayesain analysis
out <- df1 %>%
     group_by(year) %>%
     summarize(meanCond = round(mean(rel.cond),4), stdCond= round(sd(rel.cond),4), medCond = round(median(rel.cond), 4))
#write_csv(out, "data/condition_ag1_out.csv")
write_csv(out, "data/condition_ag2_out.csv")

#nafo
ggplot(data=df1) +
     geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
     facet_wrap(~nafo_div)
# there appear to be a few outliers here but they should have minimal influence on the final outcome
nrow(subset(df1, rel.cond > 1.5))
nrow(subset(df1, rel.cond < 0.75))

#month
ggplot(data=df1) +
     geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
     facet_wrap(~month)
#project
ggplot(data=df1) +
     geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
     facet_wrap(~project)

ggplot(data=df1) +
     geom_boxplot(aes(x = month, y = rel.cond, group = month)) + 
     facet_wrap(~nafo_div)

ggplot(data=df1) +
     geom_point(aes(x = log10(length), y = log10(weight), colour = nafo_div)) + 
     facet_wrap(~year)


## compare observed/expected to residuals - is there a difference
condResids <- read.csv('data/archive/condition_resids.csv')
str(condResids)
condResids[19:20,2] <- NA
condResids[c(19, 20), 1] <- matrix(c(2016, 2017), ncol = 1) 

str(condResids)
View(condResids)
condResids$resids_lag <- lag(condResids$resids, 1)

# compare cond/exp and resids
comp <- left_join(out, condResids, by = "year")
plot(comp$meanCond, comp$resids)
summary(lm(meanCond ~ resids, data = comp))