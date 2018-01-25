#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-01-25, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to:
# 1) Calculate a condition index for capelin: relative condition = observed/predicted
# 2) turn this into a markdown document


## libraries------
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)

library(purrr)
library(plotly)
library(plotrix)

rm(list=ls())

## load data----
# this from the ice-capelin-covariates file
df <- read_csv('data/capelin_condition_maturation_v1.csv')
glimpse(df)
head(df)

#create one filter for ice-capelin project and one for Fran and the markdown doc
df1 <- df %>%
     filter(year > 1992 & sex == 1 & age == 1 & maturity !=6) %>% #one-year males after 1992
     filter(!is.na(weight)) %>%
     filter(!is.na(length)) 

glimpse(df1)
levels(as.factor(df1$project))
levels(as.factor(df1$sample_number))
levels(as.factor(df1$year))
levels(as.factor(df1$month))
levels(as.factor(df1$nafo_div))
levels(as.factor(df1$maturity))
range(df1$weight)
range(df1$length)
levels(as.factor(df1$age))

head(df1$age, 100)
## DATA EXPLORATION - ZUUR 2010
## Step 1 Are there outliers in X and Y?


ggplot(df1, aes(x = weight)) + geom_dotplot() + facet_wrap(~year)

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(project), y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(year), y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(year), y = weight))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(month), y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(nafo_div), y = length))

ggplot(data=df1) +
     geom_boxplot(aes(x = as.factor(maturity), y = length))

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

low_vals <- filter(df1, weight < 1)

# Conc: still don't get the normality thing.  There are 0 values < 1, 1090 < 5: see Fran for what to filter

##Step 5 Is their collinearity among covariates
##Step 6 What are the relationships between Y and X variables

ggplot(data=df1) + geom_point(aes(x=length, y = weight))

ggplot(data=df1) + geom_point(aes(x=log10(length), y = log10(weight)))

#Conclusion: collinearity not an issue unless we divide by some factor.  Exponential relationship between weight and length but log-log takes care of that.

##Step 7 Are there interactions?  None unless we add terms


## Step 8 Are the observations of the response variable independent?
### No reason to expect unless there is observer error

## LINEAR MODELLING
m1 <- lm(log10(weight) ~ log10(length), data= df1)
summary(m1)
op <- par(mfrow = c(2,2))
plot(m1, add.smooth=F)
par(op)
par(mfrow = c(1,1))
# problems with residuals but this may be due to outliers
# 3 outliers; homogeneity looks good other than than.  Normality stinks but probably not a problem - very few values.

glimpse(df1)
df1$fits <- fitted(m1)
df1$rel.cond <- df1$weight/df1$fits
p <- ggplot(data=df1) +
     geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
     facet_wrap(~nafo_div)
p
