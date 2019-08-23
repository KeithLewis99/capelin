# this is my rewrite of Keith's relative condition index script for capelin. Doing this to better deal with some recent changes to the database and to enable better read-in of newer data in the future...

#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2018-01-25, R version 3.3.3 (2017-03-06)             #

# The purpose of this file is to:
# 1) Calculate a condition index for capelin: relative condition = observed/predicted for the capelin project
# 2) turn this into a markdown document for Fran
# From "Fall acoustic and Campelen strat samples 1981-2015 from 2J3KL_v1.xlsx" from F. Mowbray
# Uses "capelin_condition_maturation_v1.csv" - copied and pasted from above
# generates "condition_out.csv - 
# output used in ice-capelin-jagps_sequential.R and ice-caplein_jags_knock_out.R

# there has been a lot of debate about this variable.  What age to use and what sex to use.  Age 1 capelin seem to make the most sense as that is the condition of capelin going into the second year.

## libraries------
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(tidyselect)
library(magrittr)

rm(list=ls())

## load data----
# this from the ice-capelin-covariates file
df <- read_csv('data/capelin_condition_maturation_v1.csv') # read in full lsm data
df <- df[c(1:5, 7:15)]
glimpse(df)
head(df)

# pulling in data from a mix of sources...
# want to use as much data as possible from database and then the file just provided by Megan

# first is pulling in database data...
# folder with the database files
#dbasefolder<-"E:\\Database backups\\Active Mar 3 - 2019\\"
dbasefolder <- "D:\\Keith\\capelin\\2019-code\\data\\capelin_condition"
dbasefiles<-dir(dbasefolder)

# get the set of LSM files
lsmfiles<-grep("^LSM[0-9]{2}.txt$",dbasefiles)
lsm<-read.csv(paste(dbasefolder,dbasefiles[lsmfiles[1]],sep="//"), stringsAsFactors = FALSE)
for(i in 2:length(lsmfiles)){
  message(i)
  tempfile<-read.csv(paste(dbasefolder,dbasefiles[lsmfiles[i]],sep="//"), stringsAsFactors = FALSE)
  lsm<-rbind(lsm,tempfile[,colnames(tempfile) %in% colnames(lsm)])
}

# get the set of Master files
masterfiles<-grep("^Master[0-9]{2}.txt$",dbasefiles)
master<-read.csv(paste(dbasefolder,dbasefiles[masterfiles[1]],sep="//"), stringsAsFactors = FALSE)
for(i in 2:length(masterfiles)){
  message(i)
  tempfile<-read.csv(paste(dbasefolder,dbasefiles[masterfiles[i]],sep="//"), stringsAsFactors = FALSE)
  master<-rbind(master,tempfile[,colnames(tempfile) %in% colnames(master)])
}

# cleaning up master here to make life easier later...
master<-master[,c("Sample.Number", "Year", "Month", "Day", "ICNAF.Div", "Locality", "Area.square.", "Gear", "Country", "InShore.OffShore", "Ship", "Trip.Number", "Set.Number", "Strat.Protocol", "Project.Code")]

# get the set of strat files
stratfiles<-grep("^Strat[0-9]{2}.txt$",dbasefiles)
strat<-read.csv(paste(dbasefolder,dbasefiles[stratfiles[1]],sep="//"), stringsAsFactors = FALSE)
for(i in 2:length(stratfiles)){
  message(i)
  tempfile<-read.csv(paste(dbasefolder,dbasefiles[stratfiles[i]],sep="//"), stringsAsFactors = FALSE)
  if(!("Stomach_analysis_flag" %in% names(tempfile))){
    tempfile$Stomach_analysis_flag<-NA
  }
  strat<-rbind(strat,tempfile[,colnames(tempfile) %in% colnames(strat)])
}
# clean up strat a bit to make life easier later on...
strat<-strat[,c("Sample.Number", "Year", "Otolith.Number", "LSM.Number", "Length", "Sex", "Maturity", "Weight", "Stomach.Fullness", "Gonad", "Age", "Fecundity")]

# start combining data sets...
lsm.master<-merge(lsm,master,all.x=TRUE,all.y=TRUE)
lsm.master.strat<-merge(lsm.master,strat,all.x=TRUE)

# create a new lsm.master.strat archive lsm number
lsm.master.strat$new.tracking.lsm <- lsm.master.strat$Year*1000000 + 
  lsm.master.strat$Trip.Number*1000 +
  lsm.master.strat$Sample.Number +
  lsm.master.strat$LSM.Number/1000

# select the appropriate columns, and only keep project 23
lsm.master.strat2 <- lsm.master.strat %>% 
  select(Project.Code,Sample.Number,Year, Month, ICNAF.Div, new.tracking.lsm, Length, Sex, Maturity, Weight, Stomach.Fullness, Age, Fecundity, Gonad) %>% 
  filter(Project.Code==23) %>% 
  filter(!is.na(Weight))

colnames(lsm.master.strat2)<-c("project", "sample_number", "year", "month", "nafo_div", "archive_lsm_number", "length", "sex", "maturity", "weight", "stomach_fullness", "age", "fecundity", "gonad")

# I can't match older data, so I'll just use the values after the data that Keith has from the database...

lsm.master.strat3 <- lsm.master.strat2 %>% filter(year>2015)


## finally load data from Megan - can ignore this for now
#newest.data<-read.csv(file="E:\\Capelin\\2019 capelin assessment\\Condition\\all_eri_capelin_2018_cleaned.csv",header=TRUE)
newest.data<-read.csv(file="D:\\Keith\\capelin\\2019-code\\data\\all_eri_capelin_2018_cleaned.csv",header=TRUE)

head(newest.data)

# combine the two data sets
df <- rbind(df, lsm.master.strat3)
glimpse(df)
#View(tail(df, 100))

# save the combined data
write_csv(df, "E:\\Capelin\\2019 capelin assessment\\Condition\\condition_1979_2018.csv")
write_csv(df, "D:\\Keith\\capelin\\2019-code\\data\\condition_1979_2018.csv")
# 

##Age1 M/F comparisons----
# there is considerable worry that using the males only will cause questions during the review process.  The intent of the following analysis is to test for differences between males and females.

#subset data
df3 <- df %>%
  filter(year > 1992 &  age == 1 & maturity != 6 & project != 10 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(23, 31, 32)) %>% #one-year old capelin after 1992, just project 23 sex == 1 &
  filter(!is.na(weight)) %>%
  filter(!is.na(length)) 

#make these variables factors
cols <- c("project", "nafo_div", "sex", "maturity")
df3 %<>%
  mutate_each_(funs(factor(.)),cols)
glimpse(df3)

#change values to something interpretable
df3$nafo_div <- as.factor(df3$nafo_div)
levels(df3$nafo_div)[levels(df3$nafo_div) == "23"] <- "2J"
levels(df3$nafo_div)[levels(df3$nafo_div) == "31"] <- "3K"
levels(df3$nafo_div)[levels(df3$nafo_div) == "32"] <- "3L"
levels(df3$sex)[levels(df3$sex) == "1"] <- "Male"
levels(df3$sex)[levels(df3$sex) == "2"] <- "Female"

# check levels of data to ensure that subset worked
levels(df3$project)
levels(as.factor(df3$sample_number))
levels(as.factor(df3$year))
levels(as.factor(df3$month))
levels(df3$nafo_div)
levels(df3$maturity)
range(df3$weight)
range(df3$length)
levels(as.factor(df3$age))


# see email from H. Murphy and Methods of ResDoc - these are the values for filtering the dataset that Hannah and Maegen suggested based on biology
df3 <- df3 %>%
  group_by(year) %>%
  filter(weight > 2) %>% 
  filter(length > 80)

# run models for male and females seperately
#males
df3.m <- filter(df3, sex == "Male")
m2.m <- lm(log10(weight) ~ log10(length), data= df3.m)
df3.m$fits <- fitted(m2.m)

m2.m$fitted.values
str(m2.m)
df3.m$rel.cond <- df3.m$weight/10^df3.m$fits
df3.m$resids <- df3.m$weight-10^df3.m$fits


ggplot(data=df3.m) + geom_point(aes(x=log10(length), y = log10(weight), colour=nafo_div)) + 
  facet_wrap(~nafo_div)

# test for outliers of condition
quantile(df3.m$weight, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999))
plot(density(df3.m$weight))
round(quantile(df3.m$length, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999)), 1)
plot(density(df3.m$length))

# Look at capelin with very high or low condition
# I checked with Magen and all these seem good
df3.m[df3.m$rel.cond > 1.4, c("sex", "maturity", "rel.cond", "length", "weight")]
df3.m[df3.m$rel.cond < 0.7, c("sex", "maturity", "rel.cond", "length", "weight")]

plot(df3.m$rel.cond, df3.m$resids)
plot(boxplot(df3.m$rel.cond))

#females
df3.f <- filter(df3, sex == "Female")
m2.f <- lm(log10(weight) ~ log10(length), data= df3.f)

df3.f$fits <- fitted(m2.f)
glimpse(df3.f)

df3.f$rel.cond <- df3.f$weight/10^df3.f$fits
df3.f$resids <- df3.f$weight-10^df3.f$fits
plot(density(df3.f$rel.cond))

quantile(df3.f$weight, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999))
plot(density(df3.f$weight))
round(quantile(df3.f$length, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999)), 1)
plot(density(df3.f$length))

#df3.f[df3.f$rel.cond > 1.4, c("sex", "maturity", "rel.cond", "length", "weight")]
#df3.f[df3.f$rel.cond < 0.7, c("sex", "maturity", "rel.cond", "length", "weight")]

plot(df3.f$rel.cond, df3.f$resids)
plot(boxplot(df3.f$rel.cond))

df4 <- rbind(df3.m, df3.f)

# produce table of rel.cond +/- SD by year and nafo_div
df4 %>%
  group_by(year, nafo_div) %>%
  summarize(meanCond = round(mean(rel.cond),2), stdCond= round(sd(rel.cond),2)) %>% 
  unite(mean, meanCond:stdCond, sep = " +/- ") %>%
  spread(key = nafo_div, value = mean)


#nafo by sex - this graph indicates that M/F condition is very similar across years.
ggplot(data=df4) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  ylab("Relative condition") +
  xlab("Year") +
  facet_wrap(~sex, ncol=1) + 
  geom_hline(aes(yintercept = 1), colour = 'red') +
  theme_bw()
levels(as.factor(df4$year))

op <- par(mfrow = c(2,2))
plot(m2.f, add.smooth=F)
par(op)
par(mfrow = c(1,1))

# there appear to be a few outliers here but they should have minimal influence on the final outcome

# rel.cond by sex
ggplot(data=df4) +
  geom_boxplot(aes(x = sex, y = rel.cond, group = sex), notch=T)

# rel.cond by sex and year
ggplot(data=df4) +
  geom_boxplot(aes(x = sex, y = rel.cond, group = sex), notch=T) + 
  facet_wrap(~ year) + 
  geom_hline(aes(yintercept = 1), colour = 'red')

#month by year
ggplot(data=df4) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  facet_wrap(~month)

ggplot(data=df4) +
  geom_boxplot(aes(x = sex, y = rel.cond, group = sex)) + 
  facet_wrap(~maturity)

# produce output for Bayesain analysis 1999-2017
out <- df4 %>%
  group_by(year) %>%
  summarize(meanCond = round(mean(rel.cond),4), stdCond= round(sd(rel.cond),4), medCond = round(median(rel.cond), 4))

#write_csv(out, "E:\\Capelin\\2019 capelin assessment\\Condition\\condition_ag1_MF_out.csv")     
write_csv(out, "D:\\Keith\\capelin\\2019-code\\data\\condition_ag1_MF_out.csv")  


##Age2 M/F comparisons----
# first, consolidate the age 1 comparisions as above
glimpse(df3)
m2.a1 <- lm(log10(weight) ~ log10(length), data= df3)
m2.a1$fits <- fitted(m2.a1)


df.a1 <- df3
df.a1$rel.cond <- df.a1$weight/10^m2.a1$fits
df.a1$resids <- df.a1$weight-10^m2.a1$fits

# now do the same for age 2
#subset data- 
df.a2 <- df %>%
  filter(year > 1992 &  age == 2 & maturity != 6 & project != 10 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(23, 31, 32)) %>% #two-year old capelin after 1992, just project 23 sex == 1 &
  filter(!is.na(weight)) %>%
  filter(!is.na(length)) 

#make these variables factors
cols <- c("project", "nafo_div", "sex", "maturity")
df.a2 %<>%
  mutate_each_(funs(factor(.)),cols)
glimpse(df.a2)

#change values to something interpretable
df.a2$nafo_div <- as.factor(df.a2$nafo_div)
levels(df.a2$nafo_div)[levels(df.a2$nafo_div) == "23"] <- "2J"
levels(df.a2$nafo_div)[levels(df.a2$nafo_div) == "31"] <- "3K"
levels(df.a2$nafo_div)[levels(df.a2$nafo_div) == "32"] <- "3L"
levels(df.a2$sex)[levels(df.a2$sex) == "1"] <- "Male"
levels(df.a2$sex)[levels(df.a2$sex) == "2"] <- "Female"

# check levels of data to ensure that subset worked
levels(df.a2$project)
levels(as.factor(df.a2$sample_number))
levels(as.factor(df.a2$year))
levels(as.factor(df.a2$month))
levels(df.a2$nafo_div)
levels(df.a2$maturity)
range(df.a2$weight)
range(df.a2$length)
levels(as.factor(df.a2$age))

# Based on these plots and discussions with Megan Boucher, there is no need to subset.  All of these values are possible and reasonable.  Likely fish born late in the year and they never catch up.
plot(density(df.a2$weight))
plot(density(df.a2$length))


glimpse(df.a2)
m2.a2 <- lm(log10(weight) ~ log10(length), data= df.a2)
m2.a2$fits <- fitted(m2.a2)

df.a2$rel.cond <- df.a2$weight/10^m2.a2$fits
df.a2$resids <- df.a2$weight-10^m2.a2$fits

# Combine Age 1 and Age 2
df.a1_2 <- bind_rows(df.a1, df.a2)
#nafo by sex - this graph indicates that M/F condition is very similar across years.
ggplot(data=df.a1_2) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  ylab("Relative condition") +
  xlab("Year") +
  facet_wrap(~age, ncol=1) + 
  geom_hline(aes(yintercept = 1), colour = 'red') +
  theme_bw()


# rel.cond by sex
ggplot(data=df.a1_2) +
  geom_boxplot(aes(x = age, y = rel.cond, group = age), notch=T)

# rel.cond by sex and year
ggplot(data=df.a1_2) +
  geom_boxplot(aes(x = age, y = rel.cond, group = age), notch=T) + 
  facet_wrap(~ year) + 
  geom_hline(aes(yintercept = 1), colour = 'red')

#month by year
ggplot(data=df.a1_2) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  facet_grid(age~month)

ggplot(data=df.a1_2) +
  geom_boxplot(aes(x = age, y = rel.cond, group = age)) + 
  facet_wrap(~maturity)

out.a1_2 <- df.a1_2 %>%
  group_by(year, age) %>%
  summarize(meanCond = round(mean(rel.cond),4), stdCond= round(sd(rel.cond),4), medCond = round(median(rel.cond), 4))
temp <- spread(out.a1_2[c(1:3)], key = age, value = meanCond) %>% 
  setNames(c("year", "age1", "age2"))

# how related are age 1 and age 2 conditions?
plot(temp$age1 ~ temp$age2)
summary(lm(temp$age1 ~ temp$age2))
str(temp)

out.a1_2 <- out.a1_2 %>% 
  summarize(meanCond = mean(meanCond))

#write_csv(out.a1_2, "E:\\Capelin\\2019 capelin assessment\\Condition\\condition_ag1_2_MF_out.csv")  
write_csv(out.a1_2, "D:\\Keith\\capelin\\2019-code\\data\\condition_ag1_2_MF_out.csv")  

## prepare condition figure for age 1 male capelin----

#subset data
dfm1 <- df %>%
  filter(year > 1992 &  age == 1 & maturity != 6 & sex ==1 & project != 10 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(23, 31, 32)) %>% #one-year old capelin after 1992, just project 23 sex == 1 &
  filter(!is.na(weight)) %>%
  filter(!is.na(length)) 

#make these variables factors
cols <- c("project", "nafo_div", "sex", "maturity")
dfm1 %<>%
  mutate_each_(funs(factor(.)),cols)
glimpse(dfm1)

#change values to something interpretable
dfm1$nafo_div <- as.factor(dfm1$nafo_div)
levels(dfm1$nafo_div)[levels(dfm1$nafo_div) == "23"] <- "2J"
levels(dfm1$nafo_div)[levels(dfm1$nafo_div) == "31"] <- "3K"
levels(dfm1$nafo_div)[levels(dfm1$nafo_div) == "32"] <- "3L"
levels(dfm1$sex)[levels(dfm1$sex) == "1"] <- "Male"
levels(dfm1$sex)[levels(dfm1$sex) == "2"] <- "Female"

# check levels of data to ensure that subset worked
levels(df3$project)
levels(as.factor(df3$sample_number))
levels(as.factor(df3$year))
levels(as.factor(df3$month))
levels(df3$nafo_div)
levels(df3$maturity)
range(df3$weight)
range(df3$length)
levels(as.factor(df3$age))


# see email from H. Murphy and Methods of ResDoc - these are the values for filtering the dataset that Hannah and Maegen suggested based on biology
dfm1 <- dfm1 %>%
  group_by(year) %>%
  filter(weight > 2) %>% 
  filter(length > 80)

# run models for male and females seperately
#males
m1.m <- lm(log10(weight) ~ log10(length), data= dfm1)
dfm1$fits <- fitted(m1.m)

dfm1$rel.cond <- dfm1$weight/10^dfm1$fits
dfm1$resids <- dfm1$weight-10^dfm1$fits

male1.plot<-ggplot(data=dfm1) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  ylab("Relative condition") +
  xlab("Year") +
  facet_wrap(~nafo_div, ncol=1) + 
  geom_hline(aes(yintercept = 1), colour = 'red') +
  theme_bw() + 
  scale_y_continuous(limits=c(0.7,1.45),
                     expand=c(0,0),breaks=round(seq(0.7,1.4,0.2),1)) +
  scale_x_continuous(limits=c(1994.5,2018.5),expand=c(0,0)) + 
  theme(axis.text = element_text(face="bold",size=16),axis.title=element_text(face="bold",size=18),legend.title=element_text(face="bold",size=14),legend.text=element_text(size=12),axis.ticks=element_line(size=1),strip.text.y=element_text(size=12,face="bold"),strip.text.x=element_text(size=12,face="bold"))  


#png(file="E:\\Capelin\\2019 capelin assessment\\Condition\\male_age1_condition.png",width = 20,height=26,units="cm",res=300)
png(file="D:\\Keith\\capelin\\2019-code\\data\\male_age1_condition.png",width = 20,height=26,units="cm",res=300)
print(male1.plot)
dev.off()


# prepare condition figure for age 2 male capelin...

#subset data
dfm2 <- df %>%
  filter(year > 1992 &  age == 2 & maturity != 6 & sex ==1 & project != 10 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(23, 31, 32)) %>% #one-year old capelin after 1992, just project 23 sex == 1 &
  filter(!is.na(weight)) %>%
  filter(!is.na(length)) 

#make these variables factors
cols <- c("project", "nafo_div", "sex", "maturity")
dfm2 %<>%
  mutate_each_(funs(factor(.)),cols)
glimpse(dfm2)

#change values to something interpretable
dfm2$nafo_div <- as.factor(dfm2$nafo_div)
levels(dfm2$nafo_div)[levels(dfm2$nafo_div) == "23"] <- "2J"
levels(dfm2$nafo_div)[levels(dfm2$nafo_div) == "31"] <- "3K"
levels(dfm2$nafo_div)[levels(dfm2$nafo_div) == "32"] <- "3L"
levels(dfm2$sex)[levels(dfm2$sex) == "1"] <- "Male"
levels(dfm2$sex)[levels(dfm2$sex) == "2"] <- "Female"

# check levels of data to ensure that subset worked
levels(dfm2$project)
levels(as.factor(dfm2$sample_number))
levels(as.factor(dfm2$year))
levels(as.factor(dfm2$month))
levels(dfm2$nafo_div)
levels(dfm2$maturity)
range(dfm2$weight)
range(df3$length)
levels(as.factor(dfm2$age))


# see email from H. Murphy and Methods of ResDoc - these are the values for filtering the dataset that Hannah and Maegen suggested based on biology
dfm2 <- dfm2 %>%
  group_by(year) %>%
  filter(weight > 2) %>% 
  filter(length > 80)

# run models for male and females seperately
#males
m2.m <- lm(log10(weight) ~ log10(length), data= dfm2)
dfm2$fits <- fitted(m2.m)

dfm2$rel.cond <- dfm2$weight/10^dfm2$fits
dfm2$resids <- dfm2$weight-10^dfm2$fits

male2.plot<-ggplot(data=dfm2) +
  geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
  ylab("Relative condition") +
  xlab("Year") +
  facet_wrap(~nafo_div, ncol=1) + 
  geom_hline(aes(yintercept = 1), colour = 'red') +
  theme_bw() +
  scale_y_continuous(limits=c(0.7,1.45),
                     expand=c(0,0),breaks=round(seq(0.7,1.4,0.2),1)) +
  scale_x_continuous(limits=c(1994.5,2018.5),expand=c(0,0)) + 
  theme(axis.text = element_text(face="bold",size=16),axis.title=element_text(face="bold",size=18),legend.title=element_text(face="bold",size=14),legend.text=element_text(size=12),axis.ticks=element_line(size=1),strip.text.y=element_text(size=12,face="bold"),strip.text.x=element_text(size=12,face="bold")) 


#png(file="E:\\Capelin\\2019 capelin assessment\\Condition\\male_age2_condition.png",width = 20,height=26,units="cm",res=300)
png(file="D:\\Keith\\capelin\\2019-code\\data\\male_age2_condition.png",width = 20,height=26,units="cm",res=300)
print(male2.plot)
dev.off()
