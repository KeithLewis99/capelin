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
library(RODBC)
library(forcats)

rm(list=ls())

# developing code to directly pull capelin data from the database...
# open connection to the capelin database
capelin_db<-odbcConnectAccess2007("J:\\Access Database Project\\Capelin Access database and program (limited access)\\capelin.mdb")
# get directory of tables
dir.tables<-sqlTables(capelin_db, tableType="TABLE")
# need to find out which tables are LSM
LSM.index<-grep("LSM.*", dir.tables$TABLE_NAME)
# need to find out which tables are Master
Master.index<-grep("Master.*", dir.tables$TABLE_NAME)
# need to find out which tables are Strat but not Frankie
# first find all the Frankie tables...
strat.all<-grep("Strat[0-9]{2}", dir.tables$TABLE_NAME)
strat.frankie<-grep("Frankie", dir.tables$TABLE_NAME)
Strat.index<-strat.all[which(!(strat.all %in% strat.frankie))]

# we have the indices, now pull in tables, combine them, then merge, then limit to Project Code 10/23
# pull in LSM data
LSM.table.list<-lapply(LSM.index, function(x) sqlFetch(capelin_db, dir.tables$TABLE_NAME[x])) %>% lapply(function(x) x%>% select("Sample Number", Year, "LSM Number", Length, Sex, Maturity))
LSM.table<-bind_rows(LSM.table.list)
# pull in Master data
Master.table.list<-lapply(Master.index, function(x) sqlFetch(capelin_db, dir.tables$TABLE_NAME[x])) %>% lapply(function(x) x%>% select("Sample Number", Year, Month, Day, "ICNAF Div", Gear, Country, Ship, "Trip Number", "Set Number", "Project Code"))
Master.table<-bind_rows(Master.table.list)
# pull in Strat data
Strat.table.list<-lapply(Strat.index, function(x) sqlFetch(capelin_db, dir.tables$TABLE_NAME[x])) %>% lapply(function(x) x%>% select("Sample Number", Year, "Otolith Number", "LSM Number", Length, Sex, Maturity, Weight, "Stomach Fullness", Age, Gonad))
Strat.table<-bind_rows(Strat.table.list)
close(capelin_db)

# clean up column names of each of the big tables...
colnames(Master.table)<-c("sample.number", "year", "month", "day", "nafo.div", "gear", "country", "ship", "trip.number", "set.number", "project.code")
colnames(LSM.table)<-c("sample.number", "year", "lsm.number", "length", "sex", "maturity")
colnames(Strat.table)<-c("sample.number", "year", "otolith.number", "lsm.number", "length", "sex", "maturity", "weight", "stomach.fullness", "age", "gonad")

lsm.master<-merge(Master.table, LSM.table, all.x=TRUE, all.y=TRUE)
lsm.master$year2<-lsm.master$year
lsm.master.strat<-merge(lsm.master, Strat.table, all.x=TRUE, all.y=TRUE) 
# I know at this point that there are bad matches in this dataset. I will clean this below when the dataset gets cleaned up for further analysis after the archive data gets pulled in. .

# developing code to directly pull capelin archive data from the database...
# open connection to the capelin database
capelin_arch_db<-odbcConnectAccess2007("J:\\Access Database Project\\Capelin archive\\capelin archive.mdb")
# get directory of tables to see what the table names look like for import below...
dir.tables.arch<-sqlTables(capelin_arch_db, tableType="TABLE")


# we have the indices, now pull in tables, combine them, then merge, then limit to Project Code 10/23
# pull in LSM data
LSM.table.arch<-sqlFetch(capelin_arch_db, "LSM_archive")
Master.table.arch<-sqlFetch(capelin_arch_db, "MASTER_archive") %>% select(Sample_Number, Year, Month, Day, NAFO_Div, Gear, Country, Ship, Trip_Number, Set_Number, Strat_Protocol, Project_Code)
Strat.table.arch<-sqlFetch(capelin_arch_db, "STRAT_archive")
close(capelin_arch_db)

# from data cleaning, I've seen that the outlier column causes mismatches. Exclude this column from the LSM.table.arch and Strat.table.arch data.frames

# addign in the second year column in order to help exclude mismatches
Master.table.arch$year2<-Master.table.arch$Year

lsm.master.arch<-merge(LSM.table.arch %>% select(-Outlier), Master.table.arch, all.x=TRUE, all.y=TRUE)
# during data cleaning, I determined that the 5 cases where the datasets mismatched there are 4 sample_numbers where capelin were not caught (or at least aren't in LSM) - 2005-> sample_number 26, 1979 -> sample_numbers 68, 98, and 2010 -> sample_number 632
# lsm.master.arch %>% filter(is.na(year2))
# lsm.master.arch %>% filter(is.na(Length))
# Additionally, we have a weird fish that comes from a sample_number that doesn't exist 2008, sample_number 2, LSM_number = 1
# clear out these records before doing next merge...
lsm.master.arch <- lsm.master.arch %>% filter(!is.na(Length)) %>% filter(!is.na(year2))
# merge in strat while taking out the Outlier column because it creates mismatches
lsm.master.strat.arch <- merge(lsm.master.arch, Strat.table.arch %>% select(-Outlier), all.x=TRUE, all.y=TRUE)

# later on I see that Year has to be greater than 1992, so can ignore 1992 and prior...
lsm.master.strat.arch <- lsm.master.strat.arch %>% filter(Year>1992)

# need to simplify both data sets so that they are more consistent with the existing data...
# at the same time, use year2 to filter out mismatches 
lsm.master.strat <- lsm.master.strat %>% filter(!is.na(year2)) %>% select(project.code, sample.number, year, month, nafo.div, length, sex, maturity, weight, stomach.fullness, age, gonad)
lsm.master.strat.arch <- lsm.master.strat.arch %>% filter(!is.na(year2)) %>% select(Project_Code, Sample_Number, Year, Month, NAFO_Div, Length, Sex, Maturity, Weight, Stomach_Fullness, Age, Gonad)

# need to get the names to match what df has for formating...
colnames(lsm.master.strat)<-c("project", "sample_number", "year", "month", "nafo_div", "length", "sex", "maturity", "weight", "stomach.fullness", "age", "gonad")
colnames(lsm.master.strat.arch)<-c("project", "sample_number", "year", "month", "nafo_div", "length", "sex", "maturity", "weight", "stomach.fullness", "age", "gonad")


# from here have to decide which dataset to use where they overlap and then progress onwards from here
head(lsm.master.strat)
head(lsm.master.strat.arch)
# table(lsm.master.strat$year) 2006 to 2019
# table(lsm.master.strat.arch$year) 1993 to 2015

# as more work has gone into cleaning up the archive data, will go with that data rather than current
lsm.master.strat<-lsm.master.strat %>% filter(year>2015)

# join the archived data and the current data
df_new<-rbind(lsm.master.strat.arch, lsm.master.strat)

# We are focused on a particualr groups of fish...
df3_new <- df_new %>% filter(year>1992 & age==1 & maturity!=6 & project==23 & month%in% c(10, 11, 12) & nafo_div%in%c(23, 31, 32)) %>% filter(!is.na(weight)) %>% filter(!is.na(length))


## load data----
# this from the ice-capelin-covariates file
df <- read_csv('data/capelin_condition_maturation_v1.csv') # read in full lsm data
# this step drops the duplicated project code variable
df <- df[c(1:5, 7:15)]
glimpse(df)
head(df)
df$weight
df_1617 <- read_csv('data/capelin_condition_maturation_v2.csv') # this is just pulling in more recent data and tehn reformatting so an rbind can be done. Probably can shrink this step...
glimpse(df_1617)
df_1617$stomach_fullness <- as.integer(df_1617$stomach_fullness)
df_1617$age <- as.integer(df_1617$age)
df_1617$gonad <- as.integer(df_1617$gonad)

head(df_1617)

df <- rbind(df, df_1617)
glimpse(df)
#View(tail(df, 100))

# save the combined data
write_csv(df, "data/condition_1979_2017.csv")

age <- 1

test<-c(0,1,0,1,0,1,0)
sum(test)

##Age1 M/F comparisons----
# there is considerable worry that using the males only will cause questions during the review process.  The intent of the following analysis is to test for differences between males and females.


#subset data
df3 <- df %>%
     filter(year > 1992 &  age == 1 & maturity != 6 & project != 10 & month %in% c(10, 11, 12) & nafo_div %in% c(23, 31, 32)) %>% filter(!is.na(weight)) %>% filter(!is.na(length)) 

#make these variables factors
cols <- c("project", "nafo_div", "sex", "maturity")
df3 %<>%
     mutate_each_(funs(factor(.)),cols)
glimpse(df3)

#change values to something interpretable
levels(df3$nafo_div)[levels(df3$nafo_div) == "23"] <- "2J"
levels(df3$nafo_div)[levels(df3$nafo_div) == "31"] <- "3K"
levels(df3$nafo_div)[levels(df3$nafo_div) == "32"] <- "3L"
levels(df3$sex)[levels(df3$sex) == "1"] <- "Male"
levels(df3$sex)[levels(df3$sex) == "2"] <- "Female"

#change values to something interpretable
df3_new2<-df3_new %>% mutate_at(vars(cols), funs(as.factor)) %>% mutate(nafo_div=fct_recode(nafo_div, "2J" = "23", "3K" = "31", "3L" = "32")) %>% mutate(sex=fct_recode(sex, "Male"="1", "Female"="2", "Unknown"="3"))


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

# check levels of data to ensure that subset worked - repeat to see that it worked for the direct access database
levels(df3_new2$project)
levels(as.factor(df3_new2$sample_number))
levels(as.factor(df3_new2$year))
levels(as.factor(df3_new2$month))
levels(df3_new2$nafo_div)
levels(df3_new2$maturity)
range(df3_new2$weight)
range(df3_new2$length)
levels(as.factor(df3_new2$age))


# see email from H. Murphy and Methods of ResDoc - these are the values for filtering the dataset that Hannah and Maegen suggested based on biology
df3 <- df3 %>%
     group_by(year) %>%
     filter(weight > 2) %>% 
     filter(length > 80)

df3_new2 <- df3_new2 %>% 
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




# run models for male and females seperately - direct from database
#males
df3.new.m <- filter(df3_new2, sex == "Male")
m2.new.m <- lm(log10(weight) ~ log10(length), data= df3.new.m)
df3.new.m$fits <- fitted(m2.new.m)

m2.new.m$fitted.values
str(m2.new.m)
df3.new.m$rel.cond <- df3.new.m$weight/10^df3.new.m$fits
df3.new.m$resids <- df3.new.m$weight-10^df3.new.m$fits


ggplot(data=df3.new.m) + geom_point(aes(x=log10(length), y = log10(weight), colour=nafo_div)) + 
  facet_wrap(~nafo_div)


# test for outliers of condition
quantile(df3.m$weight, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999))
plot(density(df3.m$weight))
round(quantile(df3.m$length, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999)), 1)
plot(density(df3.m$length))


# test for outliers of condition - direct from database
quantile(df3.new.m$weight, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999))
lines(density(df3.new.m$weight), col="red")
round(quantile(df3.new.m$length, c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999)), 1)
lines(density(df3.new.m$length), col="red")


# Look at capelin with very high or low condition
# I checked with Magen and all these seem good
df3.m[df3.m$rel.cond > 1.4, c("sex", "maturity", "rel.cond", "length", "weight")]
df3.m[df3.m$rel.cond < 0.7, c("sex", "maturity", "rel.cond", "length", "weight")]

plot(df3.m$rel.cond, df3.m$resids)
plot(boxplot(df3.m$rel.cond))



# Look at capelin with very high or low condition - direct access
# I checked with Magen and all these seem good
df3.new.m[df3.new.m$rel.cond > 1.4, c("sex", "maturity", "rel.cond", "length", "weight")]
df3.new.m[df3.new.m$rel.cond < 0.7, c("sex", "maturity", "rel.cond", "length", "weight")]

plot(df3.new.m$rel.cond, df3.new.m$resids)
plot(boxplot(df3.new.m$rel.cond))


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
View(out)

write_csv(out, "data/condition_ag1_MF_out.csv")     



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

write_csv(out.a1_2, "data/condition_ag1_2_MF_out.csv")  


######## Additional analyses----
# confirm Fran's trend that capelin size hasn't really decreased, its been the age that has.
# 
# glimpse(df)
# levels(as.factor(df$age))
# 
# df_age <- df %>%
#      filter(maturity != 6 & as.factor(month) %in% c("10", "11", "12") & project != 10 & as.factor(nafo_div) %in% c(31, 32) & age != 0 & age !=7) %>% 
#      filter(!is.na(age)) 
# levels(as.factor(df_age$age))
# 
# 
# df_age_avg <- df_age %>%
#      group_by(year) %>% 
#      summarize(mean_age = mean(age), sd_age = sd(age))
# 
# # this does not match the 2008 doc because these are sampled capelin (I think) and not commercial capelin
# p <- ggplot(data=df_age_avg, aes(x=year, y = mean_age))
# p <- p + geom_point()
# p <- p + geom_errorbar(data=df_age_avg, aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age))
# p <- p + geom_line(aes(x=year, y = mean_age))
# p <- p + theme_tufte() + xlab("Year") + ylab("Age (years)")
# p
# 
# # age from the acoustic survey
# df_prop <- df_age %>%
#      group_by(year, age) %>% 
#      summarise(n = n()) %>% 
#      mutate(freq = n/sum(n))
# 
# p <- ggplot()
# p <- p + geom_col(data = df_prop, aes(x = year, y = freq, fill = as.factor(age)))
# p <- p + scale_fill_discrete(name = "Age")
# p <- p + theme(axis.title.x = element_text(face = "bold", size = 25), 
#                axis.text.x = element_text(size = 20), 
#                axis.title.y = element_text(face = "bold", size = 25),
#                axis.text.y = element_text(size = 20),
#                legend.title = element_text(size = 25, face = "bold"),
#                legend.text = element_text(size = 20)) + guides(colour = guide_legend(override.aes = list(size=25))) + xlab("Year") + ylab("Proportion")
# #p <- p + theme_bw(base_size = 20)
# p
# ggsave("age_proportion.jpeg", width=8, height=5, units="in")
# 
# # for length
# df_length <- df %>%
#      filter(maturity != 6 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(31, 32) & age >= 1 & age < 6) %>% 
#      filter(!is.na(length)) %>% 
#      filter(!is.na(age))
# 
# df_length_avg <- df_length %>%
#      group_by(year) %>% 
#      summarize(mean_length = mean(length), sd_length = sd(length))
# 
# p <- ggplot(data=df_length_avg, aes(x=year, y = mean_length))
# p <- p + geom_point()
# p <- p + geom_errorbar(aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length))
# p <- p + geom_line(aes(x=year, y = mean_length))
# p
# 
# # length by age
# df_length_age_avg <- df_length %>%
#      group_by(year, age) %>% 
#      filter(year > 1995) %>%
#      summarize(mean_length = mean(length), sd_length = sd(length)) %>% 
#      mutate(cv = sd_length/mean_length)
# 
# p <- ggplot(data=df_length_age_avg, aes(x=year, y = mean_length, colour = as.factor(age)))
# p <- p + geom_point()
# p <- p + geom_line()
# p <- p + scale_colour_discrete(name = "Age", breaks = c(5, 4, 3, 2, 1))
# p <- p + theme_bw(base_size = 25) + xlab("Year") + ylab("Length (mm)")
# p
# ggsave("length_age.jpeg", width=10, height=8, units="in")
# 
# p <- ggplot(data=df_length_age_avg, aes(x=year, y = cv, colour = as.factor(age)))
# p <- p + geom_point()
# p <- p + geom_line()
# p <- p + scale_colour_discrete(name = "Age", breaks = c(5, 4, 3, 2, 1))
# p <- p + theme_bw() + xlab("Year") + ylab("Coefficient of variation")
# p
# 
# 
# # for weight
# df_weight <- df %>%
#      filter(year > 1991 & maturity != 6 & as.factor(month) %in% c("10", "11", "12") & as.factor(nafo_div) %in% c(31, 32) & age >= 1 & age < 6) %>% 
#      filter(!is.na(weight)) %>% 
#      filter(!is.na(age))
# 
# df_weight_avg <- df_weight %>%
#      group_by(year) %>% 
#      summarize(mean_weight = mean(weight), sd_weight = sd(weight))
# 
# p <- ggplot(data=df_weight_avg, aes(x=year, y = mean_weight))
# p <- p + geom_point()
# p <- p + geom_errorbar(aes(ymin=mean_weight-sd_weight, ymax=mean_weight+sd_weight))
# p <- p + geom_line(aes(x=year, y = mean_weight))
# p
# 
# # weight by age
# df_weight_age_avg <- df_weight %>%
#      group_by(year, age) %>% 
#      summarize(mean_weight = mean(weight), sd_weight = sd(weight)) %>% 
#      mutate(cv = sd_weight/mean_weight)
# 
# 
# p <- ggplot(data=df_weight_age_avg, aes(x=year, y = mean_weight, colour = as.factor(age)))
# p <- p + geom_point()
# p <- p + geom_line()
# p <- p + scale_colour_discrete(name = "Age", breaks = c(5, 4, 3, 2, 1))
# p <- p + theme_bw() + xlab("Year") + ylab("Weight (mm)")
# p
# 
# p <- ggplot(data=df_weight_age_avg, aes(x=year, y = cv, colour = as.factor(age)))
# p <- p + geom_point()
# p <- p + geom_line()
# p <- p + scale_colour_discrete(name = "Age", breaks = c(5, 4, 3, 2, 1))
# p <- p + theme_bw() + xlab("Year") + ylab("Coefficient of variation")
# p
# 
