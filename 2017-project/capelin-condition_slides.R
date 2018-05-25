rm(list=ls())

df <- read_csv('data/condition_1979_2017.csv')
glimpse(df)
head(df)

cols <- c("project", "nafo_div", "sex", "maturity")
df1 <- df %<>%
     mutate_each_(funs(factor(.)),cols)

df1 <- subset(df1, !is.na(weight))
df1 <- subset(df1, !is.na(length))
glimpse(df1)

levels(df1$nafo_div)[levels(df1$nafo_div) == "23"] <- "2J"
levels(df1$nafo_div)[levels(df1$nafo_div) == "31"] <- "3K"
levels(df1$nafo_div)[levels(df1$nafo_div) == "32"] <- "3L"
levels(df1$sex)[levels(df1$sex) == "1"] <- "Male"
levels(df1$sex)[levels(df1$sex) == "2"] <- "Female"

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
     summarize(uci = quantile(weight, c(0.99))) 
View(quant)
df1 <- df1 %>%
     group_by(year) %>%
     right_join(x=df1, y=quant, by ="year") %>%
     filter(weight < uci)



########################STOP - Hammer time###############
#Make changes
df2 <- df1 %>%
     filter(year > 1992 & sex == "Male" & age == 2 &
                 maturity != 6 & 
                 project != 10 & 
                 as.factor(month) %in% c("10", "11", "12") & 
                 nafo_div %in% c("2J", "3K", "3L"))
df2 <- filter(df2, rel.cond < 1.4 & rel.cond > 0.7)

levels(as.factor(df1$month))

## LINEAR MODELLING----
m1 <- lm(log10(weight) ~ log10(length), data= df2)
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


##output----
df2$fits <- fitted(m1)
m1$fitted.values
str(m1)
df2$rel.cond <- df2$weight/10^df2$fits

df2$resids <- df2$weight-10^df2$fits

plot(df2$rel.cond, df2$resids)

ggplot(data=df2) +
     geom_boxplot(aes(x = year, y = rel.cond, group = year)) + 
     ylab("Relative condition") +
     xlab("Year") +
     facet_wrap(~nafo_div, ncol=1) + 
     theme_bw()+
     geom_hline(aes(yintercept = 1), colour = 'red')


ggsave(filename = "condition/age2-fall.png", width=10, height=8, units="in")
