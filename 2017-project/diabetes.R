##https://www.kaggle.com/avehtari/bayesian-logistic-regression-with-rstanarm/notebook

# file preview shows a header row
diabetes <- read.csv("diabetes.csv", header = TRUE)

# first look at the data set using summary() and str() to understand what type of data are you working
# with
summary(diabetes)
str(diabetes)

diabetes$Outcome <- factor(diabetes$Outcome)

# removing those observation rows with 0 in any of the variables
for (i in 2:6) {
     diabetes <- diabetes[-which(diabetes[, i] == 0), ]
}
# scale the covariates for easier comparison of coefficient posteriors
for (i in 1:8) {
     diabetes[i] <- scale(diabetes[i])
}

# modify the data column names slightly for easier typing
names(diabetes)[7] <- "dpf"
names(diabetes) <- tolower(names(diabetes))

n=dim(diabetes)[1]
p=dim(diabetes)[2]
str(diabetes)
print(paste0("number of observations = ", n))
print(paste0("number of predictors = ", p))

# preparing the inputs
x <- model.matrix(outcome ~ . - 1, data = diabetes)
y <- diabetes$outcome


library(rstanarm)
t_prior <- student_t(df = 7, location = 0, scale = 2.5)
post1 <- stan_glm(outcome ~ ., data = diabetes,
                  family = binomial(link = "logit"), 
                  prior = t_prior, prior_intercept = t_prior,
                  seed = 1)

library(ggplot2)
pplot<-plot(post1, "areas", prob = 0.95, prob_outer = 1)
pplot+ geom_vline(xintercept = 0)

round(coef(post1), 2)
round(posterior_interval(post1, prob = 0.9), 2)


library(loo)
(loo1 <- loo(post1))

# Predicted probabilities
linpred <- posterior_linpred(post1)
preds <- posterior_linpred(post1, transform=TRUE)
pred <- colMeans(preds)
pr <- as.integer(pred >= 0.5)



## from Gelmen - http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
     y <- get_y(fit)
     ypred <- posterior_linpred(fit, transform = TRUE)
     if (family(fit)$family == "binomial" && NCOL(y) == 2) {
          trials <- rowSums(y)
          y <- y[, 1]
          ypred <- ypred %*% diag(trials)
     }
     e <- -1 * sweep(ypred, 2, y)
     var_ypred <- apply(ypred, 1, var)
     var_e <- apply(e, 1, var)
     var_ypred / (var_ypred + var_e)
}

# works but only when outcome is not a factor
print(median(bayes_R2(post1)))

fit <- post1

## Example
M1 <- stan_glm(y ~ x)
print(median(bayes_R2(M1)))

library(caret)
# confusion matrix
confusionMatrix(pr, y)[2:3]
# posterior classification accuracy
round(mean(xor(pr,as.integer(y))),3)
# posterior balanced classification accuracy
round((mean(xor(pr[y==0]>0.5,as.integer(y[y==0])))+mean(xor(pr[y==1]>0.5,as.integer(y[y==1]))))/2,3)


# PSIS-LOO weights
log_lik=log_lik(post1, parameter_name = "log_lik")
psis=psislw(-log_lik)
#plot(psis$pareto_k)
#plot(psis$lw_smooth[,1],linpred[,1])
# LOO predictive probabilities
ploo=colSums(preds*exp(psis$lw_smooth))
# LOO classification accuracy
round(mean(xor(ploo>0.5,as.integer(y))),3)
# LOO balanced classification accuracy
round((mean(xor(ploo[y==0]>0.5,as.integer(y[y==0])))+mean(xor(ploo[y==1]>0.5,as.integer(y[y==1]))))/2,2)

plot(pred,ploo)
