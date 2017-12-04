x <- c(-100:100)
alpha <- 1 # increaes alpha, decrease |y|, negative turns parabola upside down
beta <- 50 # increase beta, decrease  |y|
gamma <- 2 # decrease gamma, decrease y
delta <- 1
z <- rnorm(201, mean=50, sd=10)
w <- rnorm(201, mean=30, sd=20)
hist(z)


y1 <- alpha*x
y2 <- 1-x/beta
y3 <- alpha*x*(1-x/beta)
y4 <- alpha*x*(1-x/beta)
y5 <- y4 + gamma*z
y6 <- y5 + delta*w

df <- cbind(y5, y4, x, z, w)
head(df)

plot(x, y1)
plot(x, y2)
plot(x, y3)
plot(x, y4)
plot(x, y5)
lines(x, y5)
plot(x, y6)
lines(x, y6)


alp <- MaxTice4$optim_ls$`MaxTice-m1`$cdf$par[1]
bet <- MaxTice4$optim_ls$`MaxTice-m1`$cdf$par[2]

x <- seq(0, 10, 0.1)
y <- alp*x*(1-x/bet)
plot(x, y)


lx <- 1:1000

yln <- log(lx)
ylg <- log10(lx)
ylg2 <- log2(lx)

plot(lx, yln)
