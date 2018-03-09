## somatic lipid content data taken from Montevecchi & Piatt 1984


fat <- read.table('presentation/capfat.txt', header = T, sep = '\t')
for (i in 1:12) {
  fat$month[i] <- ifelse(fat$month[i] < 8, fat$month[i] + 
                           12, fat$month[i])
}
Imon <- order(fat$month)

png("presentation/Buren_fig6.png", family = "Arial Black")
par(mfrow = c(1, 1))
par(bty = "o")
par(mar = c(4, 4.5, 0, 1))  # set margins(bottom,left,top,right)
plot(fat$month[Imon], fat$lipid[Imon], type = "l", 
     lwd = 1.5, col = "red", ylim = c(0, 22.5), yaxp = c(0, 
                                                         18, 2), ylab = "Somatic lipid content (%)", 
     xlab = "Month", xaxt = "n", yaxt = "n", cex.lab = 1.5)
axis(1, at = c(8, 10, 12, 14, 16, 18), labels = c("Aug", 
                                                  "Oct", "Dec", "Feb", "Apr", "Jun"), col = "#0000ff00", 
     col.ticks = "black")

axis(2, at = c(0, 9, 18), labels = T, , col = "#0000ff00", 
     col.ticks = "black")
axis(2, at = c(4.5, 13.5), labels = F, , col = "#0000ff00", 
     col.ticks = "black", tck = -0.01)
axis(1, at = seq(9, 19, by = 2), labels = F, tck = -0.01, 
     col = "#0000ff00", col.ticks = "black")


arrows(13, 18.5, 15, 18.5, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 1.5)
text(14, 17.5, "not feeding", col = "blue", adj = 0.5, 
     cex = 1)
arrows(15, 18.5, 18, 18.5, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 1.5)
text(16.5, 17.5, "feeding", col = "blue", adj = 0.5, 
     cex = 1)
text(16.5, 16.5, "gonad development", col = "blue", 
     adj = 0.5, cex = 1)
arrows(18, 18.5, 19, 18.5, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 1.5)
text(10.5, 17.5, "feeding", col = "blue", adj = 0.5, 
     cex = 1)
arrows(8, 18.5, 13, 18.5, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 1.5)
text(18.5, 17.5, "not feeding", col = "blue", adj = 0.5, 
     cex = 1)
arrows(18, 20, 19, 20, length = 0.05, angle = 90, code = 3, 
       col = "orange", lwd = 1.5)
text(18.75, 20.5, "spawning", col = "orange", adj = 0.5, 
     cex = 0.8)
arrows(16.75, 20, 17.8, 20, length = 0.05, angle = 90, 
       code = 3, col = "black", lwd = 1.5)
# text(5.45,21.2,'acoustic',col='black',adj=0.5,cex=1)
text(17.275, 20.5, "survey", col = "black", adj = 0.5, 
     cex = 1)
arrows(13.7, 20, 15.2, 20, length = 0.05, angle = 90, 
       code = 3, col = "green4", lwd = 1.5)
text(14.45, 20.5, "spring bloom", col = "green4", adj = 0.5, 
     cex = 1)
emerge <- expression(paste(italic(Calanus), " emergence"))
arrows(13, 21.5, 15.5, 21.5, length = 0.05, angle = 90, 
       code = 3, col = "#FF6699", lwd = 1.5)
text(14.3, 22, emerge, col = "#FF6699", adj = 0.5, 
     cex = 1)
text(9, 2, "Buren et al. 2014 \n PlosOne")
dev.off()


## next plot
install.packages("extrafont")
library(extrafont)
font_import()
y
year <- c(1, 2, 3, 4)
month <- c(1:12)
month_seq <- c(1:48)

year <- sort(rep(year, 12))
month <- rep(month, 4)

x <- as.data.frame(cbind(year, month, month_seq))
x <- x[1:36, ]
lab_1 <- rep(c("Dec", "Mar", "Jun", "Sep"), 4)
lab_2 <- lab_1[1:13]

png("presentation/variables.png", family = "Arial Black")
par(mar = c(4, 4, 0, 1))  # set margins(bottom,left,top,right)

plot(x$month_seq, x$year, type='n', xaxt = 'n', yaxt = 'n', xlab = 'Month', ylab = '')
axis(1, at = seq(1, 37, 3), labels = lab_2, col = "#0000ff00", col.ticks = "black")

#year(t+2)
arrows(1, 2.5, 12, 2.5, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 3)
text(3, 2.7, "year(t-2)", col = "blue", adj = 0.5, 
     cex = 1.5)
text(9, 2.3, "larval \n emergence", col = "blue", adj = 0.5,      cex = 1.5)
text(10, 2, bquote(atop(italic("Psuedocalanus"), "emergence")), col = "blue", adj = 0.5, cex = 1.5)

#year(t-1)
arrows(12, 2.5, 24, 2.5, length = 0.1, angle = 20, 
       code = 3, col = "red", lwd = 3)
text(14, 2.7, "year(t-1)", col = "red", adj = 0.5, 
     cex = 1.5)
text(23, 2.3, "condition", col = "red", adj = 0.5, 
     cex = 1.5)

#year(t)
arrows(24, 2.5, 36, 2.5, length = 0.1, angle = 20, 
       code = 3, col = "black", lwd = 3)
text(26, 2.7, "year(t)", col = "black", adj = 0.5, 
     cex = 1.5)
text(28, 2.3, expression("t"[italic(ice)]), col = "black", adj = 0.5, 
     cex = 1.5)
text(31, 2.1, "spring \n survey", col = "black", adj = 0.5, 
     cex = 1.5)
dev.off()


#plot(1:10, xlab=expression('hi'[5]*'there'[6]^8*'you'[2]))

plot(density(rnorm(10000, 0, 10000)))
plot(density(rnorm(10000, 0, 1000)))
plot(density(rnorm(10000, 0, 10)))
plot(density(rnorm(10000, 0, 1)))
lines(density(rnorm(10000, 0.35, 0.138)))

plot(density(rgamma(10000, 5, 1/3)))
lines(density(rgamma(10000, 1, 1/10)))
lines(density(rgamma(10000, 10, 1/10)))
curve(dgamma(x, .001, .001))


