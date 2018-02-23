## somatic lipid content data taken from Montevecchi & Piatt 1984


fat <- read.table('presentation/capfat.txt', header = T, sep = '\t')
for (i in 1:12) {
  fat$month[i] <- ifelse(fat$month[i] < 8, fat$month[i] + 
                           12, fat$month[i])
}
Imon <- order(fat$month)


par(mfrow = c(1, 1))
par(bty = "o")
par(mar = c(4, 4, 0, 1))  # set margins(bottom,left,top,right)
plot(fat$month[Imon], fat$lipid[Imon], type = "l", 
     lwd = 1.5, col = "red", ylim = c(0, 22.5), yaxp = c(0, 
                                                         18, 2), ylab = "Somatic lipid content (%)", 
     xlab = "Month", xaxt = "n", yaxt = "n")
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
text(18.5, 20.5, "spawning", col = "orange", adj = 0.5, 
     cex = 1)
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



## next plot

year <- c(1, 2, 3, 4)
month <- c(1:12)
month_seq <- c(1:48)

year <- sort(rep(year, 12))
month <- rep(month, 4)

x <- as.data.frame(cbind(year, month, month_seq))
x <- x[6:44, ]
lab_1 <- rep(c("Jun", "Sep", "Dec", "Mar"), 4)
lab_2 <- lab_1[1:13]
plot(x$month_seq, x$year, type='n', xaxt = 'n', yaxt = 'n', xlab = 'Month', ylab = '')
axis(1, at = seq(6, 44, 3), labels = lab_2, col = "#0000ff00", 
     col.ticks = "black")

arrows(6, 3, 18, 3, length = 0.1, angle = 20, 
       code = 3, col = "blue", lwd = 1.5)
text(8, 3.2, "year(t)", col = "#FF6699", adj = 0.5, 
     cex = 1)
text(8, 2.5, "larval \n emergence", col = "#FF6699", adj = 0.5,      cex = 1)
text(12, 2.3, bquote(atop(italic("Psuedocalanus"), "emergece")), col = "#FF6699", adj = 0.5, cex = 1)
arrows(18, 3, 30, 3, length = 0.1, angle = 20, 
       code = 3, col = "red", lwd = 1.5)
text(20, 3.2, "year(t+1)", col = "#FF6699", adj = 0.5, 
     cex = 1)
text(28, 2.5, "condition", col = "#FF6699", adj = 0.5, 
     cex = 1)

arrows(30, 3, 43, 3, length = 0.1, angle = 20, 
       code = 3, col = "black", lwd = 1.5)
text(32, 3.2, "year(t+2)", col = "#FF6699", adj = 0.5, 
     cex = 1)
text(34, 2.5, "tice", col = "#FF6699", adj = 0.5, 
     cex = 1)
text(36, 2.3, "spring survey", col = "#FF6699", adj = 0.5, 
     cex = 1)

