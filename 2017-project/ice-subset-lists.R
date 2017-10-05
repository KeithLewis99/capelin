################################################################
#  Script written by Keith Lewis (Keith.Lewis@dfo-mpo.gc.ca)  #
#  Created 2017-08, R version 3.X.x (201X-XX-XX)             #
#  Last modified by Keith Lewis (2017-10-05) #
################################################################

## subset the ice data
# create "model sets"

ct <- c(1:10, 9.5)
sa <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "1.", "4.", "7.", "8.", "9.") # don't have "ice of land origin" or "undetermined or unknown"
sb <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "1.", "4.", "7.", "8.", "9.") # don't have "ice of land origin" or "undetermined or unknown"
m1 <- list(ct=ct, sa=sa, sb=sb)
m1

ct <- c(1:10, 9.5)
sa <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
sb <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
m2 <- list(ct=ct, sa=sa, sb=sb)
m2

ct <- c(3:10, 9.5)
sa <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
sb <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
m3 <- list(ct=ct, sa=sa, sb=sb)
m3

ct <- c(7:10, 9.5)
sa <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
sb <- c("5", "6", "7", "8", "9", "1.", "4.", "7.")
m4 <- list(ct=ct, sa=sa, sb=sb)
m4

ct <- c(3:10, 9.5)
sa <- c("7", "8", "9", "1.", "4.", "7.")
sb <- c("7", "8", "9", "1.", "4.", "7.")
m5 <- list(ct=ct, sa=sa, sb=sb)
m5

ct <- c(7:10, 9.5)
sa <- c("7", "8", "9", "1.", "4.", "7.")
sb <- c("7", "8", "9", "1.", "4.", "7.")
m6 <- list(ct=ct, sa=sa, sb=sb)
m6

# this was meant to make a mutually exclusive data set to compare m2 to m1.  However, this is not the case: see this example and work it out for m1, m2, m4
# CT SA
# 4   5
# 2   5
save(m1, m2, m3, m4, m5, m6, file = "output-processing/subset-lists.Rdata")

rm(ct)
rm(sa)
rm(sb)