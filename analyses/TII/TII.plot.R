#read in the TII values from mesquite
TII <- read.table("../results/TII Values.txt")

#plot the TII values against the number of taxa
#Inflection point 
plot(sort(TII[,2]), 
     pch=16, 
     cex=.5,
     xlab = "Taxa",
     ylab = "Taxanomic Instability Index",
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))
abline(h = 5000, lty = 3)

#export plot at 8" x 6"