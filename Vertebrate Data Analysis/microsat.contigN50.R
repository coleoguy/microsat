#read in vertebrate data
mydata <-read.csv("../data/vert.data.csv", check.names = F)

#only keep the colmns needed
dat <- mydata[,c(1:2,4:6,8,12:23)]

#only keep data that has all information
datcomp <- dat[complete.cases(dat),]

#plot the total microsatellite content to the contig N50
plot(datcomp$`Contig N50`, datcomp$`2-6mer`,
     xlab = "Contig N50",
     ylab = "Total 2-6mer Content in bp", 
     col = "maroon")
