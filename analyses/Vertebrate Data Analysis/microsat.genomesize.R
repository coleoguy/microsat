#read in vertebrate data
mydata <-read.csv("../data/vert.data.csv", check.names = F)

#only keep the columns needed
dat <- mydata[,c(1:2,4:6,8,12:23)]

#only keep the data that has all information
datcomp <- dat[complete.cases(dat),]

#plot the total microsatellite content to the genome size
plot(datcomp$`Genome Size`, datcomp$`2-6mer`,
     xlab = "Genome Size in Gbp",
     ylab = "Total 2-6mer Content in bp", 
     col = "maroon")
