#read in vertebrate data
dat <-read.csv("../data/verts/vert.csv")

#only keep the columns needed
dat <- dat[,c(1:2,4:6,8,12:23)]

#only keep the data that has all information
dat <- dat[complete.cases(dat),]

colnames(dat)[18]<-"2.6mermb"

summary(lm(dat$Genome.Size ~ dat$`2.6mermb`))
summary(lm(dat$Contig.N50  ~ dat$`2.6mermb`))
summary(lm(dat$Scaffold.N50 ~ dat$`2.6mermb`))


(par(mfcol=c(3,1)))
#plot the total microsatellite content to the genome size
plot(dat$Genome.Size, dat$`2.6mermb`,
     xlab = "Genome Size in Gbp",
     ylab = "Total 2-6mer Content in bp", 
     col = "maroon")
#plot the total microsatellite content to the contig N50
plot(datcomp$`Contig N50`, datcomp$`2-6mer`,
     xlab = "Contig N50",
     ylab = "Total 2-6mer Content in bp", 
     col = "maroon")

#plot the total microsatellite content to the scaffold N50
plot(datcomp$`Scaffold N50`, datcomp$`2-6mer`,
     xlab = "Scaffold N50",
     ylab = "Total 2-6mer Content in bp", 
     col = "maroon")
