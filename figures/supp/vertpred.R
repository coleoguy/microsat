#read in vertebrate data
dat <-read.csv("../../data/verts/vert.csv")

#only keep the columns needed
dat <- dat[,c(1:2,4:6,8,12:23)]

#only keep the data that has all information
dat <- dat[complete.cases(dat),]

colnames(dat)[18]<-"2.6mermb"

fit1 <- summary(lm(dat$`2.6mermb` ~ dat$Genome.Size))
fit2 <- summary(lm(dat$`2.6mermb` ~ dat$Contig.N50))
fit3 <- summary(lm(dat$`2.6mermb` ~ dat$Scaffold.N50))


par(mfcol=c(3,1))
#par(oma=c(2,1,1,1))
par(mar=c(4,4,1,1))
alpha <- 100
ylab="2-6mer (bp/mb)"
#plot the total microsatellite content to the genome size
plot(dat$`2.6mermb` ~ dat$Genome.Size,
     xlab = "Genome Size in Gbp",
     ylab = ylab, 
     col = rgb(250, 159, 181,alpha, maxColorValue = 255), pch=16)
pval <- round(fit1$coefficients[2,4], digits=3)
text(x=max(dat$Genome.Size), y=max(dat$`2.6mermb`)-50,
     labels= paste("p-value=", pval, sep=""), pos=2,cex=.7)

#plot the total microsatellite content to the contig N50
plot(dat$`2.6mermb` ~ dat$Contig.N50,
     xlab = "Contig N50",
     ylab = ylab, 
     col = rgb(250, 159, 181,alpha, maxColorValue = 255), pch=16)
pval <- round(fit2$coefficients[2,4], digits=3)
text(x=max(dat$Contig.N50), y=max(dat$`2.6mermb`)-50,
     labels= paste("p-value=", pval, sep=""), pos=2,cex=.7)

#plot the total microsatellite content to the scaffold N50
plot(dat$`2.6mermb` ~ dat$Scaffold.N50,
     xlab = "Scaffold N50",
     ylab = ylab, 
     col = rgb(250, 159, 181,alpha, maxColorValue = 255), pch=16)
pval <- round(fit3$coefficients[2,4], digits=3)
text(x=max(dat$Scaffold.N50), y=max(dat$`2.6mermb`)-50,
     labels= paste("p-value=", pval, sep=""), pos=2,cex=.7)

# export 3"x8"