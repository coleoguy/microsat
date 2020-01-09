library(beeswarm)
#name
mydata <-read.csv("../data/vert.data.csv", check.names = F)
dat <- mydata[,c(1:2,4:6,12:23)]
datcomp <- dat[complete.cases(dat),]
result <- c()
counter<-1
par(mfrow=c(2,3))
fit <- list()
for(i in 6:11){
  y <- datcomp[,i]
  A <- datcomp$`Type of Sequencing`
  fit <- aov(y ~ A)
  result[[counter]] <- anova(fit)
  beeswarm(y ~ A, las=1,
       xlab='Type of Sequencing',
       ylab="",
       main=paste(colnames(datcomp)[i], "count"),
       col = c("#fbb4ae", "#b3cde3", "#decbe4"), pch=16,
       method="swarm", cex.main=.9,
       cex=.8,spacing=.5,cex.lab=.8
       )
  top <- max(y)-0.1*max(y)
  text(y = top, 
       x = 1.1, cex=.8,
       paste("P-Value =",
             round(result[[counter]][1,5],digits = 3)))
  counter <- counter + 1
}
#export as pdf 7", 4"


