library(geiger)
library(phytools)

#read in tree
trees <- read.nexus("../data/trees/post.nex")
trees <- trees[[sample(1:100, 1)]]
trees <- drop.tip(phy=trees, tip="B.terrestris") 
trees <- drop.tip(phy=trees, tip="Plutella_xylostella")

# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/ssr.inference/micRocounter_results_TII_typecentromere.csv", 
                          row.names = 4)

#set what models we want to run in fitContinuous
models <-c("BM","OU","EB","lambda","kappa","delta","white")

#call out what data we need for each trait

#bp.Mbp vector
bp.Mbp <- microsat.cent$bp.Mbp
# make named vector for bpMbp coontent
names(bp.Mbp) <- row.names(microsat.cent)

#twomer vector
twomers <- microsat.cent$twomers
# make named vector for bpMbp coontent
names(twomers) <- row.names(microsat.cent)

#threemer vector
threemers <- microsat.cent$threemers
# make named vector for bpMbp coontent
names(threemers) <- row.names(microsat.cent)

#fourmer vector
fourmers <- microsat.cent$fourmers
# make named vector for bpMbp coontent
names(fourmers) <- row.names(microsat.cent)

#fivemer vector
fivemers <- microsat.cent$fivemers
# make named vector for bpMbp coontent
names(fivemers) <- row.names(microsat.cent)

#sixmer vector
sixmers <- microsat.cent$sixmers
# make named vector for bpMbp coontent
names(sixmers) <- row.names(microsat.cent)

#totalmer vector
totalmers <- microsat.cent$all
# make named vector for bpMbp coontent
names(totalmers) <- row.names(microsat.cent)

#how to ouput all of the data for each model for each trait
#fitmodels for bp.Mbp
results.bp.Mbp <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.bp.Mbp[[i]] <- fitContinuous(phy=trees, 
                                       dat=bp.Mbp, 
                                       model=models[i])
}

#fitmodels for twomers
results.twomers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.twomers[[i]] <- fitContinuous(phy=trees, 
                                       dat=twomers, 
                                       model=models[i])
}

#fitmodels for threemers
results.threemers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.threemers[[i]] <- fitContinuous(phy=trees, 
                                       dat=threemers, 
                                       model=models[i])
}

#fitmodels for fourmers
results.fourmers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.fourmers[[i]] <- fitContinuous(phy=trees, 
                                       dat=fourmers, 
                                       model=models[i])
}

#fitmodels for fivemers
results.fivemers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.fivemers[[i]] <- fitContinuous(phy=trees, 
                                       dat=fivemers, 
                                       model=models[i])
}

#fitmodels for sixmers
results.sixmers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.sixmers[[i]] <- fitContinuous(phy=trees, 
                                       dat=sixmers, 
                                       model=models[i])
}

#fitmodels for totalmers
results.totalmers <- list()
for(i in 1:length(models)){
  cat(paste("working on",models[i], "model\n"))
  results.totalmers[[i]] <- fitContinuous(phy=trees, 
                                       dat=totalmers, 
                                       model=models[i])
}


#make a table output with all of the relevant information
result.comp <-as.data.frame(matrix(,7,42))
colnames(result.comp)<- c("model","results.bp.Mbp", "aic","par1","par1","par1","model","results.twomers", "aic","par1","par1","par1","model","results.threemers", "aic","par1","par1","par1" ,"model", "results.fourmers", "aic","par1","par1","par1", "model","results.fivemers", "aic","par1","par1","par1", "model","results.sixmers", "aic","par1","par1","par1", "model","results.totalmers", "aic","par1","par1","par1")

#fill in bp.Mbp model data
for(i in 1:7){
  result.comp[i,1]<-models[i]
  result.comp[i,2]<-"trees"
  result.comp[i,3]<- results.bp.Mbp[[i]]$opt$aicc
  result.comp[i,4]<- results.bp.Mbp[[i]]$opt[1]
  result.comp[i,5]<- results.bp.Mbp[[i]]$opt[2]
  result.comp[i,6]<- results.bp.Mbp[[i]]$opt[3]
}

#fill in twomers model data
for(i in 1:7){
  result.comp[i,7]<-models[i]
  result.comp[i,8]<-"trees"
  result.comp[i,9]<- results.twomers[[i]]$opt$aicc
  result.comp[i,10]<- results.twomers[[i]]$opt[1]
  result.comp[i,11]<- results.twomers[[i]]$opt[2]
  result.comp[i,12]<- results.twomers[[i]]$opt[3]
}

#fill in threemers model data
for(i in 1:7){
  result.comp[i,13]<-models[i]
  result.comp[i,14]<-"trees"
  result.comp[i,15]<- results.threemers[[i]]$opt$aicc
  result.comp[i,16]<- results.threemers[[i]]$opt[1]
  result.comp[i,17]<- results.threemers[[i]]$opt[2]
  result.comp[i,18]<- results.threemers[[i]]$opt[3]
}

#fill in fourmers model data
for(i in 1:7){
  result.comp[i,19]<-models[i]
  result.comp[i,20]<-"trees"
  result.comp[i,21]<- results.fourmers[[i]]$opt$aicc
  result.comp[i,22]<- results.fourmers[[i]]$opt[1]
  result.comp[i,23]<- results.fourmers[[i]]$opt[2]
  result.comp[i,24]<- results.fourmers[[i]]$opt[3]
}

#fill in fivemers model data
for(i in 1:7){
  result.comp[i,25]<-models[i]
  result.comp[i,26]<-"trees"
  result.comp[i,27]<- results.fivemers[[i]]$opt$aicc
  result.comp[i,28]<- results.fivemers[[i]]$opt[1]
  result.comp[i,29]<- results.fivemers[[i]]$opt[2]
  result.comp[i,30]<- results.fivemers[[i]]$opt[3]
}

#fill in sixmers model data
for(i in 1:7){
  result.comp[i,31]<-models[i]
  result.comp[i,32]<-"trees"
  result.comp[i,33]<- results.sixmers[[i]]$opt$aicc
  result.comp[i,34]<- results.sixmers[[i]]$opt[1]
  result.comp[i,35]<- results.sixmers[[i]]$opt[2]
  result.comp[i,36]<- results.sixmers[[i]]$opt[3]
}

#fill in totalmers model data
for(i in 1:7){
  result.comp[i,37]<-models[i]
  result.comp[i,38]<-"trees"
  result.comp[i,39]<- results.totalmers[[i]]$opt$aicc
  result.comp[i,40]<- results.totalmers[[i]]$opt[1]
  result.comp[i,41]<- results.totalmers[[i]]$opt[2]
  result.comp[i,42]<- results.totalmers[[i]]$opt[3]
}

write.csv(result.comp, file="phylo.sig.csv")
