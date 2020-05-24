#load in necessary packages
library(geiger)


#read in insect phylogeny
trees <- read.nexus("../../data/trees/post.nex")

#read in the microsatellite data
dat.mic <- read.csv("../../results/ssr.inference/micRocounter_results_TII.csv",
                    as.is = T, row.names = 4)

#subset order data
coleoptera <- dat.mic[which(dat.mic$order == "Coleoptera"),]
diptera <- dat.mic[which(dat.mic$order == "Diptera"),]
hemiptera <- dat.mic[which(dat.mic$order %in% c("Homoptera","Hemiptera")),]
hymenoptera <- dat.mic[which(dat.mic$order == "Hymenoptera"),]
lepidoptera <- dat.mic[which(dat.mic$order == "Lepidoptera"),]

#loop through to create coleoptera trees with matching data and tips
trees.coleoptera <- c()
for(i in 1:100){
  trees.coleoptera[[i]] <- treedata(phy = trees[[i]], data=coleoptera)[[1]]
}

#loop through to create diptera trees with matching data and tips
trees.diptera <- c()
for(i in 1:100){
  trees.diptera[[i]] <- treedata(phy = trees[[i]], data=diptera)[[1]]
}

#loop through to create hemiptera trees with matching data and tips
trees.hemiptera <- c()
for(i in 1:100){
  trees.hemiptera[[i]] <- treedata(phy = trees[[i]], data=hemiptera)[[1]]
}

#loop through to create hymneoptera trees with matching data and tips
trees.hymenoptera <- c()
for(i in 1:100){
  trees.hymenoptera[[i]] <- treedata(phy = trees[[i]], data=hymenoptera)[[1]]
}

#loop through to create lepidoptera trees with matching data and tips
trees.lepidoptera <- c()
for(i in 1:100){
  trees.lepidoptera[[i]] <- treedata(phy = trees[[i]], data=lepidoptera)[[1]]
}


ace.diptera <- ace.coleoptera <- ace.hemiptera <-
  ace.hymenoptera <- ace.lepidoptera <- list()

#ancestral character estimations for coleoptera
mics.coleoptera <- coleoptera$bp.Mbp
names(mics.coleoptera) <- row.names(coleoptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.coleoptera[[i]] <- fitContinuous(phy=trees.coleoptera[[i]],
                                       dat=mics.coleoptera, ncores=14)
}
for(i in 1:100){
  print(sum(ace.coleoptera[[i]]$res[,3]))
}

coleoptera.rates <-c()
for(i in 1:100){
  coleoptera.rates[i] <- ace.coleoptera[[i]]$opt$sigsq
}

#ancestral character estimations for diptera
mics.diptera <- diptera$bp.Mbp
names(mics.diptera) <- row.names(diptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.diptera[[i]] <- fitContinuous(phy=trees.diptera[[i]],
                                    dat=mics.diptera, ncores = 14)
}
diptera.rates <-c()
for(i in 1:100){
  diptera.rates[i] <- ace.diptera[[i]]$opt$sigsq
}


#ancestral character estimations for hemiptera
mics.hemiptera <- hemiptera$bp.Mbp
names(mics.hemiptera) <- row.names(hemiptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.hemiptera[[i]] <- fitContinuous(phy=trees.hemiptera[[i]],
                                      dat=mics.hemiptera, ncores = 14)
}
hemiptera.rates <-c()
for(i in 1:100){
  hemiptera.rates[i] <- ace.hemiptera[[i]]$opt$sigsq
}

#ancestral character estimations for hymenoptera
mics.hymenoptera <- hymenoptera$bp.Mbp
names(mics.hymenoptera) <- row.names(hymenoptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.hymenoptera[[i]] <- fitContinuous(phy=trees.hymenoptera[[i]],
                                        dat=mics.hymenoptera, ncores = 14)
}
hymenoptera.rates <- c()
for(i in 1:100){
  hymenoptera.rates[i] <- ace.hymenoptera[[i]]$opt$sigsq
}

#ancestral character estimations for lepidoptera
mics.lepidoptera <- lepidoptera$bp.Mbp
names(mics.lepidoptera) <- row.names(lepidoptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.lepidoptera[[i]] <- fitContinuous(phy=trees.lepidoptera[[i]],
                                        dat=mics.lepidoptera, ncores = 14)
}
lepidoptera.rates <-c()
for(i in 1:100){
  lepidoptera.rates[i] <- ace.lepidoptera[[i]]$opt$sigsq
}

#make a data frame with all the order data
order.rates <- data.frame(coleoptera.rates,diptera.rates,
                          hemiptera.rates, hymenoptera.rates,
                          lepidoptera.rates)
#store in a file
write.csv(order.rates, file = "../../results/order.rates/order.rates.csv")
