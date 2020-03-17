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
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.coleoptera <- ace(x=coleoptera$bp.Mbp, 
                        phy=trees.coleoptera[[i]], 
                        type = "continuous",
                        model="BM",
                        method="REML")
}

#ancestral character estimations for diptera
mics <- diptera$bp.Mbp
names(mics) <- row.names(diptera)
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.diptera[[i]] <- fitContinuous(phy=trees.diptera[[i]], dat=mics/1000, ncores=10)
}
dip.rates <-c()
for(i in 1:100){
  dip.rates[i] <- ace.diptera[[1]]$opt$sigsq
}


#ancestral character estimations for hemiptera
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.hemiptera <- ace(x=hemiptera$bp.Mbp, 
                       phy=trees.hemiptera[[i]],
                       type = "continuous",
                       model="BM",
                       method= "REML")
}

#ancestral character estimations for hymenoptera
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.hymenoptera <- ace(x=hymenoptera$bp.Mbp, 
                         phy=trees.hymenoptera[[i]], 
                         type="continuous",
                         model="BM",
                         method="REML")
}

#ancestral character estimations for lepidoptera
for(i in 1:100){
  print(i)
  # estimate ancestral states
  ace.lepidoptera <- ace(x=lepidoptera$bp.Mbp, 
                         phy=trees.lepidoptera[[i]],
                         type="continuous",
                         model="BM",
                         method="REML")
}



