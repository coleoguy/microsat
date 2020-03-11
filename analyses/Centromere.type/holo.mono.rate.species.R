#load in necessary packages
library(phytools)
library(geiger)

#read in insect phylogeny
trees <- read.nexus("../../data/trees/post.nex")

#read in centromere data 
#row names sets the species columns as row names for the data
dat.centromere <- read.csv("../../data/traits/centromere.type.csv", 
                           as.is=T,
                           row.names = 4)

# converting data to a named vector
#load in the centromere type
cent.type <- dat.centromere[,4]
#name the centromere type vector
names(cent.type) <- row.names(dat.centromere)
rm(dat.centromere)

# match up tree and data

trees.pruned <- list()
#loop through to drop any unmatching data or tree tips
for(i in 1:100){
  trees.pruned[[i]] <- treedata(phy = trees[[i]], data=cent.type)[[1]]
}
rm(trees, i)

# remove "Timema_cristinae" this species is not place
# in the tree correctly

foo <- list()
for(i in 1:100){
  foo[[i]] <- drop.tip(trees.pruned[[i]], tip="Timema_cristinae")
}
trees.pruned <- foo
rm(foo)
# also remove this data from the data
cent.type <- cent.type[names(cent.type)!="Timema_cristinae"]
#make a vector to store the simmaps
histories <- list()

#loop through making simmaps for each of the 100 trees
for(i in 1:100){
  print(paste("working on tree", i))
  histories[[i]] <- make.simmap(trees.pruned[[i]], 
                                cent.type, 
                                model="ARD", 
                                pi="estimated")
}

#make the class of the simmaps of type Phylo
class(histories) <- "simmap"


#read in the microsatellite data
dat.mic <- read.csv("../../results/ssr.inference/micRocounter_results_TII.csv", 
                    as.is = T, row.names = 4)

bpMbp <- dat.mic$all/(dat.mic$gsz/1000)
names(bpMbp) <- row.names(dat.mic)

# setup dataframe for results
results <- as.data.frame(matrix(NA,100, 4))
colnames(results) <- c("pval","rate.mon","rate.hol","conv")

#loop through 
for(i in 1:100){
  print(paste("Running sample", i))
  working <- T
  while(working){
    brownie.fit <- brownieREML(tree = histories[[i]], 
                               x = bpMbp, maxit = 100000)
    if(brownie.fit$convergence=="Optimization has converged."){
      working <- F
    }
  }
  results$pval[i] <- 1-pchisq(2 * (-brownie.fit$logL1 + 
                                     brownie.fit$logL.multiple),1)
  results$rate.hol[i] <- brownie.fit$sig2.multiple[1]
  results$rate.mon[i] <- brownie.fit$sig2.multiple[2]
  results$conv[i] <- brownie.fit$convergence
}

#write a csv with the brownieREML data
write.csv(results, file = "../../results/cent.vs.rate.csv", row.names = F)





