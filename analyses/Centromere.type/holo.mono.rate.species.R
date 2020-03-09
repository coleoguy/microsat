#load in necessary packages
library(phytools)
library(geiger)

#read in insect phylogeny
tree <- read.nexus("data/tree/tree.nex")

#read in centromere data 
#row names sets the species columns as row names for the data
dat.centromere <- read.csv("data/centromere.type.csv", 
                           as.is=T,
                           row.names = 4)

# match up tree and data
#load in the centromere type
cent.type <- dat.centromere[,4]
#name the centromere type vector
names(cent.type) <- row.names(dat.centromere)
# make a vector ready for treedata
tree.drop <- list()
#loop through to drop any unmatching data or tree tips
for(i in 1:100){
  tree.drop[[i]] <- treedata(phy = tree[[i]], data=cent.type)[[1]]
}

#set up lists to store likelihood ratio tests and 
lrtest.store <- species.test <- list()
#make a vector to store the simmaps
histories <- list()
#loop through making simmaps for each of the 100 posterior distribution trees
for(i in 1:100){
  histories[[i]] <- make.simmap(tree.drop[[i]], cent.type, model="ARD", pi="estimated")
}

#read in the microsatellite data
dat.mic <- read.csv("results/micRocounter_results_TII.csv", 
                    as.is = T, row.names = 4)
#make the class of the simmaps of type Phylo
class(histories) <- "Phylo"
#loop through 
for(i in 1:100){
  print(paste("Running sample", i))
  bpMbp <- dat.mic$all/(dat.mic$gsz/1000000)
  names(bpMbp) <- row.names(dat.mic)
  species.test[[i]] <- brownieREML(tree = histories[[i]], x = bpMbp)
  #lrtest.store[[j]] <- lrtest(order.test[[j]])
}
#write a csv with the brownieREML data
write.csv(species.test, file = "brownie.data")

#create a vector for the p-values
pvals <- c()
#loop through calculating the p-value for each posterior distribution
for(i in 1:100){
  pvals[i] <- 1-pchisq(2 * (-species.test[[i]]$logL1 + species.test[[i]]$logL.multiple),1)}
#write a csv with the p-values
write.csv(pvals, "pvals.holo.mono.new")
#STORE A FILE WITH ALL PVALUES AND RATE ESTIMATES
#FOR TIME BEING MAYBE STORE YOUR ENVIRONMENT AS WELL

#create vectors for the monocentric and holocentric data
rates.mon <- rates.hol <- c()
#loop through testing each of the rates for holocentric and monocentric species
for(i in 1:100){
  rates.mon[i] <- species.test[[i]]$sig2.multiple[2]
  rates.hol[i] <- species.test[[i]]$sig2.multiple[1]
}
#write csv files for both monocentric and holocentric values
write.csv(rates.mon, "rates.mon.new")
write.csv(rates.hol, "rates.hol.new")

#read in the monocentric and holocentric data
rates <- read.csv("rate.holo.mono.new.csv", as.is = T)
#adjust the format of the rates column
rates$Rates <- as.numeric(gsub(",", "", rates$Rates, fixed = T))
#change the unit of rates to mbp/million years
rates$Rates <- rates$Rates/1000000
#create a boxplot with the data
boxplot(Rates~Type.of.Centromere, 
        data = rates, 
        xlab = "Type of Centromere",
        ylab = 'Rate of Microsatellite Evolution (mbp/my)',
        outline =F)
#add the inidividual points to the boxplot
stripchart(rates$Rates ~ rates$Type.of.Centromere, vertical = TRUE, data = rates, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
#save as pdf 5" x 6"
#ggplot a density plot of the data as well
library(ggplot2)
ggplot(data = rates, aes(x=Rates)) + geom_density(aes(color=Type.of.Centromere))
#save as pdf 6" x 4"

