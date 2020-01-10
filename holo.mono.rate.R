#setwd as r scripts

#read in phytools
library(phytools)

#read in the insect phylogeny
tree <- read.nexus("../data/tree/tree.nex")

# read microsat content in with the species data
dat.mic <- read.csv("../results/bp.results.good.csv", as.is=T)

#creates a vector of unique orders
orders <- unique(dat.mic$order)

#organizes the orders into the same order they should appear on the tree
new.orders <- orders[c(11, 12, 1, 14, 9, 4, 15, 2, 8, 13, 7, 6, 16, 10, 5, 3)]

#creates a vector with matching tree tips
keep <- c("Blaberus", "Tribolium", "Drosophila", "Campodea","Ephemera", "Aphis",
          "Acyrthosiphon", "Bombus", "Bombyx", "Calopteryx", "Tetrix", "Timema",
          "Menopon", "Perla", "Mengenilla", "Frankliniella")

#makes a new tree with only the orders that are on the plot to 
#match the tree tips
new.tree <- keep.tip(tree, keep)

#label the tips with the orders instead of the genus information
orders <- c("Hemiptera", "Homoptera", "Lepidoptera", "Blattodea",
            "Hymenoptera", "Odonata", "Diplura", "Diptera", 
            "Ephemeroptera","Thysanoptera", "Strepsiptera", "Phthiraptera",
            "Plecoptera", "Coleoptera", "Orthoptera", "Phasmatodea")
new.tree$tip.label<-orders

#plots the new tree
plot(new.tree)

#make a vector with the holocentric or metacentric designation for use in the 
#simmap
chrom.type <- c("H", "H", "H", "M", "M", "H", "M", "H", "H", "H", "M", "M", "M")
names(chrom.type) <- c("Hemiptera", "Homoptera", "Lepidoptera", "Blattodea", 
                       "Hymenoptera", "Odonata", "Diptera", "Ephemeroptera", 
                       "Thysanoptera", "Phthiraptera", "Coleoptera", "Orthoptera",
                       "Phasmatodea")
#drop the tips where we do not know 
drop <- c("Plecoptera", "Strepsiptera", "Diplura")
drop.tree <- drop.tip(new.tree, drop)
plot(drop.tree)

#make.simmap for use in the brownieREML function

orders.new <- c("Hemiptera", "Homoptera", "Lepidoptera", "Blattodea",
            "Hymenoptera", "Odonata", "Diptera", 
            "Ephemeroptera","Thysanoptera", "Phthiraptera",
            "Coleoptera", "Orthoptera", "Phasmatodea")

lrtest.store <- order.test <- trees <- list()
for(i in 1:100){
  trees[[i]] <- make.simmap(drop.tree, chrom.type, model="ARD", pi="estimated")
  
}


for(j in 1:100){
  print(paste("Running sample", j))
  order.data <- c()
  bpMbp <- dat.mic$all/(dat.mic$gsz/1000000)
  for(i in 1:length(orders.new)){
    order.data[i] <- sample(bpMbp[which(dat.mic$order == orders.new[i])], 1)
  }
  names(order.data) <- orders.new 
  #order.data<-order.data/1000000
  order.test[[j]] <- brownieREML(tree = trees[[j]], x = order.data)
  #lrtest.store[[j]] <- lrtest(order.test[[j]])
}
pvals <- c()
for(i in 1:100){
 pvals[i] <- 1-pchisq(2 * (-order.test[[i]]$logL1 + order.test[[i]]$logL.multiple),1)
}


#WITHOUT PTHIRAPTERA

#make a vector with the holocentric or metacentric designation for use in the 
#simmap
chrom.type <- c("H", "H", "H", "M", "M", "H", "M", "H", "H", "M", "M", "M")
names(chrom.type) <- c("Hemiptera", "Homoptera", "Lepidoptera", "Blattodea", 
                       "Hymenoptera", "Odonata", "Diptera", "Ephemeroptera", 
                       "Thysanoptera", "Coleoptera", "Orthoptera",
                       "Phasmatodea")
#drop the tips where we do not know 
drop <- c("Plecoptera", "Strepsiptera", "Diplura", "Phthiraptera")
drop.tree <- drop.tip(new.tree, drop)
plot(drop.tree)

#make.simmap for use in the brownieREML function

orders.new <- c("Hemiptera", "Homoptera", "Lepidoptera", "Blattodea",
                "Hymenoptera", "Odonata", "Diptera", 
                "Ephemeroptera","Thysanoptera",
                "Coleoptera", "Orthoptera", "Phasmatodea")

lrtest.store <- order.test <- trees <- list()
for(i in 1:100){
  trees[[i]] <- make.simmap(drop.tree, chrom.type, model="ARD", pi="estimated")
  
}


for(j in 1:100){
  print(paste("Running sample", j))
  order.data <- c()
  bpMbp <- dat.mic$all/(dat.mic$gsz/1000000)
  for(i in 1:length(orders.new)){
    order.data[i] <- sample(bpMbp[which(dat.mic$order == orders.new[i])], 1)
  }
  names(order.data) <- orders.new 
  order.test[[j]] <- brownieREML(tree = trees[[j]], x = order.data)
}
pvals <- c()
for(i in 1:100){
  pvals[i] <- 1-pchisq(2 * (-order.test[[i]]$logL1 + order.test[[i]]$logL.multiple),1)
}

rates.mon <- rates.hol <- c()
for(i in 1:100){
  rates.mon[i] <- order.test[[i]]$sig2.multiple[2]
  rates.hol[i] <- order.test[[i]]$sig2.multiple[1]
}


plot(x=rates.mon,y=rates.hol, pch=16, col=rgb(1,0,0,.5),
     xlim=c(0,34000),ylim=c(0,34000),
     xlab="rate for monocentric",
     ylab="rate for holocentric")
lines(x=c(0,34000),y=c(0,34000),lty=3)



