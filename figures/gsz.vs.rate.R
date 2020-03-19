library(geiger)
library(phylolm)

#read in csv with rates of evolution
rates <- read.csv("../analyses/tip.rates/tip.rates.csv",
                  row.names = 1)

#store the average rate in a named vector by species name
rates.species <- rates$Average
names(rates.species) <- row.names(rates)


#read in nmicrosatellite data
str <- read.csv("../data/traits/micro.vs.chrom.csv")
str <- str[,1:3]
str$rates <- NA
str$species <- gsub(pattern = " ", replacement = "_", str$species)

#fix the names in the data that are different from the trees
str$species <- sub(pattern = "Harpergnathos_saltator", 
                   replacement = "Harpegnathos_saltator", str$species)
str$species <- gsub(pattern = "helicoverpa_zea", 
                    replacement = "Helicoverpa_zea", str$species)
str$species <- gsub(pattern = "papilio_polytes", 
                    replacement = "Papilio_polytes", str$species)
str$species <- gsub(pattern = "Pogomyrmex_barbatus", 
                   replacement = "Pogonomyrmex_barbatus", str$species)
str$species <- gsub(pattern = "Scaptodrosophila_lebanonesis", 
                    replacement = "Scaptodrosophila_lebanonensis", str$species)



row.names(str) <- str$species
for(i in 1:nrow(str)){
    hit <- which(names(rates.species) == str$species[i])[1]
    #fill in rates for those species that have a match in species column 
    str$rates[i] <- rates.species[hit]
}

#clean environment
rm(list = c("rates", "hit", "i", "rates.species"))

#read in trees
trees <- read.nexus("../data/trees/post.nex")


#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- trees[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=str)
  #stores current trees data
  tree.cur <- foo[[1]]
  #stores p-value on phylolm analysis
  pvals.rates[i] <- summary(phylolm(gsz ~ all, 
                                    data = str, 
                                    phy = tree.cur, 
                                    model = "BM",
                                    boot = 100))$coefficients[2,6]
}

#plot the microsatellite evolution rates and genome size
plot(str$rates~str$gsz,
     xlab = "Genome Size",
     ylab = "Microsatellite Evolution Rates",
     pch = 16, 
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

#plot the microsatellite content and genome size
plot(str$gsz~str$all,
     xlab = "Microsatellite Evolution Rates",
     ylab = "Microsatellite Content (bp)",
     pch = 16, 
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))
