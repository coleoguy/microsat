#set figrues as working directory
#load in libraries
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


#assign rownames
row.names(str) <- str$species
#fill in the rates into the other data frame
for(i in 1:nrow(str)){
    hit <- which(names(rates.species) == str$species[i])[1]
    #fill in rates for those species that have a match in species column
    str[i, 4:104] <- rates[hit, 1:100]
}

#clean environment
rm(list = c("rates", "hit", "i", "rates.species"))

#read in trees
trees <- read.nexus("../data/trees/post.nex")


#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
pvals.rates <- beta.rates <- intercept <- c()
str$gsz <- str$gsz/1000000
for(i in 1:100){
  print(i)
  #stores tree number
  tree.test <- trees[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=str)
  #stores current trees data
  tree.cur <- foo[[1]]
  #stores p-value on phylolm analysis
  cur.results <- summary(phylolm(str[,(i+3)] ~ gsz,
                                    data = str,
                                    phy = tree.cur,
                                    model = "BM",
                                    boot = 100))
  pvals.rates[i] <- cur.results$coefficients[2,6]
  beta.rates[i] <- cur.results$coefficients[2,1]
  intercept[i] <- cur.results$coefficients[1,1]

}
#make the results into a data frame
results <- data.frame(pvals.rates,beta.rates, intercept)
#write the results to a csv
write.csv(pvals.rates, "../results/genome.size/gsz.rates.csv")