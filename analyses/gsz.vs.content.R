#set figures as working directory
#load in libraries
library(geiger)
library(phylolm)
library(rr2)

#read in nmicrosatellite data
str <- read.csv("../data/micro.vs.chrom.csv")
str <- str[,1:3]
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

#write clean data to a csv
write.csv(str, "../data/str.content")

#read in trees
trees <- read.nexus("../data/post.nex")

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
pvals.content <- beta.content <- intercept <- c()
str$gsz <- str$gsz/1000000
prop.var.exp <- c()
for(i in 1:100){
  print(i)
  #stores tree number
  tree.test <- trees[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=str)
  #stores current trees data
  tree.cur <- foo[[1]]
  #stores p-value on phylolm analysis
  fit <- phylolm(all ~ gsz,
                 data = str,
                 phy = tree.cur,
                 model = "BM",
                 boot = 100)
  simp.mod <- lm(str$all ~ str$gsz)
  cur.results <- summary(fit)
  prop.var.exp[i] <- R2(fit, simp.mod, phy = tree.cur)[3]
  pvals.content[i] <- cur.results$coefficients[2,6]
  beta.content[i] <- cur.results$coefficients[2,1]
  intercept[i] <- cur.results$coefficients[1,1]
}

#make the results into a data frame
results <- data.frame(pvals.content,beta.content, prop.var.exp, intercept)

#write the results to a csv
write.csv(results, "../results/gsz.content.csv")
