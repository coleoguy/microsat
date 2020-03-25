#set figures as working directory
#load in libraries
library(geiger)
library(phylolm)
library(rr2)

#read in nmicrosatellite data
str <- read.csv("../data/traits/micro.vs.chrom.csv")
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

#read in trees
trees <- read.nexus("../data/trees/post.nex")

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
pvals.content <- beta.content <- c()
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
  R2(fit, simp.mod, phy = tree.cur)
  cur.results[[16]]
  pvals.content[i] <- cur.results$coefficients[2,6]
  beta.content[i] <- cur.results$coefficients[2,1]
}
#make the results into a data frame
results <- data.frame(pvals.content,beta.content)
#write the results to a csv
write.csv(pvals.content, "../results/genome.size/gsz.content.csv")



#put gsz and microsatellite data into different units
str$all/1000000 ->msat
#plot the microsatellite content and genome size
plot(msat~str$gsz,
     ylab = "Microsatellite Content (Mbp)",
     xlab = "Genome Size (Mbp)",
     pch = 16,
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

#save as pdf 4.3" x 4.3"
