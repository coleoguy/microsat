# Read data
files <- list.files("../data/busco", full.names = F)

#Pulls the species name for each file to create a vector of row names
species<- c()
for(i in 1:length(files)){
  X <- strsplit(files[i],split=".",fixed=T)[[1]][2]
  species[i] <- paste(X[1], "_microsat.txt", collapse="",sep="")
}

#Creates the matrix with the column names
results.busco <- as.data.frame(matrix(, length(files), 9))
colnames(results.busco) <- c("species", t(read.table(paste("../data/busco/", files[1], sep=""), 
                            skip = 9, sep="\t", as.is=T,
                            check.names = F, header = F)[,3]), "score", "good.qual")

#Fills in the matrix row by row for each species
for(i in 1:length(files)){
  results.busco[i, 1] <- species[i]
  results.busco[i, 2:7] <- read.table(paste("../data/busco/", 
                                         files[i], 
                                         sep=""),
                                   skip = 9, 
                                   sep="\t", 
                                   as.is=T,
                                   check.names = F, 
                                   header = F)[,2]
#This calculates the BUSCO score for each species  
    results.busco$score[i] <- results.busco$`Complete BUSCOs (C)`[i] / results.busco$`Total BUSCO groups searched`[i]
#This makes a column "good.qual" which will indicate based on a 90% cutoff threshold for BUSCO scores, which 
#species meet the criteria necessary to move on with further analyses (<90% is false, >90% is true)
  if(results.busco$score[i] < .90) results.busco$good.qual[i] <- F
  if(results.busco$score[i] >= .90) results.busco$good.qual[i] <- T
}

#Lets yo uknow how many species need to be removed from downstream analyses based on the cutoff of 90% for the 
#BUSCO score
sum(results.busco$good.qual == FALSE)

#identify which species are those with quality under the threshold set at 90%
results.busco$good.qual == FALSE

#Writes CSV file to store data
write.csv(results.busco, file="../results/busco.results.csv", row.names = F)



