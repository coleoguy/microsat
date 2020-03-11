library(phytools)
library(geiger) 
library(ape)
library(BAMMtools)
library(coda)
library(evobiR)

# read in sampled trees
trees <- read.nexus("../data/trees/post.nex")

# select a tree at random for making a continuous trait map
tree.BAMM1 <- trees[[1]]

# read in  microsatellite data 
dat.microsat <- read.csv("../results/ssr.inference/micRocounter_results_TII.csv",
                         as.is=T,
                         row.names = 4)

# fill in the column with all microsatellite content from dat.microsat
totalmers <- dat.microsat[,9]
# adds names to the bpMbp vector from dat.microsat
names(totalmers) <- row.names(dat.microsat)

#drop the tips
tree.pruned <- drop.tip(phy=tree.BAMM1, tip= c("B.terrestris", "Plutella_xylostella")) 

#write the new tree into the BAMM folder
write.tree(tree.pruned, file = "tree(1).new")

#load in tree
tree.BAMM1 <- read.tree("tree(1).new")
plot(tree.BAMM1)

#tree needs to be ultrametric and binary
is.ultrametric(tree.BAMM)
is.binary(tree.BAMM)

#does not work because data is not in correct format
setBAMMpriors(phy = tree.BAMM1, traits = totalmers)

#sets the priors needed for the output file
#you then run the command file in windows bamm -c BAMM.bp.Mbp.controlfile.txt

#ran for 10,000,000 generations

#plot the loglike to see if it's reached convergence
mcmcout <- read.csv("totalmers.microsat.mcmc_out.txt", header =T)
plot(mcmcout$logLik ~ mcmcout$generation)

#print the plot
tiff("loglikoutput.tiff", height = 6, width = 8, units = 'in', res = 300)
plot(mcmcout$logLik ~ mcmcout$generation)
dev.off()

#specify 10% burnin
burnstart <- floor(0.1*nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout),]

#came up with effective size of 735.2577
effectiveSize(postburn$N_shifts)
#came up with effective size of 415.7432
effectiveSize(postburn$logLik)

edata <- getEventData(tree.BAMM1, eventdata = "event_data.txt", burnin = 0.1, type = "trait")

post_probs <- table(postburn$N_shifts)/nrow(postburn)

names(post_probs)

shift_probs <- summary(edata)

#makes a credible shift set; calculates the expected number of shifts to be most
#likely
css <- credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold = 15, set.limit = 0.95)
css$number.distinct


summary(css)

write.table(summary(css), "microsatsummary.txt")

plot.credibleshiftset(css)
tiff("credibleshiftset.tiff", height = 12, width = 12, units = 'in', res = 300)
plot.credibleshiftset(css)
dev.off()

marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex=0.3)



tiff("marginalshifttree.tiff", height = 6, width = 6, units = 'in', res = 300)
plot.phylo(marg_probs, cex = 0.25)
dev.off()

plotRateThroughTime(edata, xlim = "auto", ylim = c(-1, 1000))
plotRateThroughTime(edata, xlim = "auto", ylim = c(-10, 30))

plotRateThroughTime(edata, xlim = "auto", interval = "blue", 
                    avgCol = "blue", node = XXX, ylim = (0, 50))