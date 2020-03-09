# READ IN ALL DATA NEEDED FOR ANALYSES

# set wd as figures
# load phytoools
library(phytools)

# read in sampled trees
trees <- read.nexus("data/tree/tree.nex")

# select a tree at random for making a continuous trait map
tree <- trees[[sample(1:100, 1)]]

# read in  microsatellite data 
dat.microsat <- read.csv("results/micRocounter_results_TII.csv",
                         as.is=T,
                         row.names = 4)

# drops the tip for B.terrestis
pruned.tree <- drop.tip(phy=tree, tip="B.terrestris")
pruned.tree <- drop.tip(phy=pruned.tree, tip="Plutella_xylostella")


# PREPARE DATA FOR CONT MAPS


# fill in the bpMbp with the correct column of data from dat.microsat
mscontbpMbp <- dat.microsat[,18]
# adds names to the bpMbp vector from dat.microsat
names(mscontbpMbp) <- row.names(dat.microsat)

# make vectors with the same format as bpMbp from above
mscontbp2 <- mscontbp3 <- mscontbp4 <- mscontbp5 <- mscontbp6 <- mscontbpall <- mscontbpMbp

# fill in the 2-6mers and all microsats with the correct column 
# of data from dat.microsat
mscontbp2 <-  dat.microsat[,4]
# fill in the 3mers with the correct column of data from dat.microsat
mscontbp3 <-  dat.microsat[,5]
# fill in the 4mers with the correct column of data from dat.microsat
mscontbp4 <-  dat.microsat[,6]
# fill in the 5mers with the correct column of data from dat.microsat
mscontbp5 <-  dat.microsat[,7]
# fill in the 6mers with the correct column of data from dat.microsat
mscontbp6 <-  dat.microsat[,8]
# fill in the column with all microsatellite content from dat.microsat
mscontbpall <-  dat.microsat[,9]

# adds names to the vectors for 2-6mer and all microsatellite content
names(mscontbp2) <- names(mscontbp3) <- names(mscontbp4) <- names(mscontbp5) <- names(mscontbp6) <- names(mscontbpall) <- row.names(dat.microsat)


# MAKE CONTMAPS FOR EACH TYPE OF MICROSATELLITE


# contmap bpMbp
anc.con <- contMap(pruned.tree, round(mscontbpMbp), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan")

# export as PDF 10" x 15"

#set up area for contmaps
par(mfrow = c(2,3))

# contmap 2mer
anc.con <- contMap(pruned.tree, round(mscontbp2), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan", legend = F)

# contmap 3mer
anc.con <- contMap(pruned.tree, round(mscontbp3), type="fan", plot = F)
plot(anc.con, legend=F, fsize=0.01,type="fan", legend = F)

# contmap 4mer
anc.con <- contMap(pruned.tree, round(mscontbp4), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan", legend = F)

# contmap 5mer
anc.con <- contMap(pruned.tree, round(mscontbp5), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan", legend = F)

# contmap 6mer
anc.con <- contMap(pruned.tree, round(mscontbp6), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan", legend = F)

#contmap all microsatellites
anc.con <- contMap(pruned.tree, round(mscontbpall), type="fan", plot = F)
plot(anc.con, ftype="off", fsize=c(.01,.5), type="fan", legend = F)

# export as PDF 10" x 12"
