library(geiger)

# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/ssr.inference/micRocounter_results_TII_typecentromere.csv", 
                          row.names = 4)

# run aovphylo with phylogenetic correction
# make named vector for bpMbp coontent
bp.Mbp <- microsat.cent$bp.Mbp
names(bp.Mbp) <- row.names(microsat.cent)

# make named vector for all microsat content
bp.all <- microsat.cent$all
names(bp.all) <- row.names(microsat.cent)

# plot for presentation
boxplot(log(bp.all) ~ holo.or.mono,
        data = microsat.cent,
        outpch = NA,
        xlab = "Type of Centromere",
        ylab = "log Microsatellite Content (bp)")
stripchart(log(microsat.cent$all) ~ microsat.cent$holo.or.mono, 
           vertical = TRUE, 
           data = microsat.cent, 
           method = "jitter", 
           add = TRUE, 
           pch = 20, 
           col = rgb(250, 159, 181, 100,
                     maxColorValue = 255))


# export pdf at 4.3" x 4.3"

holo.or.mono <- microsat.cent$holo.or.mono
names(holo.or.mono) <- row.names(microsat.cent)

#make named vector for 2mer content
bp.2 <- microsat.cent$twomers
names(bp.2) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp2 <- aov.phylo(bp.2 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

#make named vector for 3mer content
bp.3 <- microsat.cent$threemers
names(bp.3) <- row.names(microsat.cent)
# run phyloANOVA for 3mers and centromere type
aovphylo.bp3 <- aov.phylo(bp.3 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 4mer content
bp.4 <- microsat.cent$fourmers
names(bp.4) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp4 <- aov.phylo(bp.4 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 5mer content
bp.5 <- microsat.cent$fivemers
names(bp.5) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp5 <- aov.phylo(bp.5 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 6mer content
bp.6 <- microsat.cent$sixmers
names(bp.6) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp2 <- aov.phylo(bp.6 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)
