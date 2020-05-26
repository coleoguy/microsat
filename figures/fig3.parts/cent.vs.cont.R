# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/micRocounter_results_TII_typecentromere.csv",
                          row.names = 4)

# plot for presentation
boxplot(log(all) ~ holo.or.mono,
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
