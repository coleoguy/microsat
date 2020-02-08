genomes <- list.files("../genomes2")
for(i in 1:length(genomes)){
  sys.com <- paste("python ~/BUSCOVM/busco3/scripts/run_BUSCO.py -i ../genomes2/",
                   genomes[i],
                   " -o result.",
                   genomes[i],
                   " -l ~/BUSCOVM/lineages/insecta_odb9/ -m geno -c 8", sep="")
  system(command=sys.com)
}
