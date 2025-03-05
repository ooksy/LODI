# data
mouse <- read.csv("microbiome_modified.csv")
mouse$Genotype <- factor(mouse$Genotype, levels = c("WT", "ACKG", "ACK", "HACKG", "HACK"))
mouse$Grouped.genotype <- factor(mouse$Grouped.genotype, levels = c("WT", "ACK/G", "HACK/G"), ordered = T)
mouse$Timepoint <- factor(mouse$Timepoint, levels = c("W0", "W1", "W3", "W7", "W11", "W15", "W19"))

uniq_grp <- unique(mouse$Grouped.genotype)
uniq_timept <- unique(mouse$Timepoint)

par(mfrow = c(1,length(uniq_timept)))
for(i in 1:length(uniq_grp)){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == uniq_grp[i]) & (mouse$Timepoint == uniq_timept[j])
    vec_dat <- unname(unlist(mouse[idx ,6:59]))
    plot(density(vec_dat), main = paste("grp :", uniq_grp[i], "timept", uniq_timept[j]))
  }
}


