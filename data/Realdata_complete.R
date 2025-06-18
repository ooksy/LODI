# data
mouse <- read.csv("data/microbiome_modified.csv")
mouse$Genotype <- factor(mouse$Genotype, levels = c("WT", "ACKG", "ACK", "HACKG", "HACK"))
mouse$Grouped.genotype <- factor(mouse$Grouped.genotype, levels = c("WT", "ACK/G", "HACK/G"), ordered = T)
mouse$Timepoint <- factor(mouse$Timepoint, levels = c("W0", "W1", "W3", "W7", "W11", "W15", "W19"))

uniq_grp <- unique(mouse$Grouped.genotype)
uniq_grp3 <- unique(mouse$Grouped.genotype)
uniq_grp5 <- unique(mouse$Genotype)
uniq_timept <- unique(mouse$Timepoint)

View(mouse)

# Make true data list
# 같은 mouse가 같은 자리에 있어야함? - time point 마다 샘플들이 다르다
Data <- mouse[, 6:59]
Data_true <- list()

for(i in uniq_timept){
  Data_true[[i]] <- Data[mouse$Timepoint == i, ]
}

Data_true


# density plot aloOTU1# density plot aloOTU1# density plot along the timepoint : real data
par(mfrow = c(1,length(uniq_timept)))
for(i in 1:length(uniq_grp)){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == uniq_grp[i]) & (mouse$Timepoint == uniq_timept[j])
    vec_dat <- unname(unlist(mouse[idx ,6:59]))
    plot(density(vec_dat), main = paste("grp :", uniq_grp[i], "timept", uniq_timept[j]))
  }
}



# number of sample size for each category
for(i in 1:length(uniq_grp5)){
  for(j in 1:length(uniq_timept)){
    n <- sum(mouse$Genotype == levels(uniq_grp5)[i] & mouse$Timepoint == uniq_timept[j])
    cat("sample size of genotype", levels(uniq_grp5)[i],"and timepoint", levels(uniq_timept)[j], "is", n,"\n")
  }
}


# number of sample size for each category
for(i in 1:length(uniq_grp3)){
  for(j in 1:length(uniq_timept)){
    n <- sum(mouse$Grouped.genotype == levels(uniq_grp3)[i] & mouse$Timepoint == uniq_timept[j])
    cat("sample size of genotype", levels(uniq_grp3)[i],"and timepoint", levels(uniq_timept)[j], "is", n,"\n")
  }
}


# summary statistics for each category
# location parameters
for(i in 1:length(uniq_grp3)){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == levels(uniq_grp3)[i]) & (mouse$Timepoint == uniq_timept[j])
    otus_grp <- unname(unlist(mouse[idx, 6:59]))
    res <- summary(otus_grp)
    cat("summary stat of genotype", levels(uniq_grp3)[i],"and timepoint", levels(uniq_timept)[j], "is", res,"\n")
  }
}

# sd
for(i in 1:length(uniq_grp3)){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == levels(uniq_grp3)[i]) & (mouse$Timepoint == uniq_timept[j])
    otus_grp <- unname(unlist(mouse[idx, 6:59]))
    res <- sd(otus_grp)
    cat("summary stat of genotype", levels(uniq_grp3)[i],"and timepoint", levels(uniq_timept)[j], "is", res, "\n")
  }
}



par(mfrow = c(1,1))
# density - transform X
for(i in 1:length(uniq_grp3)){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == levels(uniq_grp3)[i]) & (mouse$Timepoint == uniq_timept[j])
    otus_grp <- unname(unlist(mouse[idx, 6:59]))
    
    png(file = paste("C:/Users/SEC/Desktop/research/eda/density/density",ACKG,uniq_timept[j],".png", sep=""),
        width = 11,
        height = 9,
        units = "in",
        res = 1200)
    
    plot(density(otus_grp), main = paste("group :", levels(uniq_grp3)[i], ", timepoint :", levels(uniq_timept)[j]))
    
    dev.off()
  }
}

levels(uniq_grp3) <- c("WT", "ACK/G", "HACK/G")

par(mfrow = c(1,1))
# density - transform한 값으로 해볼까?
for(i in 2){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == levels(uniq_grp3)[i]) & (mouse$Timepoint == uniq_timept[j])
    otus_grp <- unname(unlist(mouse[idx, 6:59]))
    
    png(file = paste("C:/Users/SEC/Desktop/research/eda/density/densityACKG",uniq_timept[j],".png", sep=""),
        width = 11,
        height = 9,
        units = "in",
        res = 1200)
    
    plot(density(otus_grp), main = paste("group :", levels(uniq_grp3)[i], ", timepoint :", levels(uniq_timept)[j]))
    
    dev.off()
  }
}

i=2;j=1

par(mfrow = c(1,1))
# density - transform한 값으로 해볼까?
for(i in 3){
  for(j in 1:length(uniq_timept)){
    idx <- (mouse$Grouped.genotype == levels(uniq_grp3)[i]) & (mouse$Timepoint == uniq_timept[j])
    otus_grp <- unname(unlist(mouse[idx, 6:59]))
    
    png(file = paste("C:/Users/SEC/Desktop/research/eda/density/densityHACKG",uniq_timept[j],".png", sep=""),
        width = 11,
        height = 9,
        units = "in",
        res = 1200)
    
    plot(density(otus_grp), main = paste("group :", levels(uniq_grp3)[i], ", timepoint :", levels(uniq_timept)[j]))
    
    dev.off()
  }
}

graphics.off()


# covariance/correlation along the timepoint
library(corrplot)

for(t in 1:length(uniq_timept)){
  idx <- (mouse$Timepoint == uniq_timept[t])
  correlation <- cor(mouse[idx, 6:59]) 
  corrplot(correlation, method = "circle")
}

for(i in 1:length(uniq_grp3)){
  for(t in 1:length(uniq_timept)){
    idx <-  (mouse$Grouped.genotype == uniq_grp3[i]) & (mouse$Timepoint == uniq_timept[t])
    corrplot(cor(mouse[idx, 6:59]))
  }
}
