library(compositions)
library(corrplot)

# data
mouse <- read.csv("data/microbiome_complete.csv")
mouse$Genotype.f <- factor(mouse$Genotype, levels = c("WT", "ACKG", "ACK", "HACKG", "HACK"), ordered = T)
mouse$Grouped.genotype.f <- factor(mouse$Grouped.genotype, levels = c("WT", "ACK/G", "HACK/G"), ordered = T)
mouse$Timepoint.f <- factor(mouse$Timepoint, levels = paste0("W", c(0, 1, 3, 7, 11, 15, 19)), ordered = T)

uniq_grp3 <- unique(mouse$Grouped.genotype)
uniq_grp5 <- unique(mouse$Genotype)
uniq_timept <- unique(mouse$Timepoint)


# Make true data array
Data <- mouse[, 6:59]
Data_true <- array(0, dim=c(20, 7, 54))
nrow(Data)
ncol(Data)

for(i in seq_along(uniq_timept)){
  idx <- uniq_timept[i]
  Data_true[,i,] <- as.matrix(Data[mouse$Timepoint == idx, ])
}

Data_true
log(sum(Data_true[1,1,]))
log(sum(Data_true[2,1,]))
log(sum(Data_true[3,1,]))


sum(log(Data_true + 0.01)[1,1,])
sum(log(Data_true + 0.01)[2,1,])
sum(log(Data_true + 0.01)[3,1,])


Data_true_clr <- array(0, dim=c(20, 7, 54))
for(t in 1:T){
  Data_true_clr[,t,] <- clr(Data_true[,t,])
}


setwd("C:/Users/SEC/Desktop/research/25summer/may19th")
t=7
for(t in 1:T){
png(paste0("cov",t,".png"), width=6, height = 5, units="in",res=300)
cov1 <- cov(Data_true_clr[,t,])
rownames(cov1) <- paste("row", 1:J)
colnames(cov1) <- paste("col", 1:J)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(cov1, 
              column_order = colnames(cov1), 
              row_order = rownames(cov1),
              row_title = "OTU 1-54", 
              column_title = "OTU 1-54",
              show_row_names = FALSE,
              show_column_names = FALSE,
              name = "Covariance")
draw(ht, column_title = paste0("True Cov : time point", t))
dev.off()
}


# density plot along the timepoint : real data
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
