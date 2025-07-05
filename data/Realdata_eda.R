library(compositions)
library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggalt)

# data
mouse <- read.csv("data/microbiome_complete.csv")
mouse$Genotype.f <- factor(mouse$Genotype, levels = c("WT", "ACKG", "ACK", "HACKG", "HACK"), ordered = T)
mouse$Grouped.genotype.f <- factor(mouse$Grouped.genotype, levels = c("WT", "ACK/G", "HACK/G"), ordered = T)
mouse$Timepoint.f <- factor(mouse$Timepoint, levels = paste0("W", c(0, 1, 3, 7, 11, 15, 19)), ordered = T)

uniq_grp3 <- unique(mouse$Grouped.genotype)
uniq_grp5 <- unique(mouse$Genotype)
uniq_timept <- unique(mouse$Timepoint)

whole_pred <- paste0("OTU",1:54)

# number of zeros
# by group and time
# WT
tbl1 <- mouse %>% 
          group_by(Grouped.genotype.f, Timepoint.f) %>% 
          dplyr::select(all_of(whole_pred)) %>% 
          summarize(across(whole_pred, function(x) sum(x==0)))

# 1) 긴 형식으로 변환 (OTU 컬럼 모두 pivot_longer)
df_long <- tbl1 %>%
  pivot_longer(cols = starts_with("OTU"),
               names_to = "OTU",
               values_to = "value") %>% 
  filter(Grouped.genotype.f == "WT")

df_long$Timepoint.f <- factor(df_long$Timepoint.f, levels = c("W0", "W1", "W3", "W7", "W11", "W15", "W19"))
df_long$OTU <- factor(df_long$OTU, levels = unique(df_long$OTU))

# 2) 마지막 타임포인트 추출
df_labels <- df_long %>%
  group_by(OTU) %>%
  filter(Timepoint.f == max(Timepoint.f))  # 마지막 시점 선택

# 3) 스파게티 플롯 그리기
ggplot(df_long, aes(x = Timepoint.f, y = value, group = OTU, color = OTU)) +
  geom_line() +
  geom_point() +
  geom_text_repel(data = df_labels,
                  aes(label = OTU),
                  nudge_x = 0.3,        # 라벨을 오른쪽으로 약간 이동
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  show.legend=FALSE) +  # 라벨과 점을 잇는 선 없앰
  labs(title = "Spaghetti plot of group WT",
       x = "Timepoint", 
       y = "Number of 0s", 
       color = "OTU") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ACK/G
# 1) 긴 형식으로 변환 (OTU 컬럼 모두 pivot_longer)
df_long <- tbl1 %>%
  pivot_longer(cols = starts_with("OTU"),
               names_to = "OTU",
               values_to = "value") %>% 
  filter(Grouped.genotype.f == "ACK/G")

df_long$Timepoint.f <- factor(df_long$Timepoint.f, levels = c("W0", "W1", "W3", "W7", "W11", "W15", "W19"))
df_long$OTU <- factor(df_long$OTU, levels = unique(df_long$OTU))

# 2) 마지막 타임포인트 추출
df_labels <- df_long %>%
  group_by(OTU) %>%
  filter(Timepoint.f == max(Timepoint.f))  # 마지막 시점 선택

# 3) 스파게티 플롯 그리기
ggplot(df_long, aes(x = Timepoint.f, y = value, group = OTU, color = OTU)) +
  geom_line() +
  geom_point() +
  geom_text_repel(data = df_labels,
                  aes(label = OTU),
                  nudge_x = 0.3,        # 라벨을 오른쪽으로 약간 이동
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  show.legend=FALSE) +  # 라벨과 점을 잇는 선 없앰
  labs(title = "Spaghetti plot of group ACK/G",
       x = "Timepoint", 
       y = "Number of 0s", 
       color = "OTU") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# HACK/G
# 1) 긴 형식으로 변환 (OTU 컬럼 모두 pivot_longer)
df_long <- tbl1 %>%
  pivot_longer(cols = starts_with("OTU"),
               names_to = "OTU",
               values_to = "value") %>% 
  filter(Grouped.genotype.f == "HACK/G")

df_long$Timepoint.f <- factor(df_long$Timepoint.f, levels = c("W0", "W1", "W3", "W7", "W11", "W15", "W19"))
df_long$OTU <- factor(df_long$OTU, levels = unique(df_long$OTU))

# 2) 마지막 타임포인트 추출
df_labels <- df_long %>%
  group_by(OTU) %>%
  filter(Timepoint.f == max(Timepoint.f))  # 마지막 시점 선택

# 3) 스파게티 플롯 그리기
ggplot(df_long, aes(x = Timepoint.f, y = value, group = OTU, color = OTU)) +
  geom_line() +
  geom_point() +
  geom_text_repel(data = df_labels,
                  aes(label = OTU),
                  nudge_x = 0.3,        # 라벨을 오른쪽으로 약간 이동
                  direction = "y",
                  hjust = 0,
                  segment.color = NA,
                  show.legend=FALSE) +  # 라벨과 점을 잇는 선 없앰
  labs(title = "Spaghetti plot of group HACK/G",
       x = "Timepoint", 
       y = "Number of 0s", 
       color = "OTU") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# # by subject and time
# tbl2 <- mouse %>% 
#           group_by(Mouse_ID, Timepoint.f) %>% 
#           select(all_of(whole_pred)) %>% 
#           summarize(across(whole_pred, function(x) sum(x==0)))

# Timepoint 별 0의 개수
predictors <- mouse[,whole_pred]
apply(predictors, 2, function(x) sum(x==0))
for(i in seq_along(uniq_timept)){
  res_0 <- apply(predictors[mouse$Timepoint==uniq_timept[i],], 2, function(x) sum(x==0))
  cat("Number of 0s at Timepoint", uniq_timept[i], ":",  "\n") 
  print(res_0)
}

# 0의 비율 
sum(predictors==0)/(nrow(predictors)*ncol(predictors))

for(i in seq_along(uniq_timept)){
  mat <- predictors[mouse$Timepoint==uniq_timept[i],]
  res_0 <- sum(mat == 0)/(nrow(mat)*ncol(mat))
  cat("Ratio of 0s at Timepoint", uniq_timept[i], ":", round(res_0, 2), "\n") 
}


# histograms


# Data transformation 후에 분석
# transformation - centered log ratio/modified centered log ratio
mclr  <- function(X, c){
  X <- as.matrix(X)
  p <- ncol(X)
  for(i in 1:nrow(X)) {
    q <- sum(X[i, ] == 0)
    geo_mean <- exp(sum(log(X[i, X[i, ] != 0])))^(1/(p-q))
    eps <- abs(min(log(X[i, X[i, ] != 0]/geo_mean))) + c
    X[i, X[i, ] != 0] <- log(X[i, X[i, ] != 0]/geo_mean) + eps
  }
  return(X)
}

predictors_trns <- mclr(predictors, c = 1)
predictors_trns <- clr(predictors)


# mean 
# variance/covariance