install.packages("matrixcalc")
library(matrixcalc)

A <- matrix(c(1, 2, 3, 4), 2, 2)
B <- matrix(c(2, 3, 4, 5), 2, 2)


frobenius.norm(A)
frobenius.norm(B)
frobenius.norm(A-B)

# normalized frobenius norm?
obs <- 1:n
otus <- 1:J



original_grp_mat <- est_grp_mat <- expand.grid(obs=obs, otus=otus)

original_grp_mat[,3] <- as.vector(alpha0_grp)
est_grp_mat[,3] <- as.vector(alpha0_grp_update[,,n_sim])

original_grp_mat  
est_grp_mat


table(original_grp_mat[,3]) # original group
table(est_grp_mat[,3]) # estimated group, but this is estimated based on original group...
# 1 3 2 5 6 4 / 1 3 2 6 5 4 추정
est_grp_mat[,3]

# change the number of group

# clustering adjusted rand index
install.packages("fossil")
library(fossil)
adj.rand.index(original_grp_mat[,3], est_grp_mat[,3]) # 원 데이터

est_grp_mat[,4] <- est_grp_mat[,3]
est_grp_mat[est_grp_mat[,3] == 3, 4] <- 2 # 3 group을 2 group으로
est_grp_mat[est_grp_mat[,3] == 2, 4] <- 3 # 2 group을 3 group으로
est_grp_mat[est_grp_mat[,3] == 5, 4] <- 4 # 5 group을 4 group으로
est_grp_mat[est_grp_mat[,3] == 6, 4] <- 5 # 6 group을 5 group으로
est_grp_mat[est_grp_mat[,3] == 4, 4] <- 6 # 4 group을 6 group으로
table(est_grp_mat[,4])

adj.rand.index(original_grp_mat[,3], est_grp_mat[,4]) # sample 숫자별로 분류

est_grp_mat[,5] <- est_grp_mat[,3]
est_grp_mat[est_grp_mat[,3] == 3, 5] <- 2 # 3 group을 2 group으로
est_grp_mat[est_grp_mat[,3] == 2, 5] <- 3 # 2 group을 3 group으로
est_grp_mat[est_grp_mat[,3] == 6, 5] <- 4 # 5 group을 4 group으로
est_grp_mat[est_grp_mat[,3] == 5, 5] <- 5 # 6 group을 5 group으로
est_grp_mat[est_grp_mat[,3] == 4, 5] <- 6 # 4 group을 6 group으로
table(est_grp_mat[,5])

adj.rand.index(original_grp_mat[,3], est_grp_mat[,5]) # sample 숫자별로 분류


# normalized mutual information
install.packages("aricode")
library(aricode)

NMI(original_grp_mat[,3], est_grp_mat[,3])
NMI(original_grp_mat[,3], est_grp_mat[,4])
NMI(original_grp_mat[,3], est_grp_mat[,5])
