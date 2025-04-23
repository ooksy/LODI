# 결과 : posterior mean/last iteration과 true값을 비교 #########################
post_latent
r_it
alpha_ijt
lambda_j
eta_it
rho_est
sigma_t

graphics.off()
par("mar")
# posterior of latent variable check ###########################################
plot(post_latent[,,1,n_sim], log(data[,,1]+0.01))
plot(post_latent[,,2,n_sim], log(data[,,2]+0.01))
plot(post_latent[,,3,n_sim], log(data[,,3]+0.01))


# r_it check ###################################################################
post_mean_r <- matrix(NA, nrow = n, ncol = T)
for(i in 1:n){
  for(t in 1:T){
    post_mean_r[i,t] <- mean(r_it[i,t,])
  }
}
post_mean_r ; r_it[,,n_sim] ; r[,1,]
plot(r[,1,], r_it[,,n_sim] , xlab = "true values", ylab = "last samples", main = "r_it")
abline(0, 1, col = "red")

par(mar = c(1, 1, 1, 1))  


for(t in 1:T){
  
  png(file = paste("C:/Users/SEC/Desktop/research/25-1/week8/r_i",t,".png", sep=""),
      width = 6,
      height = 5,
      units = "in",
      res = 1200)
  
  plot(r[,1,t], r_it[,t,n_sim] , xlab = "true values", ylab = "last samples", main = paste("r_it at time point", t))
  abline(0, 1, col = "red")
  
  dev.off()
  
}


graphics.off()
ts.plot(r_it[1,1,]) # trace plot seems good
ts.plot(r_it[2,1,]) # trace plot seems good



# alpha check ##################################################################
alpha_mean <- numeric()

# calculating posterior mean
for(t in 1:T){
  alpha_mean[t] <- mean(alpha_ijt[1,1,t,])  
}
alpha_mean
cbind("post mean"=alpha_mean, "last value" = alpha_ijt[1,1,,n_sim], "true" = alpha[1,1,])


alpha_ijt[,,,n_sim] ; alpha

# creating image
par(mar = c(1,1,1,1))

for(t in 1:T){
  
  png(file = paste("C:/Users/SEC/Desktop/research/25-1/week8/alpha_ij",t,".png", sep=""),
      width = 6,
      height = 5,
      units = "in",
      res = 1200)
  
  plot(alpha[,,t], alpha_ijt[,,t,n_sim], xlab = "true alpha",  ylab = "last samples", main = paste("normalized abundance at time point", t))
  abline(0, 1, col = "red")
  
  dev.off()
  
}


graphics.off()

# check trace plot
ts.plot(alpha_ijt[1,1,1,])




# alpha0
library(ComplexHeatmap)
alpha0 ; alpha0_update[,,n_sim]

rownames(alpha0) <- paste0("obs",1:20)
colnames(alpha0) <- paste0("otu",1:30)
rownames(alpha0_update[,,n_sim]) <- paste0("obs",1:20)
colnames(alpha0_update[,,n_sim]) <- paste0("otu",1:30)

ht1 <- Heatmap(alpha0, column_order = colnames(alpha0), row_order = rownames(alpha0), row_title = "observation 1:20"
               , column_title = "true alpha0")
ht2 <- Heatmap(alpha0_update[,,n_sim], column_order = colnames(alpha0_update[,,n_sim]), row_order = rownames(alpha0_update[,,n_sim]),
               , column_title = "estimated alpha0")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "comparing true and estimated alpha0")


# latent variable (group member)
alpha0_grp ; alpha0_grp_update[,,n_sim]

rownames(alpha0_grp) <- paste0("row",1:20)
colnames(alpha0_grp) <- paste0("column",1:30)
rownames(alpha0_grp_update[,,n_sim]) <- paste0("row",1:20)
colnames(alpha0_grp_update[,,n_sim]) <- paste0("column",1:30)

ht1 <- Heatmap(alpha0_grp, column_order = colnames(alpha0_grp), row_order = rownames(alpha0_grp),
               , column_title = "alpha0_group")
ht2 <- Heatmap(alpha0_grp_update[,,n_sim], column_order = colnames(alpha0_grp_update[,,n_sim]), row_order = rownames(alpha0_grp_update[,,n_sim]),
               , column_title = "estimated alpha0_group")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "comparing true and estimated alpha0_grp")


# weight
round(w, digits = 2) ; round(w_update[n_sim,], digits = 2)

# 별로 좋은 visualization이 아님
# par(mfrow=c(1,2))
# pie(round(w[1:7], digits = 2), labels = round(w[1:7], digits = 2), main= "true weight")
# pie(round(w_update[n_sim,1:7], digits = 2), labels = round(w_update[n_sim,1:7], digits = 2), main= "estimated weight")


# mu_j0
mu0_j ; mu_j0_update[n_sim,]


# mu_jl
mu_jl ; mu_jl_update[,,n_sim]

rownames(mu_jl) <- paste0("cluster",1:30)
colnames(mu_jl) <- paste0("otu",1:30)
rownames(mu_jl_update[,,n_sim]) <- paste0("cluster",1:30)
colnames(mu_jl_update[,,n_sim]) <- paste0("otu",1:30)

ht1 <- Heatmap(mu_jl, column_order = colnames(mu_jl), row_order = rownames(mu_jl), row_title = "cluster 1:30",
               , column_title = "true mu_jl")
ht2 <- Heatmap(mu_jl_update[,,n_sim], column_order = colnames(mu_jl_update[,,n_sim]), row_order = rownames(mu_jl_update[,,n_sim]),
               , column_title = "estimated mu_jl")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "comparing true and estimated mu_jl")



# alpha over time ##############################################################

par(mfrow = c(1,4))

# make array to matrix
alpha0
alpha1 <- alpha[,,1]
alpha2 <- alpha[,,2]
alpha3 <- alpha[,,3]
table(alpha0_grp)


# whole samples
par(mfrow = c(1,4))
hist(alpha0, prob = T); lines(density(alpha0_update[,,n_sim]), col = "red")
hist(alpha1, prob = T); lines(density(alpha_ijt[,,1,n_sim]), col = "red")
hist(alpha2, prob = T); lines(density(alpha_ijt[,,2,n_sim]), col = "red")
hist(alpha3, prob = T); lines(density(alpha_ijt[,,3,n_sim]), col = "red")
# alpha0 group 제대로 추정을 못하는 모습 X
# alpha group 시간이 갈수록 좀 불명확해지는 모습

# creating image
# par(mar = c(1,1,1,1))
par(mar = c(5, 4, 4, 2))
otu_sam <- sample(1:30, 5)


for(j in otu_sam){
  # png(file = paste("C:/Users/SEC/Desktop/research/25-1/week10/density",j,".png", sep=""),
  #     width = 10,
  #     height = 8,
  #     units = "in",
  #     res = 1200)
  
  par(mfrow = c(1,4))
  max_alpha <- max(alpha0[,j], alpha0_update[,j,n_sim])
  min_alpha <- min(alpha0[,j], alpha0_update[,j,n_sim])
  hist(alpha0[,j], breaks = 20, probability = T, main = paste("OTU",j,", time 0"), xlim = c(min_alpha, max_alpha))
  lines(density(alpha0_update[,j,n_sim]), col = "red")
  
  max_alpha <- max(alpha1[,j], alpha_ijt[,j,1,n_sim])
  min_alpha <- min(alpha1[,j], alpha_ijt[,j,1,n_sim])
  hist(alpha1[,j], breaks = 20, probability = T, main = paste("OTU",j,", time 1"), xlim = c(min_alpha, max_alpha))
  lines(density(alpha_ijt[,j,1,n_sim]), col = "red")
  
  max_alpha <- max(alpha2[,j], alpha_ijt[,j,2,n_sim])
  min_alpha <- min(alpha2[,j], alpha_ijt[,j,2,n_sim])
  hist(alpha2[,j], breaks = 20, probability = T, main = paste("OTU",j,", time 2"), xlim = c(min_alpha, max_alpha))
  lines(density(alpha_ijt[,j,2,n_sim]), col = "red")
  
  max_alpha <- max(alpha3[,j], alpha_ijt[,j,3,n_sim])
  min_alpha <- min(alpha3[,j], alpha_ijt[,j,3,n_sim])
  hist(alpha3[,j], breaks = 20, probability = T, main = paste("OTU",j,", time 3"), xlim = c(min_alpha, max_alpha))
  lines(density(alpha_ijt[,j,3,n_sim]), col = "red")
  
  # dev.off()
}

par(mfrow = c(2,2), mar = c(4, 4, 2, 1))

overall_min <- min(alpha0, alpha1, alpha2, alpha3, alpha0_update, alpha_ijt)
overall_max <- max(alpha0, alpha1, alpha2, alpha3, alpha0_update, alpha_ijt)

hist(alpha0, xlim = c(overall_min, overall_max), prob = T); lines(density(alpha0_update[,,n_sim]), col = "red")
hist(alpha1, xlim = c(overall_min, overall_max), prob = T); lines(density(alpha_ijt[,,1,n_sim]), col = "red")
hist(alpha2, xlim = c(overall_min, overall_max), prob = T); lines(density(alpha_ijt[,,2,n_sim]), col = "red")
hist(alpha3, xlim = c(overall_min, overall_max), prob = T); lines(density(alpha_ijt[,,3,n_sim]), col = "red")

for(j in otu_sam){
  hist(alpha0[,j], breaks = 20, probability = TRUE, main = paste("OTU",j,", time 0"),
       xlim = c(overall_min, overall_max), xlab = "Alpha values", col = "gray")
  lines(density(alpha0_update[,j,n_sim]), col = "red")
  
  hist(alpha1[,j], breaks = 20, probability = TRUE, main = paste("OTU",j,", time 1"),
       xlim = c(overall_min, overall_max), xlab = "Alpha values", col = "gray")
  lines(density(alpha_ijt[,j,1,n_sim]), col = "red")
  
  hist(alpha2[,j], breaks = 20, probability = TRUE, main = paste("OTU",j,", time 2"),
       xlim = c(overall_min, overall_max), xlab = "Alpha values", col = "gray")
  lines(density(alpha_ijt[,j,2,n_sim]), col = "red")
  
  hist(alpha3[,j], breaks = 20, probability = TRUE, main = paste("OTU",j,", time 3"),
       xlim = c(overall_min, overall_max), xlab = "Alpha values", col = "gray")
  lines(density(alpha_ijt[,j,3,n_sim]), col = "red")
}

# r + alpha check ##############################################################
r_it_estimated <- array(0, dim = c(n, J, T))
r_it_estimated[,,1] <- r_it[,1,n_sim]
r_it_estimated[,,2] <- r_it[,2,n_sim]
r_it_estimated[,,3] <- r_it[,3,n_sim]

par(mfrow = c(1,1))
plot(r + alpha, r_it_estimated[,,] + alpha_ijt[,,,n_sim], xlab = "true (r + alpha)", ylab = "estimated (r + alpha)", main = "r + alpha")
abline(0, 1, col = "red")


# for(t in 1:T){
#   
#   png(file = paste("C:/Users/SEC/Desktop/research/25-1/week8/mu_ij",t,".png", sep=""),
#       width = 6,
#       height = 5,
#       units = "in",
#       res = 1200)
#   
#   plot(r[,,t] + alpha[,,t], r_it_estimated[,,t] + alpha_ijt[,,t,n_sim], xlab = "true (r + alpha)", ylab = "estimated (r + alpha)", main = paste("r + alpha at time point", t))
#   abline(0, 1, col = "red")
#   
#   dev.off()
#   
# }

graphics.off()




# phi check ####################################################################
phi1[n_sim] ; phi_g1
phi2[n_sim] ; phi_g2
ts.plot(phi1)
ts.plot(phi2)


# lambda check #################################################################
lambda_j[,,n_sim]
t(Lambda)
ts.plot(lambda_j[1,1,])


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
# browseVignettes("ComplexHeatmap")


# Identifiability issue 때문에 곱 자체를 비교해야함.. 그건 일정하기 때문.. why?
# lambda'lambda 비교
true <- Lambda %*% t(Lambda)
estimated <- t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim]

compare_mat <- true
compare_mat[upper.tri(compare_mat)] <- estimated[upper.tri(estimated)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:30)
colnames(compare_mat) <- paste0("column", 1:30)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = "upper : estimated", name = "lambda'lambda") # diagonal part is true value, lower - true, upper - posterior estimate

# lower triangular is truth, upper triangular is posterior

# lambda*eta 비교 ##############################################################
# (i, j) : i번째 otu에의 j번째 observation에의 factor의 영향 (time = 1)
true_f <- Lambda%*%eta[,,1]
estimated_f <- t(lambda_j[,,n_sim])%*%eta_it[,,1,n_sim]

rownames(estimated_f) <- paste0("row",1:30)
colnames(estimated_f) <- paste0("column",1:20)

rownames(true_f) <- paste0("row",1:30)
colnames(true_f) <- paste0("column",1:20)

ht1 <- Heatmap(true_f, column_order = colnames(true_f), row_order = rownames(true_f),
               row_title = "OTU", column_title = "True")
ht2 <- Heatmap(estimated_f, column_order = colnames(estimated_f), row_order = rownames(estimated_f),
               row_title = "OTU", column_title = "Estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda'eta")


# eta check ####################################################################
eta_it[,,1,n_sim]
eta[,,1]

true_f <- eta[,,1]
estimated_f <- eta_it[,,1,n_sim]

rownames(true_f) <- paste0("row",1:3)
colnames(true_f) <- paste0("column",1:20)

rownames(estimated_f) <- paste0("row",1:3)
colnames(estimated_f) <- paste0("column",1:20)

ht1 <- Heatmap(estimated_f, column_order = colnames(estimated_f), row_order = rownames(estimated_f),
               row_title = "num of factor", column_title = "Estimated")
ht2 <- Heatmap(true_f, column_order = colnames(true_f), row_order = rownames(true_f),
               row_title = "num of factor", column_title = "True")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Eta")


graphics.off()
ts.plot(eta_it[1,1,1,])
ts.plot(eta_it[2,1,1,])


# rho check ####################################################################
rho
rho_est[n_sim]


ts.plot(rho_est) 
mean(rho_est)

# sigma_t check ################################################################
cbind("post mean" = mean(sigma_t), "last value" = sigma_t[n_sim], "true value" = sigma_sq)

ts.plot(sigma_t)


# plot(r + alpha, mu_updated, xlab = "true (r + alpha + lambda*eta)", ylab = "estimated (r + alpha + lambda*eta)", main = "r + alpha + lambda*eta")
# abline(0, 1, col = "red") 


# number of 0s in simulated data ###############################################
# number of 0s in simulated data 
num0t1 <- sum(floor(exp(r_it_estimated[,,1] + alpha_ijt[,,1,n_sim])) == 0)
num0t2 <- sum(floor(exp(r_it_estimated[,,2] + alpha_ijt[,,2,n_sim])) == 0)
num0t3 <- sum(floor(exp(r_it_estimated[,,3] + alpha_ijt[,,3,n_sim])) == 0)

cbind("time 1" = num0t1 , "time 2" =  num0t2 , "time 3" = num0t3)

# number of 0s in true data
true0t1 <- sum(data[,,1] == 0)
true0t2 <- sum(data[,,2] == 0)
true0t3 <- sum(data[,,3] == 0)

cbind("time 1" = true0t1 , "time 2" =  true0t2 , "time 3" = true0t3)

# ratio of 0s in simulated data 
num0ratio1 <- sum(floor(exp(r_it_estimated[,,1] + alpha_ijt[,,1,n_sim])) == 0)/length(data[,,1])
num0ratio2 <- sum(floor(exp(r_it_estimated[,,2] + alpha_ijt[,,2,n_sim])) == 0)/length(data[,,2])
num0ratio3 <- sum(floor(exp(r_it_estimated[,,3] + alpha_ijt[,,3,n_sim])) == 0)/length(data[,,3])

cbind("time 1" = num0ratio1 , "time 2" =  num0ratio2 , "time 3" = num0ratio3)

# ratio of 0s in true data 
truenum0ratio1 <- sum(data[,,1] == 0)/length(data[,,1])
truenum0ratio2 <- sum(data[,,2] == 0)/length(data[,,2])
truenum0ratio3 <- sum(data[,,3] == 0)/length(data[,,3])

cbind("time 1" = truenum0ratio1 , "time 2" =  truenum0ratio2 , "time 3" = truenum0ratio3)


# covariance between time point check ##########################################
# Cov(Y_it) 
Cov_Y_tr <- (1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J)
Cov_Y_est <- (1/(1-mean(rho_est)^2))*t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim] + sigma_t[n_sim]*diag(J)
sample_cov1 <- cov(y_tilde_t[,,1])
sample_cov2 <- cov(y_tilde_t[,,2])
sample_cov3 <- cov(y_tilde_t[,,3])


rownames(Cov_Y_tr) <- paste0("row",1:30)
colnames(Cov_Y_tr) <- paste0("column",1:30)

rownames(Cov_Y_est) <- paste0("row",1:30)
colnames(Cov_Y_est) <- paste0("column",1:30)

compare_mat <- Cov_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:30)
colnames(compare_mat) <- paste0("column", 1:30)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


# compare with sample covariance
rownames(sample_cov1) <- paste0("row",1:30)
colnames(sample_cov1) <- paste0("column",1:30)

rownames(Cov_Y_est) <- paste0("row",1:30)
colnames(Cov_Y_est) <- paste0("column",1:30)

compare_mat <- sample_cov1
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:30)
colnames(compare_mat) <- paste0("column", 1:30)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : sample cov", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


# time point1 and time point2
Cov_Y_tr <- (rho/(1-rho^2))*Lambda%*%t(Lambda)
Cov_Y_est <- (rho_est[n_sim]/(1-mean(rho_est[n_sim])^2))*t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim]
# rho_est[n_sim]/mean(rho_est) use

rownames(Cov_Y_tr) <- paste0("row",1:30)
colnames(Cov_Y_tr) <- paste0("column",1:30)

rownames(Cov_Y_est) <- paste0("row",1:30)
colnames(Cov_Y_est) <- paste0("column",1:30)

compare_mat <- Cov_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:30)
colnames(compare_mat) <- paste0("column", 1:30)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cov(Y_it, Y_i(t+1)) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


# 0.7 / 0.82
1/(1-0.7^2)
1/(1-0.82^2)
(1/(1-0.82^2))/(1/(1-0.7^2))

Cov_Y_tr <- (1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J)
Cov_Y_est <- (1/(1-rho_est[n_sim]^2))*t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim] + sigma_t[n_sim]*diag(J)

Cov_Y_est[1,1]/Cov_Y_tr[1,1]
Cov_Y_est[20,20]/Cov_Y_tr[20,20]


times <- matrix(0, nrow = 30, ncol = 30)
for(i in 1:30){
  for(j in i:30){
    times[i,j] <- print(Cov_Y_est[i,j]/Cov_Y_tr[i,j])
  }
}

times
summary(times[upper.tri(times)])