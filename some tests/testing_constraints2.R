rm(list = ls())

# library 
library(mvtnorm)
library(ComplexHeatmap)
library(gridtext)
library(truncnorm)

# set seed
set.seed(124)

# Data generation
# initial setting
n <- 20
J <- 30
K <- 3
T <- 3

# generating parameter
# factor loading
lambda <- matrix(0, nrow = J, ncol = K) # lambda_sample은 K*J임
for(j in 1:J) {
  for(k in 1:min(j, K)) {
    if(j == k) {
      lambda[j, k] <- 1
    } else {
      lambda[j, k] <- rnorm(1, 0, 0.5)
    }
  }
}

sigma2 <- 0.5
rho <- 0.7
Q <- 0.2

# factor
eta <- array(0, c(K, n, T))
eta[,,1] <- t(rmvnorm(n, mean = rep(0, K), sigma = diag(1, K)))
for(t in 2:T){
  eta[,,t] <- rho*eta[,,t-1] + t(rmvnorm(n, mean = rep(0, K), sigma = diag(Q, K)))
}

# data y
y <- array(0, dim=c(n, J, T))

for(t in 1:T){
  y[,,t] <- t(eta[,,t]) %*% t(lambda) + rmvnorm(n, mean = rep(0, J), sigma = diag(sigma2, J))
}

# MCMC setting 
n_sim <- 2000 # 적은 iteration 해보고 50000으로
burn_in <- n_sim/2
thin <- 2
keep <- seq(burn_in + 1, n_sim, by = thin)

# Initialization
lambda_sam <- array(0, dim=c(K, J, n_sim))
lambda_sam[,,1] <- array(rnorm(K*J), dim = c(K, J))

eta_sam <- array(0, dim=c(K, n, T, n_sim))
eta_sam[,,,1] <- array(rnorm(K*n*T), dim = c(K,n,T))

# m_vec <- matrix(0, nrow = T, ncol = K)  # filter
# C_mat <- array(0, dim=c(K,K,T))         # filter
# a_vec <- matrix(0, nrow = T, ncol = K)  # pred
# R_mat <- array(0, dim=c(K,K,T))         # pred

rho_sam <- numeric(n_sim)
rho_sam[1] <- 0.5

sigma2_sam <- numeric(n_sim)
sigma2_sam[1] <- 1

Q_sam <- numeric(n_sim)
Q_sam[1] <- 1

# prior specificaton
sig_a_Q <- 2
sig_b_Q <- 1
sig_a <- 2
sig_b <- 1



# MCMC loop
for(s in 2:n_sim){
  # lambda_sam[,,s-1] <- lambda_sam[,,s] <- t(lambda)
  # eta_sam[,,,s-1] <- eta_sam[,,,s] <- eta
  # rho_sam[s-1] <- rho_sam[s] <- rho
  # Q_sam[s-1] <- Q_sam[s] <- Q
  # sigma2_sam[s-1] <- sigma2_sam[s] <- sigma2
  
  
  # 1. sampling lambda
  for(j in 1:J){

      idx_seq <- 1:min(K, j)
      eta_j <- aperm(eta_sam[,,,s-1], c(2,3,1))
      eta_j <- matrix(eta_j[,,idx_seq], nrow = n*T)
      y_j <- matrix(y[,j,], nrow = n*T)

      lambda_cov <- crossprod(eta_j)/sigma2_sam[s-1] + diag(1, length(idx_seq))
      lambda_cov <- solve(lambda_cov)
      lambda_mean <- lambda_cov %*% crossprod(eta_j, y_j)/sigma2_sam[s-1]


      for(idx_k in seq_along(idx_seq)){
        k <- idx_seq[idx_k]
        if(k < j){
          lambda_sam[k,j,s] <- rnorm(1, mean = lambda_mean[idx_k], sd = sqrt(lambda_cov[idx_k,idx_k]))
        } else if(k == j){
          lambda_sam[k,j,s] <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_mean[idx_k], sd = sqrt(lambda_cov[idx_k,idx_k]))
        }
      }
  }


  # 2. sampling eta (FFBS)
  # forward filtering
  for(i in 1:n){
    
    m_vec <- matrix(0, nrow = T, ncol = K)  # filter
    C_mat <- array(0, dim=c(K,K,T))         # filter
    a_vec <- matrix(0, nrow = T, ncol = K)  # pred
    R_mat <- array(0, dim=c(K,K,T))         # pred
    
    
    # initialization
    m_vec[1,] <- rep(0, K)
    C_mat[,,1] <- diag(K)

    for(t in 1:T){
      if(t > 1){
        a_vec[t,] <- rho_sam[s-1]*m_vec[(t-1),]
        R_mat[,,t] <- (rho_sam[s-1]^2)*C_mat[,,(t-1)] + Q_sam[s-1]*diag(K)
      } else {
        a_vec[t,] <- m_vec[1,]
        R_mat[,,t] <- C_mat[,,1]
      }


      # update filter
      Q_t <- t(lambda_sam[,,s]) %*% R_mat[,,t] %*% lambda_sam[,,s] + sigma2_sam[s-1] * diag(J)
      Q_t_inv <- solve(Q_t)
      K_t <- t(R_mat[,,t]) %*% lambda_sam[,,s] %*% Q_t_inv
      m_vec[t,] <- a_vec[t,] + K_t %*% (y[i,,t] - t(lambda_sam[,,s]) %*% a_vec[t,])

      C_mat[,,t] <- R_mat[,,t] - K_t %*% t(lambda_sam[,,s]) %*% R_mat[,,t]

    }


    eta_sam[,i,T,s] <- rmvnorm(1, mean = m_vec[T,], sigma = C_mat[,,T]) # eta sampling for time point T


  # backward sampling
  for(t in (T-1):1){
      A_t <- rho_sam[s-1] * C_mat[,,t] %*% solve(rho_sam[s-1]^2 * C_mat[,,t] + Q_sam[s-1]*diag(K))
      sm_mean <- m_vec[t,] + A_t %*% (eta_sam[,i,(t+1),s] - rho_sam[s-1]*m_vec[t,])
      sm_cov <- C_mat[,,t] - A_t %*% (rho_sam[s-1]*C_mat[,,t])

      eta_sam[,i,t,s] <- rmvnorm(1, mean = sm_mean, sigma = sm_cov)
    }


  }




  # 3. sampling rho
  # eta time varying equation의 Q도 고려해줘야함

  rho_mean <- 0
  rho_var <- 0

  for(i in 1:n){
    for(t in 2:T){
      rho_var <- rho_var + as.numeric((t(eta_sam[,i,t-1,s]) %*% eta_sam[,i,t-1,s]))
    }
  }

  for(i in 1:n){
    for(t in 2:T){
      rho_mean <- rho_mean + as.numeric((t(eta_sam[,i,t-1,s]) %*% eta_sam[,i,t,s]))
    }
  }

  var_rho <- 1/rho_var
  mean_rho <- rho_mean * var_rho
  var_rho <- Q_sam[s-1]/rho_var


  rho_sam[s] <- rtruncnorm(1, a = 0, b = Inf, mean = mean_rho, sd = sqrt(var_rho))


  # 4. sampling Q
  resid <- eta_sam[,,2:T,s] - rho_sam[s]*eta_sam[,,1:(T-1),s]

  Q_sam[s] <- 1/rgamma(1, sig_a_Q + (n*(T-1)*K/2), sig_b_Q + (sum(resid^2)/2))

  # 5. sampling sigma2
  # method 3
  residual <- y - aperm(apply(eta_sam[,,,s], 2:3, function(x) t(lambda_sam[,,s]) %*% x), c(2,1,3))
  ssr <- sum(residual^2)
  sigma2_sam[s] <- 1/rgamma(1, sig_a + (n*J*T/2), sig_b + ssr/2)
  
  if(s %% 100 == 0) cat("Iteration:", s, "\n")

}

# save.image("C:/Users/SEC/Desktop/research/25summer/may5th/testing_constraints2.RData")

# remove burn-in period & thinning
lambda_sam <- lambda_sam[,,keep]
eta_sam <- eta_sam[,,,keep]
rho_sam <- rho_sam[keep]
Q_sam <- Q_sam[keep]
sigma2_sam <- sigma2_sam[keep]

# calculate the posterior mean
lambda_mean <- apply(lambda_sam, 1:2, mean)
eta_mean <- apply(eta_sam, 1:3, mean)
rho_mean <- mean(rho_sam)
Q_mean <- mean(Q_sam)
sigma2_mean <- mean(sigma2_sam)

# calculate the posterior median
lambda_median <- apply(lambda_sam, 1:2, median)
eta_median <- apply(eta_sam, 1:3, median)
rho_median <- median(rho_sam)
Q_median <- median(Q_sam)
sigma2_median <- median(sigma2_sam)


# compare the result
# 1. lambda
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
true <- lambda
estimated <- t(lambda_mean)
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:K)
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda")


ts.plot(lambda_sam[1,1,])
ts.plot(lambda_sam[2,2,])
ts.plot(lambda_sam[3,3,])
mean(lambda_sam[1,1,])
mean(lambda_sam[2,2,])
mean(lambda_sam[3,3,])


# 2. eta
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
true <- t(eta[,,3])
estimated <- t(eta_mean[,,3])
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:n)
colnames(compare_mat) <- paste0("column", 1:K)
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Eta")

# 3. rho
rho_mean ; rho_median ; rho

# 4. Q
Q_mean ; Q_median ; Q

# 5. sigma2
sigma2_mean ; sigma2_median ; sigma2


# Covariance matrix
cross_array <- array(NA, dim = c(J, J, dim(lambda_sam)[3]))
for (i in 1:(dim(lambda_sam)[3])) {
  cross_array[,,i] <- crossprod(lambda_sam[,,i])
}
# Compute the posterior median across iterations
lambda_median <- apply(cross_array, c(1, 2), median)
true <- tcrossprod(lambda)
estimated <- lambda_median
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "LL^T")



Cov_array <- array(NA, dim = c(J, J, dim(lambda_sam)[3]))
for (i in 1:(dim(lambda_sam)[3])) {
  Cov_array[ , , i] <- 1/(1-rho_sam[i]^2) * crossprod(lambda_sam[,,i]) + Q_sam[i] * diag(J)
}
# Compute the posterior median across iterations
Cov_median <- apply(Cov_array, c(1, 2), median)
true <- 1/(1-rho^2) *tcrossprod(lambda)+ Q * diag(J)
estimated <- Cov_median
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True")
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Cov")

