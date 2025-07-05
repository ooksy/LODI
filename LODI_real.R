rm(list = ls())

library(mvtnorm)
library(MCMCpack)
library(truncnorm) 
library(circlize)
library(ComplexHeatmap)

set.seed(123) 
N <- 20     # Number of subjects
T <- 7      # Time points
J <- 54     # Observed variables
K <- 5      # Latent factors
L <- 30     # Dirichlet mixture group


Y_count_true <- Data_true


# ----------------------
# MCMC Setup
# ----------------------
# Priors
r_prior_mean <- 0
r_prior_var <- 1^2
c_prior <- 1
alpha_prior_mean <- 0
alpha_prior_var <- 5^2
alpha_a0 <- 2; alpha_b0 <- 1
V_prior_shape <- 2
V_prior_rate <- 1
phi_prior_mean <- 0
phi_prior_var <- 1
a0 <- 2; b0 <- 1 # Inverse-Gamma for sigma2
rho_prior_mean <- 0  # Prior mean for rho (truncated normal)
rho_prior_var <- 1   # Prior variance for rho (truncated normal)
Q_prior_shape <- 2     # Inverse-Gamma shape for Q
Q_prior_rate <- 1      # Inverse-Gamma rate for Q

# Initialize parameters
Y <- array(0, c(N, T, J))
r <- array(0, c(N, T, J))
alpha <- array(0, c(N, T, J))
alpha_g <- matrix(sample(1:3, N*J, replace = T), N, J)
weight <- rdirichlet(1, rep(1, L))
v <- rep(0.5, (L-1))
alpha_mu <- rep(0, L)
alpha_sig2 <- rep(1, L)
phi <- 0.5
V <- 1
Lambda <- matrix(0, J, K)
sigma2 <- 1
eta <- array(rnorm(N*T*K), c(N, T, K))
rho <- 0.5    # AR(1) coefficient
Q <- 1      # AR(1) innovation variance

# MCMC settings
n_iter <- 20000
burn_in <- n_iter/2
thin <- 2
keep <- seq(burn_in + 1, n_iter, by = thin)


# Storage
Y_samples <- array(NA, dim = c(length(keep), N, T, J))
r_samples <- array(NA, dim = c(length(keep), N, T, J))
alpha_samples <- array(NA, dim = c(length(keep), N, T, J))
alpha_g_samples <- array(NA, dim = c(length(keep), N, J))
weight_samples <- array(NA, dim = c(length(keep), L))
alpha_mu_samples <- array(NA, dim = c(length(keep), L))
alpha_sig2_samples <- array(NA, dim = c(length(keep), L))
phi_samples <- numeric(length(keep))
V_samples <- numeric(length(keep))
Lambda_samples <- array(NA, dim = c(length(keep), J, K))
sigma2_samples <- numeric(length(keep))
eta_samples <- array(NA, dim = c(length(keep), N, T, K))
rho_samples <- numeric(length(keep))
Q_samples <- numeric(length(keep))

# ----------------------
# FFBS Implementation
# ----------------------
ffbs_alpha <- function(alpha_g, alpha_mu, alpha_sig2, phi, V, Y, r, Lambda, eta){
  Time <- length(Y)
  
  # Kalman filter
  m <- rep(0, Time)
  C <- rep(0, Time)
  a <- rep(0, Time)
  R <- rep(0, Time)
  
  # Initialize
  m[1] <- alpha_mu[alpha_g]
  C[1] <- alpha_sig2[alpha_g]
  
  # Forward Filter
  for(t in 1:Time){
    if(t > 1){
      a[t] <- phi * m[t-1]
      R[t] <- (phi^2) * C[t-1] + V
    } else {
      a[t] <- m[t]
      R[t] <- C[t]
    }
    
    k_gain <- R[t]*solve(R[t]+sigma2)
    m[t] <- a[t] + k_gain*(Y[t] - r[t] - a[t] - Lambda %*% eta[t,])
    C[t] <- R[t] - k_gain*R[t]
  }
  
  alpha_st <- rep(0, Time)
  alpha_st[Time] <- rnorm(1, mean = m[Time], sd = sqrt(C[Time]))
  
  # Backwards sampling
  for(t in (Time-1):1){
    m[t] <- solve(phi^2/V + 1/C[t]) * (phi*alpha_st[t+1]/V + m[t]/C[t])
    C[t] <- solve(phi^2/V + 1/C[t])
    
    alpha_st[t] <- rnorm(1, mean = m[t], sd = sqrt(C[t])) 
  }
  
  return(alpha_st)
}

ffbs <- function(y, r, alpha, Lambda, rho, Q, sigma2) {
  Time <- dim(y)[1]
  K <- ncol(Lambda)
  
  # Kalman filter
  mu_pred <- matrix(0, Time, K)
  Sigma_pred <- array(0, dim = c(Time, K, K))
  mu_filt <- matrix(0, Time, K)
  Sigma_filt <- array(0, dim = c(Time, K, K))
  
  # Initialize
  mu_filt[1,] <- rep(0, K)
  Sigma_filt[1,,] <- diag(K)
  
  for(t in 1:Time) {
    if(t > 1) {
      mu_pred[t,] <- rho * mu_filt[t-1,]
      Sigma_pred[t,,] <- (rho^2) * Sigma_filt[t-1,,] + diag(Q, K)
    } else {
      mu_pred[t,] <- mu_filt[1,]
      Sigma_pred[t,,] <- Sigma_filt[1,,]
    }
    
    S <- Lambda %*% Sigma_pred[t,,] %*% t(Lambda) + diag(sigma2, J)
    K_gain <- Sigma_pred[t,,] %*% t(Lambda) %*% solve(S)
    resid <- y[t,] - r[t,] - alpha[t,] - Lambda %*% mu_pred[t,]
    
    mu_filt[t,] <- mu_pred[t,] + K_gain %*% resid
    Sigma_filt[t,,] <- (diag(K) - K_gain %*% Lambda) %*% Sigma_pred[t,,]
  }
  
  # Backward sampling
  eta <- matrix(0, Time, K)
  eta[Time,] <- rmvnorm(1, mu_filt[Time,], Sigma_filt[Time,,])
  
  for(t in (Time-1):1) {
    # Corrected matrix operations
    J_smooth <- Sigma_filt[t,,] %*% (rho * diag(K)) %*% 
      solve(rho^2 * Sigma_filt[t,,] + diag(Q, K))
    mu_smooth <- mu_filt[t,] + J_smooth %*% (eta[t+1,] - rho * mu_filt[t,])
    Sigma_smooth <- Sigma_filt[t,,] - J_smooth %*% (rho * Sigma_filt[t,,])
    eta[t,] <- rmvnorm(1, mu_smooth, Sigma_smooth)
  }
  
  return(eta)
}

# ----------------------
# MCMC Loop
# ----------------------
for(iter in 1:n_iter) {
  
  # Y <- Y_true
  # r <- r_true
  # alpha <- alpha_true
  # alpha_g <- alpha_g_true
  # weight <- weight_true
  # alpha_mu <- alpha_mu_true
  # alpha_sig2 <- alpha_sig2_true
  # phi <- phi_true
  # V <- V_true
  # Lambda <- Lambda_true
  # eta <- eta_true
  # rho <- rho_true
  # Q <- Q_true
  # sigma2 <- sigma2_true
  
  
  # Sampling latent variable Y
  for(i in 1:N){
    for(j in 1:J){
      for(t in 1:T){
        if(Y_count_true[i,t,j]==0){
          Y[i,t,j] <- rtruncnorm(1, a = -Inf, b = 0,
                                 mean = r[i,t,j] + alpha[i,t,j] + Lambda[j,] %*% eta[i,t,], sd = sqrt(sigma2))
        } else {
          Y[i,t,j] <- rtruncnorm(1, a=log(Y_count_true[i,t,j]), b=log(Y_count_true[i,t,j]+1),
                                 mean = r[i,t,j] + alpha[i,t,j] + Lambda[j,] %*% eta[i,t,], sd = sqrt(sigma2))
        }
      }
    }
  }
  
  # Update r
  for(i in 1:N){
    for(t in 1:T){
      r_res <- Y[i,t,] - alpha[i,t,] - Lambda %*% eta[i,t,]
      r_res_sum <- sum(r_res)
      r_mean <- solve(1/r_prior_var + J/sigma2)*(r_prior_mean/r_prior_var + r_res_sum/sigma2)
      r_var <- solve(1/r_prior_var + J/sigma2)
      r[i,t,] <- rtruncnorm(1, a=0, b = Inf, mean = r_mean, sd = sqrt(r_var))
    }
  }
  
  # # Update alpha
  # Update alpha group
  for(i in 1:N){
    for(j in 1:J){
      group_prob <- weight * dnorm(alpha[i,1,j], mean = alpha_mu, sd = sqrt(alpha_sig2))
      group_prob <- group_prob / sum(group_prob)
      
      if(all(group_prob == rep(0, L))) {group_prob <- rep(1/L, L)}
      
      alpha_g[i,j] <- sample(1:L, 1, prob = group_prob)
    }
  }
  
  # Update alpha - FFBS
  for(i in 1:N){
    for(j in 1:J){
      alpha[i,,j] <- ffbs_alpha(alpha_g[i,j], alpha_mu, alpha_sig2, phi, V, Y[i,,j], r[i,,j], Lambda[j,], eta[i,,])
    }
  }
  
  
  # Update alpha weight
  for(l in 1:(L-1)){
    v[l] <- rbeta(1, shape1 = 1 + sum(alpha_g == l), shape2 = c_prior + sum(alpha_g > l))
  }
  
  weight[1] <- v[1]
  for(i in 2:L){
    weight[l] <- v[l]*prod(1-v[1:(l-1)])
  }
  if(sum(weight[1:(L-1)]) >= 1) weight[L] <- 0 else weight[L] <- 1-sum(weight[1:(L-1)])
  
  # Update alpha_mu
  for(l in 1:L){
    alpha_mu_var <- solve(1/alpha_prior_var + sum(alpha_g == l)/alpha_sig2[l])
    alpha_mu_mean <- alpha_mu_var * (alpha_prior_mean/alpha_prior_var + sum(alpha[,1,][alpha_g == l])/alpha_sig2[l])
    alpha_mu[l] <- rnorm(1, mean =alpha_mu_mean , sd = sqrt(alpha_mu_var))
  }
  
  # Update alpha_sig2
  for(l in 1:L){
    alpha_sig2[l] <- 1/rgamma(1, alpha_a0 + sum(alpha_g == l)/2,
                              alpha_b0 + sum((alpha[,1,][alpha_g == l] - alpha_mu[l])^2)/2)
  }
  
  # Update phi
  alpha_prev <- alpha[,1:(T-1),]
  alpha_curr <- alpha[,2:T,]
  SS_prev <- sum(alpha_prev^2)
  SS_pcur <- sum(alpha_prev * alpha_curr)
  phi_var <- 1 / (1/phi_prior_var + SS_prev/V)
  phi_mean <- phi_var * (phi_prior_mean/phi_prior_var + SS_pcur/V)
  phi <- rtruncnorm(1, a=-1, b=1, mean = phi_mean, sd = sqrt(phi_var))
  
  # Update V
  V_res <- alpha[,2:T,] - phi * alpha[,1:(T-1),]
  V_shape <- V_prior_shape + N*(T-1)*J/2
  V_rate <- V_prior_rate + sum(V_res^2)/2
  V <- 1/rgamma(1, shape=V_shape, rate=V_rate)
  
  
  # Update Lambda
  for(j in 1:J) {
    k_indices <- 1:min(j, K)
    E_j <- matrix(eta[,,k_indices], nrow = N*T)
    y_j <- c(Y[,,j] - r[,,j] - alpha[,,j])
    
    V_j_inv <- crossprod(E_j)/sigma2 + diag(1, length(k_indices))
    V_j <- solve(V_j_inv)
    m_j <- V_j %*% crossprod(E_j, y_j)/sigma2
    
    for(k_idx in seq_along(k_indices)) {
      k <- k_indices[k_idx]
      if(k < j) {
        Lambda[j,k] <- rnorm(1, m_j[k_idx], sqrt(V_j[k_idx,k_idx]))
      } else if(k == j) {
        Lambda[j,k] <- rtruncnorm(1, a=0, mean=m_j[k_idx], sd=sqrt(V_j[k_idx,k_idx]))
      }
    }
  }
  
  # Update eta (FFBS)
  for(i in 1:N) {
    eta[i,,] <- ffbs(Y[i,,], r[i,,], alpha[i,,], Lambda, rho, Q, sigma2)
  }
  
  # Update rho
  eta_prev <- eta[,1:(T-1),]
  eta_curr <- eta[,2:T,]
  SS_xx <- sum(eta_prev^2)
  SS_xy <- sum(eta_prev * eta_curr)
  rho_mean <- (SS_xy/Q + rho_prior_mean/rho_prior_var) / (SS_xx/Q + 1/rho_prior_var)
  rho_sd <- sqrt(1 / (SS_xx/Q + 1/rho_prior_var))
  rho <- rtruncnorm(1, a=-1, b=1, mean=rho_mean, sd=rho_sd)
  
  # Update Q
  resid <- eta[,2:T,] - rho * eta[,1:(T-1),]
  Q_shape <- Q_prior_shape + N*(T-1)*K/2
  Q_rate <- Q_prior_rate + sum(resid^2)/2
  Q <- 1/rgamma(1, shape=Q_shape, rate=Q_rate)
  
  # Update sigma2
  residuals <- Y - r - alpha - aperm(apply(eta, 1:2, function(x) Lambda %*% x), c(2,3,1))
  ssr <- sum(residuals^2)
  sigma2 <- 1/rgamma(1, a0 + N*T*J/2, b0 + ssr/2)
  
  # Store samples
  if(iter %in% keep) {
    idx <- which(keep == iter)
    Y_samples[idx,,,] <- Y
    r_samples[idx,,,] <- r
    alpha_samples[idx,,,] <- alpha
    alpha_g_samples[idx,,] <- alpha_g
    weight_samples[idx,] <- weight
    alpha_mu_samples[idx,] <- alpha_mu
    alpha_sig2_samples[idx,] <- alpha_sig2
    phi_samples[idx] <- phi
    V_samples[idx] <- V
    Lambda_samples[idx,,] <- Lambda
    sigma2_samples[idx] <- sigma2
    eta_samples[idx,,,] <- eta
    rho_samples[idx] <- rho
    Q_samples[idx] <- Q
  }
  
  if(iter %% 100 == 0) cat("Iteration:", iter, "\n")
}

save.image("C:/Users/SEC/Desktop/research/25summer/may19th/RealData.RData")
load("C:/Users/SEC/Desktop/research/25summer/may19th/RealData.RData")

# ----------------------
# Posterior Analysis
# ----------------------
# Parameter estimates
Y_mean <- apply(Y_samples, 2:4, mean)
r_mean <- apply(r_samples, 2:4, mean)
alpha_g_final <- alpha_g_samples[length(alpha_g_samples[,1,1]),,]
alpha_mean <- apply(alpha_samples, 2:4, mean)
weight_final <- weight_samples[length(weight_samples[,1]),]
alpha_mu_mean <- apply(alpha_mu_samples, 2, mean)
alpha_sig2_mean <- apply(alpha_sig2_samples, 2, mean)
phi_mean <- mean(phi_samples)
V_mean <- mean(V_samples)
Lambda_mean <- apply(Lambda_samples, 2:3, mean)
eta_mean <- apply(eta_samples, 2:4, mean)
rho_mean <- mean(rho_samples)
Q_mean <- mean(Q_samples)
sigma2_mean <- mean(sigma2_samples)

# 1. Parameter comparison

# 2. Trace plots
gettwd()
setwd("C:/Users/SEC/Desktop/research/25summer/may26th")
png(paste0("tsplot1.png"), width=5, height = 5, units="in",res=500)
par(mfrow=c(2,1))
plot(phi_samples, type="l", main="Trace plot: phi")
plot(V_samples, type="l", main="Trace plot: V")
dev.off()

png(paste0("tsplot2.png"), width=5, height = 5, units="in",res=500)
par(mfrow=c(2,1))
plot(rho_samples, type="l", main="Trace plot: rho")
plot(Q_samples, type="l", main="Trace plot: Q")
dev.off()

# Latent Y Plot
par(mfrow = c(1,T))
for(t in 1:T){
  plot(Y_mean[,t,], log(Y_count_true[,t,]+0.01),
       main = paste("time point", t), xlab = "post mean Y", ylab = "log(data+0.01)")
  abline(0, 1, col = "red")
}

# normalizing constant r
for(t in 1:T){
  cat(paste("normalizing constant at time point", t), round(r_mean[,t,1], digits = 2), "\n")
}

# normalized abundance alpha
# alpha group
table(alpha_g_final)

# alpha
for(t in 1:T){
  png(paste0("alpha",t,".png"), width=8, height = 5, units="in",res=500)
  alp_min <- min(alpha_mean[,t,])
  alp_max <- max(alpha_mean[,t,])
  
  col_fun <- colorRamp2(c(alp_min, median(c(alp_min, alp_max)), alp_max), c("blue", "white", "red"))
  estimated <- alpha_mean[,t,]
  
  rownames(estimated) <- paste("row", 1:N)
  colnames(estimated) <- paste("col", 1:J)
  
  ht2 <- Heatmap(estimated, 
                 column_order = colnames(estimated), 
                 row_order = rownames(estimated),
                 row_title = "Subject 1-20", 
                 column_title = "OTU 1-54", 
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 col = col_fun,
                 name = "alpha value")
  ht_list <- ht2
  draw(ht_list, column_title = paste("posterior mean of alpha at time point",t))
  dev.off()
}

# weight
round(weight_final, digits = 2)

# alpha_mu
alpha_mu_mean 

# alpha_sig2
alpha_sig2_mean

# phi
phi_mean

# V
V_mean

# 3. Heatmap
# Lambda
png(paste0("Lambda.png"), width=5, height = 6, units="in",res=500)
col_fun <- colorRamp2(c(min(Lambda_mean), median(Lambda_mean), max(Lambda_mean)), c("blue", "white", "red"))
estimated <- Lambda_mean
rownames(estimated) <- paste0("row", 1:J)
colnames(estimated) <- paste0("column", 1:K)
ht <- Heatmap(estimated, 
               column_order = colnames(estimated), 
               row_order = rownames(estimated),
               row_title = "OTU 1-54", 
               column_title = "Factor 1-5",
               show_row_names = FALSE,
               show_column_names = FALSE,
               col = col_fun,
               name = "Lambda")
draw(ht, column_title = "Posterior mean of Lambda")
dev.off()

# Eta
for(t in 1:T){
  png(paste0("Eta",t,".png"), width=5, height = 6, units="in",res=500)
  col_fun <- colorRamp2(c(min(eta_mean), median(eta_mean), max(eta_mean)), c("blue", "white", "red"))
  estimated <- eta_mean[,t,]
  rownames(estimated) <- paste0("row", 1:N)
  colnames(estimated) <- paste0("column", 1:K)
  ht <- Heatmap(estimated, 
                 column_order = colnames(estimated), 
                 row_order = rownames(estimated),
                 row_title = "Subject 1-20", 
                 column_title = "Factor 1-5", 
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 col = col_fun,
                 name = "Eta")
  draw(ht, column_title = paste("Posterior mean of eta at time point",t))
  dev.off()
}

# rho
rho_mean

# Q
Q_mean

# Sigma2
sigma2_mean


# LL^T - median
png(paste0("LL^T.png"), width=9, height = 8, units="in",res=500)
tcross_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  tcross_array[ , , i] <- tcrossprod(Lambda_samples[i, , ])
}
# Compute the posterior median across iterations
Lambda_median <- apply(tcross_array, c(1, 2), median)
estimated <- Lambda_median
rownames(estimated) <- paste0("row", 1:J)
colnames(estimated) <- paste0("column", 1:J)
col_fun = colorRamp2(c(-5,0,10), c("blue", "white", "red"))
ht <- Heatmap(estimated, 
              column_order = colnames(estimated), 
              row_order = rownames(estimated),
               row_title = "OTU 1-54", 
              column_title = "OTU 1-54", 
              show_row_names = FALSE,
              show_column_names = FALSE,
              name = "Post median",
              col = col_fun)
draw(ht, column_title = "Posterior median of LL^T")
dev.off()

# Cov - median
png(paste0("Cov.png"), width=9, height = 8, units="in",res=500)
Cov_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  Cov_array[ , , i] <- 1/(1-rho_samples[i]^2) *tcrossprod(Lambda_samples[i, , ]) + Q_samples[i] * diag(J)
}
# Compute the posterior median across iterations
Cov_median <- apply(Cov_array, c(1, 2), median)
estimated <- Cov_median
rownames(estimated) <- paste0("row", 1:J)
colnames(estimated) <- paste0("column", 1:J)
col_fun = colorRamp2(c(-20,0,100), c("blue", "white", "red"))
ht <- Heatmap(estimated, 
               column_order = colnames(estimated), 
               row_order = rownames(estimated),
               row_title = "OTU 1-54", 
               column_title = "OTU 1-54",
              show_row_names = FALSE,
              show_column_names = FALSE,
               name = "Post median")
draw(ht, column_title = "Posterior median of Cov")
dev.off()


rownames(true) <- paste0("row",1:J)
colnames(true) <- paste0("column",1:J)
rownames(Cov_median) <- paste0("row",1:J)
colnames(Cov_median) <- paste0("column",1:J)

compare_mat <- true
compare_mat[upper.tri(compare_mat)] <- Cov_median[upper.tri(Cov_median)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)

Heatmap(compare_mat, 
        column_order = colnames(compare_mat), 
        show_row_names = FALSE,
        show_column_names = FALSE,
        row_order = rownames(compare_mat),
        row_title = "lower : true", 
        column_title = gt_render("Cov(Y_it) <br> upper : estimated"), n
        ame = "mat") # diagonal part is true value, lower - true, upper - posterior estimate



# LL^T - mean
tcross_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  tcross_array[ , , i] <- tcrossprod(Lambda_samples[i, , ])
}
# Compute the posterior median across iterations
Lambda_mean <- apply(tcross_array, c(1, 2), mean)
true <- tcrossprod(Lambda_true)
estimated <- Lambda_mean
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


# Cov - mean
Cov_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  Cov_array[ , , i] <- 1/(1-rho_samples[i]^2) *tcrossprod(Lambda_samples[i, , ]) + Q_samples[i] * diag(J)
}
# Compute the posterior median across iterations
Cov_mean <- apply(Cov_array, c(1, 2), mean)
true <- 1/(1-rho_true^2) *tcrossprod(Lambda_true)+ Q_true * diag(J)
estimated <- Cov_mean
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

rownames(true) <- paste0("row",1:J)
colnames(true) <- paste0("column",1:J)
rownames(Cov_mean) <- paste0("row",1:J)
colnames(Cov_mean) <- paste0("column",1:J)

compare_mat <- true
compare_mat[upper.tri(compare_mat)] <- Cov_mean[upper.tri(Cov_mean)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate
