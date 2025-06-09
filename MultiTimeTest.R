rm(list = ls())

library(mvtnorm)   
library(truncnorm) 
library(circlize)
library(ComplexHeatmap)
set.seed(123) # seed = 124, 125 : there is a similar problem(sign-changing problem)
N <- 20    # Number of subjects
T <- 3      # Time points
J <- 30     # Observed variables
K <- 3      # Latent factors

# True parameters
r_true <- array(0, dim = c(N, T, J))
r_true[,,] <- matrix(runif(N*T, min = 0, max = 2), nrow = N, ncol = T)

# alpha simple case
alpha_g_true <- matrix(0, nrow = N, ncol = J)
alpha_true <- array(0, dim = c(N, T, J))

alpha_g_true[,] <- sample(1:3, size = N*J, replace = T, prob = c(0.5, 0.3, 0.2))
alpha_mean <- c(-4, 2, 8)
alpha_var <- c(0.5, 1, 2)
phi_true <- 0.8
V <- 0.2

for(i in 1:N){
  for(j in 1:J){
    g_idx <- alpha_g_true[i,j]
    alpha_true[i,1,j] <- rnorm(1, alpha_mean[g_idx], alpha_var[g_idx])
  }
}

for(t in 2:T){
  alpha_true[,t,] <- phi_true * alpha_true[,t-1,] + rmvnorm(N, mean = rep(0,J), sigma = diag(V, J))
}



Lambda_true <- matrix(0, nrow = J, ncol = K)
for(j in 1:J) {
  for(k in 1:min(j, K)) {
    if(j == k) {
      Lambda_true[j, k] <- 1
    } else {
      Lambda_true[j, k] <- rnorm(1, 0, 0.5)
    }
  }
}
sigma2_true <- 0.5  # Observation error variance
rho_true <- 0.7     # AR(1) coefficient
Q_true <- 0.2       # Innovation variance


# Simulate latent factors (AR(1) process)
eta_true <- array(0, dim = c(N, T, K))
for(i in 1:N) {
  eta_true[i, 1, ] <- rnorm(K)
  for(t in 2:T) {
    eta_true[i, t, ] <- rho_true * eta_true[i, t-1, ] + rnorm(K, sd = sqrt(Q_true))
  }
}

# Simulate observed data
Y <- array(0, dim = c(N, T, J))
for(i in 1:N) {
  for(t in 1:T) {
    Y[i, t, ] <- r_true[i, t,] + alpha_true[i, t,] + Lambda_true %*% eta_true[i, t, ] + rnorm(J, sd = sqrt(sigma2_true))
  }
}

Y_count <- floor(exp(Y))
Y_count[,1,]
Y_count[,2,]
Y_count[,3,]


# ----------------------
# MCMC Setup
# ----------------------
# Priors
tau2 <- 1          # Prior variance for Lambda
a0 <- 2; b0 <- 1 # Inverse-Gamma for sigma2
rho_prior_mean <- 0  # Prior mean for rho (truncated normal)
rho_prior_var <- 1   # Prior variance for rho (truncated normal)
Q_prior_shape <- 2     # Inverse-Gamma shape for Q
Q_prior_rate <- 1      # Inverse-Gamma rate for Q

# Initialize parameters
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
Lambda_samples <- array(NA, dim = c(length(keep), J, K))
sigma2_samples <- numeric(length(keep))
eta_samples <- array(NA, dim = c(length(keep), N, T, K))
rho_samples <- numeric(length(keep))
Q_samples <- numeric(length(keep))

# ----------------------
# FFBS Implementation
# ----------------------
ffbs <- function(y, Lambda, rho, Q, sigma2) {
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
    resid <- y[t,] - Lambda %*% mu_pred[t,]
    
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
  
  # 1. Update Lambda
  for(j in 1:J) {
    k_indices <- 1:min(j, K)
    E_j <- matrix(eta[,,k_indices], nrow = N*T)
    y_j <- c(Y[,,j])

    V_j_inv <- crossprod(E_j)/sigma2 + diag(1/tau2, length(k_indices))
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
  
  # 2. Update eta (FFBS)
  for(i in 1:N) {
    eta[i,,] <- ffbs(Y[i,,], Lambda, rho, Q, sigma2)
  }
  
  # 3. Update rho
  eta_prev <- eta[,1:(T-1),]
  eta_curr <- eta[,2:T,]
  SS_xx <- sum(eta_prev^2)
  SS_xy <- sum(eta_prev * eta_curr)
  rho_mean <- (SS_xy + rho_prior_mean/rho_prior_var) / (SS_xx + 1/rho_prior_var)
  rho_sd <- sqrt(Q / (SS_xx + 1/rho_prior_var))
  rho <- rtruncnorm(1, a=-1, b=1, mean=rho_mean, sd=rho_sd)
  
  # 4. Update Q
  resid <- eta[,2:T,] - rho * eta[,1:(T-1),]
  Q_shape <- Q_prior_shape + N*(T-1)*K/2
  Q_rate <- Q_prior_rate + sum(resid^2)/2
  Q <- 1/rgamma(1, shape=Q_shape, rate=Q_rate)
  
  # 5. Update sigma2
  residuals <- Y - aperm(apply(eta, 1:2, function(x) Lambda %*% x), c(2,3,1))
  ssr <- sum(residuals^2)
  sigma2 <- 1/rgamma(1, a0 + N*T*J/2, b0 + ssr/2)
  
  # Store samples
  if(iter %in% keep) {
    idx <- which(keep == iter)
    Lambda_samples[idx,,] <- Lambda
    sigma2_samples[idx] <- sigma2
    eta_samples[idx,,,] <- eta
    rho_samples[idx] <- rho
    Q_samples[idx] <- Q
  }
  
  if(iter %% 100 == 0) cat("Iteration:", iter, "\n")
}

# save.image("C:/Users/SEC/Desktop/research/25summer/may5th/MultiTimeTest.RData")

# ----------------------
# Posterior Analysis
# ----------------------
# Parameter estimates
Lambda_mean <- apply(Lambda_samples, 2:3, mean)
sigma2_mean <- mean(sigma2_samples)
rho_mean <- mean(rho_samples)
Q_mean <- mean(Q_samples)
eta_mean <- apply(eta_samples, 2:4, mean)

# 1. Parameter comparison
cat("\nParameter Recovery:\n")
cat("True sigma²:", sigma2_true, 
    "| Estimated sigma²:", round(mean(sigma2_samples), 3), "\n")
cat("True rho:", rho_true, 
    "| Estimated rho:", round(mean(rho_samples), 3), "\n")
cat("True Q:", Q_true, 
    "| Estimated Q:", round(mean(Q_samples), 3), "\n")

# 2. Trace plots
par(mfrow=c(2,1))
plot(rho_samples, type="l", main="Trace plot: rho")
abline(h = rho_true, col = 'red')
plot(Q_samples, type="l", main="Trace plot: Q")
abline(h = Q_true, col = 'red')
par(mfrow=c(1,1))
ts.plot(Lambda_samples[,3,3])

# 3. Heatmap
# Lambda
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
true <- Lambda_true
estimated <- Lambda_mean
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:K)
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda")

# Eta
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
true <- eta_true[,,2]
estimated <- eta_mean[,,2]
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:N)
colnames(compare_mat) <- paste0("column", 1:K)
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Eta")

# LL^T - median
tcross_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  tcross_array[ , , i] <- tcrossprod(Lambda_samples[i, , ])
}
# Compute the posterior median across iterations
Lambda_median <- apply(tcross_array, c(1, 2), median)
true <- tcrossprod(Lambda_true)
estimated <- Lambda_median
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


# Cov - median
Cov_array <- array(NA, dim = c(J, J, dim(Lambda_samples)[1]))
for (i in 1:(dim(Lambda_samples)[1])) {
  Cov_array[ , , i] <- 1/(1-rho_samples[i]^2) *tcrossprod(Lambda_samples[i, , ]) + Q_samples[i] * diag(J)
}
# Compute the posterior median across iterations
Cov_median <- apply(Cov_array, c(1, 2), median)
true <- 1/(1-rho_true^2) *tcrossprod(Lambda_true)+ Q_true * diag(J)
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



rownames(true) <- paste0("row",1:J)
colnames(true) <- paste0("column",1:J)

rownames(Cov_median) <- paste0("row",1:J)
colnames(Cov_median) <- paste0("column",1:J)

compare_mat <- true
compare_mat[upper.tri(compare_mat)] <- Cov_median[upper.tri(Cov_median)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)

Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate






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
