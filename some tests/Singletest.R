library(mvtnorm)   
library(truncnorm) 
library(circlize)

set.seed(123)

# Dimensions
n <- 100
J <- 30 
K <- 5

# True parameters
Lambda_true <- matrix(0, nrow = J, ncol = K)
for (j in 1:J) {
  for (k in 1:min(j, K)) {
    if (j == k) {
      Lambda_true[j, k] <- runif(1, min = 0.5, max = 1.5)
    } else {
      Lambda_true[j, k] <- rnorm(1, mean = 0, sd = 0.5)
    }
  }
}
sigma2_true <- 0.5

# Simulate latent factors
eta_true <- matrix(rnorm(n * K), nrow = n, ncol = K)

# Simulate observed data
Y <- eta_true %*% t(Lambda_true) + matrix(rnorm(n * J, sd = sqrt(sigma2_true)), n, J)



# Prior parameters
tau2 <- 1.0  # Prior variance for Lambda
a0 <- 2.0    # Prior shape for sigma^2
b0 <- 1.0    # Prior rate for sigma^2

# Initialization
Lambda <- matrix(0, nrow = J, ncol = K)
sigma2 <- 1.0
eta <- matrix(0, nrow = n, ncol = K)

# MCMC settings
n_iter <- 2000
burn_in <- n_iter/2

# Storage for posterior samples
Lambda_samples <- array(NA, dim = c(n_iter - burn_in, J, K))
sigma2_samples <- numeric(n_iter - burn_in)
eta_samples <- array(NA, dim = c(n_iter - burn_in, n, K))



for (iter in 1:n_iter) {
  # Update Lambda
  for (j in 1:J) {
    k_indices <- 1:min(j, K)
    E_j <- eta[, k_indices, drop = FALSE]
    y_j <- Y[, j]
    
    # Posterior calculations
    V_j_inv <- crossprod(E_j) / sigma2 + diag(1 / tau2, length(k_indices))
    V_j <- solve(V_j_inv)
    m_j <- V_j %*% (crossprod(E_j, y_j) / sigma2)
    
    # Sampling
    for (k_idx in seq_along(k_indices)) {
      k <- k_indices[k_idx]
      if (k < j) {
        Lambda[j, k] <- rnorm(1, mean = m_j[k_idx], sd = sqrt(V_j[k_idx, k_idx]))
      } else if (k == j) {
        Lambda[j, k] <- rtruncnorm(1, a = 0, mean = m_j[k_idx], sd = sqrt(V_j[k_idx, k_idx]))
      }
    }
  }
  
  # Update eta
  Lambda_t_Lambda <- crossprod(Lambda)
  V_eta_inv <- Lambda_t_Lambda / sigma2 + diag(K)
  V_eta <- solve(V_eta_inv)
  
  for (i in 1:n) {
    y_i <- Y[i, ]
    m_eta <- V_eta %*% (t(Lambda) %*% y_i / sigma2)
    eta[i, ] <- rmvnorm(1, mean = as.vector(m_eta), sigma = V_eta)
  }
  
  # Update sigma2
  residuals <- Y - eta %*% t(Lambda)
  ssr <- sum(residuals^2)
  a_post <- a0 + (n * J) / 2
  b_post <- b0 + ssr / 2
  sigma2 <- 1 / rgamma(1, shape = a_post, rate = b_post)
  
  # Store samples after burn-in
  if (iter > burn_in) {
    idx <- iter - burn_in
    Lambda_samples[idx, , ] <- Lambda
    sigma2_samples[idx] <- sigma2
    eta_samples[idx, , ] <- eta
  }
  
  # Optional: print progress
  if (iter %% 100 == 0) {
    cat("Iteration:", iter, "\n")
  }
}


# Posterior means
Lambda_mean <- apply(Lambda_samples, c(2, 3), mean)

Lambda_mean <- apply(Lambda_samples, c(2, 3), mean)
sigma2_mean <- mean(sigma2_samples)
eta_mean <- apply(eta_samples, c(2, 3), mean)

# Display results
print("Posterior mean of Lambda:")
print(Lambda_mean)


library(ComplexHeatmap)
true <- Lambda_true
estimated <- Lambda_mean
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:K)
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True")
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda")



library(ComplexHeatmap)
true <- Lambda_true
estimated <- Lambda_mean
compare_mat <- true
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:K)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht1 <- Heatmap(true, column_order = colnames(true), row_order = rownames(true),
               row_title = "OTU", column_title = "True", col = col_fun)
ht2 <- Heatmap(estimated, column_order = colnames(estimated), row_order = rownames(estimated),
               row_title = "OTU", column_title = "Estimated", col = col_fun)
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda")



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
draw(ht_list, column_title = "Lambda")






print("Posterior mean of sigma^2:")
print(sigma2_mean)
sigma2_true



