rm(list = ls())

# library 
library(mvtnorm)
library(ComplexHeatmap)
library(gridtext)
library(truncnorm)

# set seed
set.seed(123)

# Data generation
# initial setting
n <- 20
J <- 30
K <- 3
T <- 3

# generating parameter
# factor loading
lambda <- matrix(0, nrow = J, ncol = K)
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
  eta[,,t] <- rho*eta[,,t-1] + t(rmvnorm(n, mean = rep(0, K), sigma = diag(sqrt(Q), K)))
}

# data y
y <- array(0, dim=c(n, J, T))

for(t in 1:T){
  y[,,t] <- t(eta[,,t]) %*% t(lambda) + rmvnorm(n, mean = rep(0, J), sigma = diag(sqrt(sigma2), J))
}

# MCMC setting 
n_sim <- 50000
burn_in <- n_sim/2
thin <- 2
keep <- seq(burn_in + 1, n_sim, by = thin)

# MCMC loop
for(s in 1:n_sim){
  # sampling lambda
  
  # sampling eta
  
  # sampling rho
  
  # sampling sigma2
}


# compare the result





