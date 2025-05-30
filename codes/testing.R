rm(list = ls())

# library #####################################################################
library(mvtnorm)
library(ComplexHeatmap)
library(gridtext)


# Data generation ##############################################################
set.seed(10000)
g1 <- 1:10
g2 <- 11:20
n <- max(g2)
J <- 30
T <- 3
y <- r <- alpha <- factor <- error <- array(0, dim=c(n, J, T))


# generating normalizing parameter r
# sigma_g1 = 1, sigma_g1 = 1
sig_g1 <- 1 # variance
sig_g2 <- 1
for(t in 1:T){
  for(i in g1){
    r[i, , t] <- runif(1, min = 0, max = 2)
  }
  for(i in g2){
    r[i, , t] <- runif(1, min = 0, max = 2)
  }
}

r

# generating normalized abundance alpha
# define truncate the infinite mixture to 30, DP concentration parameter c
# matrices and vectors
mu0 <- 0 # prior mean of mu_j0
sigma0 <- 5 # prior standard deviation of mu_jl
sigmaj0 <- 5 # prior standard deviation of mu_j0
sigmajl <- 5 # prior standard deviation of mu_jl
sigma_base <- 1 # same role with the sigmajl
c <- 1
L <- 30
alpha0 <- alpha0_grp <- matrix(0, nrow = n, ncol = J)


# group mean
grp_mean <- c(-4, 2, 8)
# exp(-4) ; exp(2) ; exp(8)

# group var
grp_var <- 1/2

# weight
w1 <- c(0.6, 0.2, 0.2)
w2 <- c(0.3, 0.5, 0.2)
w3 <- c(0.1, 0.2, 0.7)

# alpha0_grp
# J를 weight 그룹 별로 묶고, 그 weight 따라 안의 샘플들 그룹핑
wgrp1 <- sample(1:J,size = J/3)
wgrp2 <- sample(setdiff(1:J,wgrp1), size= J/3)
wgrp3 <- setdiff(setdiff(1:J,wgrp1), wgrp2)

wgrp1 ; wgrp2 ; wgrp3

for(i in wgrp1){
    alpha0_grp[,i] <- sample(1:3, size = n, replace = T, prob = w1) 
}

for(i in wgrp2){
    alpha0_grp[,i] <- sample(1:3, size = n, replace = T, prob = w2) 
}

for(i in wgrp3){
    alpha0_grp[,i] <- sample(1:3, size = n, replace = T, prob = w3) 
}

alpha0_grp

# alpha0
for(i in 1:n){
  for(j in 1:J){
    index <- alpha0_grp[i,j]
    alpha0[i,j] <- rnorm(1, mean = grp_mean[index], sd = grp_var)
  }
}

alpha0


phi_g1 <- 0.9 ; phi_g2 <- -0.3 
sig_w_t <- 1 # just give 1 to all time points


for(i in g1){
  for(j in 1:J){
    alpha[i,j,1] <- phi_g1*alpha0[i,j] + rnorm(1, 0, sqrt(sig_w_t))
  }
}

for(i in g2){
  for(j in 1:J){
    alpha[i,j,1] <- phi_g2*alpha0[i,j] + rnorm(1, 0, sqrt(sig_w_t))
  }
}


for(t in 2:T){
  for(i in g1){
    for(j in 1:J){
      alpha[i,j,t] <- phi_g1*alpha[i,j, t-1] + rnorm(1, 0, sqrt(sig_w_t))  
    }
  }
  for(i in g2){
    for(j in 1:J){
      alpha[i,j,t] <- phi_g2*alpha[i,j, t-1] + rnorm(1, 0, sqrt(sig_w_t))  
    }
  }
}


# generating factor eta
# eta 는 관찰값, 시간마다 다른 값을 가진다
K <- J/3
rho <- 0.7

eta <- array(0, dim=c(K, n, T))

eta0 <- rmvnorm(n, mean = rep(0, K), sigma = diag(K))
eta0 <- t(eta0)

eta[,,1] <- rho*eta0 + t(rmvnorm(n, mean = rep(0, K), sigma = diag(K)))


for(t in 2:T){
  eta[,,t] <- rho*eta[,,t-1] + t(rmvnorm(n, mean = rep(0, K), sigma = diag(K)))
}


# generating factor loading lambda
# Lambda는 관찰값, 시간별로 모두 같다고 가정
# Lambda <- rmvnorm(J, mean = rep(0, K), sigma = diag(K))
Lambda <- rmvnorm(J, mean = rep(0, K), sigma = diag(1, nrow = K))


# generating factor array
for(t in 1:T){
  factor[,,t] <- t(eta[,,t])%*%t(Lambda)
}

factor

# generating the error
# variance of error term is 1
sig_a = 2 ; sig_b = 2 # hyperparameter
sigma_sq <- 1/rgamma(1, sig_a, sig_b)
sigma_sq

# final data
mu <- r + alpha  # mean of the log of the latent variable(factor부분 안 더함?)
sigma_sq         # variance of the latent variable

for(i in 1:n){
  for(j in 1:J){
    for(t in 1:T){
      y[i,j,t] <- rnorm(1, mean = mu[i,j,t] + factor[i,j,t], sd = sqrt(sigma_sq))
    }
  }
}

y # log of the latent variable, follows normal distribution
y_exp <- exp(y) # original latent variable

y_count <- floor(y_exp) # observed count data

data <- y_count

cov_y <- array(0, dim = c(n,J,T))

for(i in 1:n){
  for(j in 1:J){
    for(t in 1:T){
      cov_y[i,j,t] <- y[i,j,t] - mu[i,j,t]
    }
  }
}

true_samcov1 <- cov(cov_y[,,1])
true_samcov2 <- cov(cov_y[,,2])
true_samcov3 <- cov(cov_y[,,3])


################################################################################
# Gibbs sampling ###############################################################
# parameter setting ############################################################
n_sim <- 1000
mu_updated <- array(0, dim=c(n, J, T))
post_latent <- array(0, dim=c(n, J, T, n_sim))
r_it <- array(0, dim=c(n, T, n_sim))
alpha_ijt <- array(0, dim=c(n, J, T, n_sim))
lambda_j <- array(0, dim=c(K, J, n_sim))
eta_it <- array(0, dim=c(K, n, T, n_sim))
sigma_t <- matrix(NA, nrow = n_sim, ncol = 1) # variance. norm function 안에 들어갈때는 sqrt 해줘야함
sigma_t[1] <- 2 # initial value
w_update <- matrix(0, nrow = n_sim, ncol = L)
w_update[1,] <- rep(1/L,L) # initial value
alpha0_grp_update <- array(data = 0, dim = c(n, J, n_sim))
alpha0_update <- array(data = 0, dim = c(n, J, n_sim))
mu_j0_update <- matrix(0, nrow = n_sim, ncol = J)
mu_jl_update <- array(data = 0, dim = c(L, J, n_sim))
mu_jl_update[,,1] <- rnorm(L*J, mean = 2, sd = 1) # initialize mu_jl with reasonable value # need to be updated with estimated value from the data
m <- rep(0, (T+1))
C <- rep(0, (T+1))
m[1] <- 1 ; C[1] <- 0.8
a <- rep(0, T)
R <- rep(0, T)
K # prespecified?
phi1 <- rep(0, n_sim)
phi2 <- rep(0, n_sim)
phi1[1] <- 1 # initial value
phi2[1] <- 1 # initial value
m_vec <- matrix(0, nrow = (T+1), ncol = K)
C_mat <- array(0, dim = c(K, K, (T+1)))
m_vec[1,] <- rep(0, K)
C_mat[,,1] <- diag(K)
a_vec <- matrix(0, nrow = T, ncol = K)
R_mat <- array(0, dim = c(K, K, T))
rho_est <- rep(0, n_sim)
rho_est[1] <- rbeta(1, shape1 = 1, shape2 = 0.5) # initial value
sig_phi <- 1
sig_w <- 1

# parameter initialize at the true value #######################################
# n_sim <- 1000
# mu_updated <- array(0, dim=c(n, J, T))
# post_latent <- array(0, dim=c(n, J, T, n_sim))
# r_it <- array(0, dim=c(n, T, n_sim))
# alpha_ijt <- array(0, dim=c(n, J, T, n_sim))
# phi1 <- rep(0, n_sim)
# phi2 <- rep(0, n_sim)
# lambda_j <- array(0, dim=c(K, J, n_sim))
# eta_it <- array(0, dim=c(K, n, T, n_sim))
# sigma_t <- matrix(NA, nrow = n_sim, ncol = 1) # variance. norm : sqrt 
# w_update <- matrix(0, nrow = n_sim, ncol = L)
# alpha0_grp_update <- array(data = 0, dim = c(n, J, n_sim))
# alpha0_update <- array(data = 0, dim = c(n, J, n_sim))
# mu_j0_update <- matrix(0, nrow = n_sim, ncol = J)
# mu_jl_update <- array(data = 0, dim = c(L, J, n_sim))
# mu_jl_update[,,1] <- rnorm(L*J, mean = 2, sd = 1) # initialize mu_jl with reasonable value # need to be updated with estimated value from the data
# m <- rep(0, (T+1))
# C <- rep(0, (T+1))
# a <- rep(0, T)
# R <- rep(0, T)
# m_vec <- matrix(0, nrow = (T+1), ncol = K)
# C_mat <- array(0, dim = c(K, K, (T+1)))
# a_vec <- matrix(0, nrow = T, ncol = K)
# R_mat <- array(0, dim = c(K, K, T))
# rho_est <- rep(0, n_sim)
# 
# 
# # initial value
# post_latent[,,,1] <- y
# r_it[,,1] <- r[,1,]
# alpha_ijt[,,,1] <- alpha
# phi1[1] <- phi_g1
# phi2[1] <- phi_g2
# lambda_j[,,1] <- t(Lambda)
# eta_it[,,,1] <- eta
# w_update[1,] <- rep(1/L,L)
# m[1] <- 1 ; C[1] <- 0.8
# m_vec[1,] <- rep(0, K) ; C_mat[,,1] <- diag(K)
# rho_est[1] <- rho # initial value
# sigma_t[1] <- sigma_sq # initial value
# 
# sig_phi <- 1
# sig_w <- 1

startTm <- Sys.time()
for(s in 2:n_sim){
  
  # post_latent[,,,s] <- y
  # r_it[,,s] <- r[,1,]
  # alpha0_grp_update[,,s] <- alpha0_grp
  # alpha_ijt[,,,s] <- alpha
  # alpha0_update[,,s] <- alpha0
  # alpha0_grp_update[,,s-1] <- alpha0_grp
  # alpha_ijt[,,,s-1] <- alpha
  # alpha0_update[,,s-1] <- alpha0
  # w_update[s,] <- w
  # mu_j0_update[s,] <- mu0_j
  # mu_jl_update[,,s-1] <- mu_jl
  # phi1[s-1] <- phi_g1
  # phi2[s-1] <- phi_g2
  # lambda_j[,,s-1] <- t(Lambda)
  # eta_it[,,,s-1] <- eta[,,]
  # rho_est[s-1] <- rho
  # rho_est[s] <- rho
  # sigma_t[s-1] <- sigma_sq
  
  # sampling posterior latent variable ########################################
  for(i in 1:n){
    for(j in 1:J){
      for(t in 1:T){
        mu_updated[i,j,t] <- r_it[i,t,s-1] + alpha_ijt[i,j,t,s-1] + lambda_j[,j,s-1]%*%eta_it[,i,t,s-1]
        if(data[i,j,t] == 0){
          post_latent[i,j,t,s] <- -abs(rnorm(1, mu_updated[i,j,t], sqrt(sigma_t[s-1])))  #change? absolute value of normal
          #u0  <- runif(1, min = pnorm(-5, 0, 1), max = pnorm(0, 0, 1))
          #post_latent[i,j,t,s] <- sqrt(sigma_t[s-1])*qnorm(u0, 0, 1) + mu_updated[i,j,t]
        } else {
          u <- runif(1, min = pnorm((log(data[i,j,t])-mu_updated[i,j,t])/sqrt(sigma_t[s-1]), 0, 1), max = pnorm((log(data[i,j,t]+1)-mu_updated[i,j,t])/sqrt(sigma_t[s-1]), 0, 1))
          ifelse(u == 0, u <- 1e-10, ifelse(u == 1, u<- 1-1e-10, u <- u))
          post_latent[i,j,t,s] <- sqrt(sigma_t[s-1])*qnorm(u, 0, 1) + mu_updated[i,j,t]
        }
      }
    }
  }
  
  
  
  
  
  # r update - 1iter 당 nT번, group별로 sigma 나눠야함! variance/sd인지 명확히 구분하기 ####
  for(i in g1){
    for(t in 1:T){
      # posterior mean
      r_mean <- 0
      for(j in 1:J){
        r_mean <- r_mean + (post_latent[i,j,t,s] - alpha_ijt[i,j,t,s-1] - lambda_j[,j,s-1]%*%eta_it[,i,t,s-1])
      }
      
      # posterior variance - var_part는 norm함수에서 sqrt 해주고, sig_g1는 variance임 그냥
      r_mean <- r_mean/sigma_t[s-1]
      r_var <- (J/sigma_t[s-1] + 1/sig_g1)
      r_mean <- r_mean/r_var
      r_var <- 1/r_var
      
      r_it[i,t,s] <- rnorm(1, mean = r_mean, sd = sqrt(r_var))
    }
  }
  
  for(i in g2){
    for(t in 1:T){
      # posterior mean
      r_mean2 <- 0
      for(j in 1:J){
        r_mean2 <- r_mean2 + (post_latent[i,j,t,s] - alpha_ijt[i,j,t,s-1] - lambda_j[,j,s-1]%*%eta_it[,i,t,s-1])
      }
      
      # posterior variance - var_part는 norm함수에서 sqrt 해주고, sig_g1는 variance임 그냥
      r_mean2 <- r_mean2/sigma_t[s-1]
      r_var2 <- (J/sigma_t[s-1] + 1/sig_g2)
      r_mean2 <- r_mean2/r_var2
      r_var2 <- 1/r_var2
      
      r_it[i,t,s] <- rnorm(1, mean = r_mean2, sd = sqrt(r_var2))
    }
  }
  
  
  
  
  # alpha update #############################################################
  # update group membership
  grp_prob <- c()
  for(i in 1:n){
    for(j in 1:J){
      for(l in 1:L){
        grp_prob[l] <- w_update[s-1,l]*dnorm(alpha0_update[i,j,s-1], mean = mu_jl_update[l,j,s-1], sd = sigma_base)
      }

      if(all(grp_prob == rep(0, L))) {grp_prob <- rep(1/L, L)}

      alpha0_grp_update[i,j,s] <- sample(x=1:L, size = 1, replace = T, prob = grp_prob/sum(grp_prob))
    }
  }

  # update alpha0 | group
  # time point T update
  for(i in 1:n){
    for(j in 1:J){
      # draw alpha0
      grp_idx <- alpha0_grp_update[i,j,s]

      m[1] <- mu_jl_update[grp_idx,j,s-1]
      C[1] <- sigma_base

      for(t in 1:T){
        if(i %in% g1){
          a[t] <- phi1[s-1]*m[t] ; R[t] <- (phi1[s-1]^2)*C[t] + sig_w_t
        } else {
          a[t] <- phi2[s-1]*m[t] ; R[t] <- (phi2[s-1]^2)*C[t] + sig_w_t
        }

        m[t+1] <- (a[t]/R[t] + (post_latent[i,j,t,s] - r_it[i,t,s] - lambda_j[,j,s-1]%*%eta_it[,i,t,s-1])/sigma_t[s-1])/(1/R[t] + 1/sigma_t[s-1])
        C[t+1] <- 1/(1/R[t] + 1/sigma_t[s-1])
      }

      alpha_ijt[i,j,T,s] <- rnorm(1, mean = m[T+1], sd = sqrt(C[T+1])) # first draw a_ijT

      # draw a_ij(T-1)
      for(t in (T-1):1){
        if(i %in% g1){
          alpha_mean <- (m[t+1]/C[t+1] + phi1[s-1]*alpha_ijt[i,j,t+1,s]/sig_w_t) # sig_w_t is known, fixed
          alpha_var <- (1/C[t+1] + phi1[s-1]^2/sig_w_t)
          alpha_mean <- alpha_mean/alpha_var
          alpha_var <- 1/alpha_var
        } else {
          alpha_mean <- (m[t+1]/C[t+1] + phi2[s-1]*alpha_ijt[i,j,t+1,s]/sig_w_t) # sig_w_t is known, fixed
          alpha_var <- (1/C[t+1] + phi2[s-1]^2/sig_w_t)
          alpha_mean <- alpha_mean/alpha_var
          alpha_var <- 1/alpha_var
        }

        alpha_ijt[i,j,t,s] <- rnorm(1, mean = alpha_mean, sd = sqrt(alpha_var))
      }

      # draw alpha0
      if(i %in% g1){
        alpha0_mean <- (m[1]/C[1] + phi1[s-1]*alpha_ijt[i,j,1,s]/sig_w_t)
        alpha0_var <- (1/C[1] + phi1[s-1]^2/sig_w_t)
        alpha0_mean <- alpha0_mean/alpha0_var
        alpha0_var <- 1/alpha0_var
      } else {
        alpha0_mean <- (m[1]/C[1] + phi2[s-1]*alpha_ijt[i,j,1,s]/sig_w_t)
        alpha0_var <- (1/C[1] + phi2[s-1]^2/sig_w_t)
        alpha0_mean <- alpha0_mean/alpha0_var
        alpha0_var <- 1/alpha0_var
      }

      alpha0_update[i,j,s] <- rnorm(1, mean = alpha0_mean, sd = sqrt(alpha0_var))



    }
  }

  # update weight
  v <- c()

  for(l in 1:(L-1)){
    v[l] <- rbeta(1, shape1 = 1 + sum(alpha0_grp_update[,,s]==l), shape2 = c + sum(alpha0_grp_update[,,s]>l))
  }

  w_update[s,1] <- v[1]
  for(l in 2:(L-1)){
    w_update[s,l] <- v[l]*prod(1-v[1:(l-1)])
  }
  if(sum(w_update[s,1:(L-1)]) >= 1) w_update[s,L] <- 0 else w_update[s,L] <- 1-sum(w_update[s,1:(L-1)])

  # # update mu_jl
  # for(j in 1:J){
  #   for(l in 1:L){
  #     mu_sigma <- (1/(sigma0^2) + sum(alpha0_grp_update[,j,s]==l)/(sigma_base^2))
  #     mu_sigma <- 1/mu_sigma
  #     alpha_current <- alpha0_update[,,s]
  #     mu_mean <- mu_sigma*(mu0/(sigma0^2) + sum(alpha_current[alpha0_grp_update[,j,s]==l])/(sigma_base^2))
  #     mu_jl_update[l,j,s] <- rnorm(1, mean = mu_mean, sd = sqrt(mu_sigma))
  #   }
  # }

  # new update mu_j0 and mu_jl
  for(j in 1:J){
    mu0_sigma <- (1/(sigma0^2)) + (L/sigmaj0^2)
    mu0_sigma <- 1/mu0_sigma
    mu0_mean <- (mu0/(sigma0^2)) + (sum(mu_jl_update[,j,s-1])/sigmaj0^2)
    mu0_mean <- mu0_sigma*mu0_mean

    mu_j0_update[s,j] <- rnorm(1, mean = mu0_mean, sd = sqrt(mu0_sigma))

    for(l in 1:L){
      mu_jl_sigma <- (1/(sigmaj0^2)) + (sum(alpha0_grp_update[,j,s]==l)/(sigmajl^2))
      mu_jl_sigma <- 1/mu_jl_sigma
      mu_jl_mean <- (mu_j0_update[s,j]/(sigmaj0^2)) + (sum(alpha0_update[alpha0_grp_update[,j,s]==l,j,s])/(sigmajl^2))
      mu_jl_mean <- mu_jl_sigma*mu_jl_mean

      mu_jl_update[l,j,s] <- rnorm(1, mean = mu_jl_mean, sd = sqrt(mu_jl_sigma))
    }
  }
  
  
  
  # phi1 update ##############################################################
  alpha_term <- 0
  alpha_sq <- 0
  mu_phi1 <- 0
  sig_phi1 <- 0
  for(i in g1){
    for(j in 1:J){
      for(t in 2:T){
        alpha_term <- alpha_term + alpha_ijt[i,j,t,s]*alpha_ijt[i,j,t-1,s]
        alpha_sq <- alpha_sq + alpha_ijt[i,j,t-1,s]^2
      }
    }
  }
  mu_phi1 <- (alpha_term/sig_w)/(alpha_sq/sig_w + sig_phi)
  sig_phi1 <- 1/(alpha_sq/sig_w + sig_phi)
  
  u <- runif(1, min = pnorm((-1 - mu_phi1)/sqrt(sig_phi1), 0, 1), max = pnorm((1 - mu_phi1)/sqrt(sig_phi1), 0, 1))
  phi1[s] <- sqrt(sig_phi1)*qnorm(u, 0, 1) + mu_phi1
  
  # phi2 update ##############################################################
  alpha_term2 <- 0
  alpha_sq2 <- 0
  mu_phi2 <- 0
  sig_phi2 <- 0
  for(i in g2){
    for(j in 1:J){
      for(t in 2:T){
        alpha_term2 <- alpha_term2 + alpha_ijt[i,j,t,s]*alpha_ijt[i,j,t-1,s]
        alpha_sq2 <- alpha_sq2 + alpha_ijt[i,j,t-1,s]^2
      }
    }
  }
  mu_phi2 <- (alpha_term2/sig_w)/(alpha_sq2/sig_w + sig_phi)
  sig_phi2 <- 1/(alpha_sq2/sig_w + sig_phi)
  
  u <- runif(1, min = pnorm((-1 - mu_phi2)/sqrt(sig_phi2), 0, 1), max = pnorm((1 - mu_phi2)/sqrt(sig_phi2), 0, 1))
  phi2[s] <- sqrt(sig_phi2)*qnorm(u, 0, 1) + mu_phi2
  
  
  # lambda update ############################################################
  for(j in 1:J){
    
    eta_term <- 0
    mean_term <- 0
    
    for(t in 1:T){
      eta_term <- eta_term + eta_it[,,t,s-1]%*%t(eta_it[,,t,s-1])
      mean_term <- mean_term + eta_it[,,t,s-1]%*%(post_latent[,j,t,s] - r_it[,t,s] - alpha_ijt[,j,t,s])
    }
    eta_term <- eta_term/sigma_t[s-1]
    mean_term <- mean_term/sigma_t[s-1]
    
    lambda_mean <- solve(diag(K) + eta_term)%*%mean_term
    lambda_cov <- solve(diag(K) + eta_term)
    lambda_cov[lower.tri(lambda_cov)] <- t(lambda_cov)[lower.tri(lambda_cov)]
    
    lambda_j[,j,s] <- rmvnorm(1 , mean =  lambda_mean, sigma = lambda_cov)
  }
  
  
  
  
  
  # eta update #################################################################
  # forward filtering
  for(i in 1:n){
    # m, C, a, R update - no save for whole i, j
    # m, C는 initialization value가 1번째 row에 들어가서 index가 +1됨
    for(t in 1:T){
      a_vec[t,] <- rho_est[s-1]*m_vec[t,]
      R_mat[,,t] <- (rho_est[s-1]^2)*C_mat[,,t] + diag(K)
      
      R_mat_inv <- matrix(0, nrow = K, ncol = K)
      
      R_mat_inv <- solve(R_mat[,,t])
      # R_mat_inv[lower.tri(R_mat_inv)] <- t(R_mat_inv)[lower.tri(R_mat_inv)]
      
      m_vec[t+1,] <- solve(R_mat_inv + lambda_j[,,s]%*%t(lambda_j[,,s])/sigma_t[s-1])%*%
        ((lambda_j[,,s]%*%(post_latent[i,,t,s] - r_it[i,t,s]*rep(1,J) - alpha_ijt[i,,t,s])/sigma_t[s-1]) +
           R_mat_inv%*%a_vec[t,])
      
      C_mat[,,t+1] <- solve(R_mat_inv + lambda_j[,,s]%*%t(lambda_j[,,s])/sigma_t[s-1])
      C_mat[,,t+1][lower.tri(C_mat[,,t+1])] <- t(C_mat[,,t+1])[lower.tri(C_mat[,,t+1])]
    }
    
    
    eta_it[,i,T,s] <- rmvnorm(1, mean = m_vec[T+1,], sigma = C_mat[,,T+1]) # first draw eta_iT
    
    # draw eta_i(T-1)...
    for(t in (T-1):1){
      C_mat_inv <- matrix(0, nrow = K, ncol = K)
      
      C_mat_inv <- solve(C_mat[,,t+1])
      # C_mat_inv[lower.tri(C_mat_inv)] <- t(C_mat_inv)[lower.tri(C_mat_inv)]
      
      eta_mean <- solve((rho_est[s-1]^2)*diag(K) + C_mat_inv)%*%(eta_it[,i,t+1,s]/rho_est[s-1] + C_mat_inv%*%m_vec[t+1,])
      eta_cov <- solve((rho_est[s-1]^2)*diag(K) + C_mat_inv)
      eta_cov[lower.tri(eta_cov)] <- t(eta_cov)[lower.tri(eta_cov)]
      
      eta_it[,i,t,s] <- rmvnorm(1, mean = eta_mean, sigma = eta_cov)
    }
  }
  
  
  # rho update - gibbs sampling (overestimate) #################################
  rho_mean <- 0
  rho_var <- 0
  
  for(i in 1:n){
    for(t in 2:T){
      rho_var <- rho_var + as.numeric((t(eta_it[,i,t-1,s]) %*% eta_it[,i,t-1,s]))
    }
  }
  
  for(i in 1:n){
    for(t in 2:T){
      rho_mean <- rho_mean + as.numeric((t(eta_it[,i,t-1,s]) %*% eta_it[,i,t,s]))
    }
  }
  
  var_rho <- 1/rho_var
  mean_rho <- rho_mean * var_rho
  
  u2 <- runif(1, pnorm(0, mean = mean_rho, sd = sqrt(var_rho)), pnorm(1, mean = mean_rho, sd = sqrt(var_rho)))
  
  rho_est[s] <- sqrt(var_rho)*qnorm(u2, 0, 1) + mean_rho
  

  # sigma_t update ###########################################################
  lik_term <- 0
  for(i in 1:n){
    for(j in 1:J){
      for(t in 1:T){
        lik_term <- lik_term + (post_latent[i,j,t,s] - r_it[i,t,s] - alpha_ijt[i,j,t,s] - lambda_j[,j,s]%*%eta_it[,i,t,s])^2
      }
    }
  }
  
  
  sigma_t[s] <- 1/rgamma(1, sig_a + (n*J*T/2), sig_b + (lik_term/2))
  
}
endTm <- Sys.time()

endTm - startTm



burn_lef = (n_sim/2+1):n_sim
mean(rho_est[burn_lef]); median(rho_est[burn_lef]); rho
ts.plot(rho_est)



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
# post latent variable로 만든 data의 log값, data의 log값 비교
# latent variable의 log값이 0보다 클 때, 일직선을 따름
# latent variable의 log값이 -5에 가까울 때(아마 데이터값이 0인듯), 실제 데이터의 로그값은 -5에서 10까지 따름..
# log 데이터 값이 -5에서 10인데 latent variable이 -5라는 것은... posterior가 잘못 derive된 것은 아닌지..



# r_it check ###################################################################
r_it[,,n_sim] ; r[,1,]
plot(r[,1,], r_it[,,n_sim] , xlab = "true values", ylab = "last samples", main = "r_it")
abline(0, 1, col = "red")

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
ts.plot(lambda_j[1,2,])
ts.plot(lambda_j[1,3,])
ts.plot(lambda_j[1,4,])


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


# covariance between time point check ##########################################
Cov_Y_tr <- (1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J)
Cov_Y_est <- matrix(0, J, J)
for(s in burn_lef){
  Cov_Y_est <- Cov_Y_est + (1/(1-rho_est[s]^2))*t(lambda_j[,,s])%*%lambda_j[,,s] + sigma_t[s]*diag(J)
}
Cov_Y_est = Cov_Y_est/length(burn_lef)
rownames(Cov_Y_tr) <- paste0("row",1:J)
colnames(Cov_Y_tr) <- paste0("column",1:J)
rownames(Cov_Y_est) <- paste0("row",1:J)
colnames(Cov_Y_est) <- paste0("column",1:J)
compare_mat <- Cov_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


teta_Y_tr <- Lambda%*%t(Lambda)
teta_Y_est <- matrix(0, J, J)
for(s in burn_lef){
  teta_Y_est <- teta_Y_est + t(lambda_j[,,s])%*%lambda_j[,,s]
}
teta_Y_est = teta_Y_est/length(burn_lef)
rownames(teta_Y_tr) <- paste0("row",1:J)
colnames(teta_Y_tr) <- paste0("column",1:J)
rownames(teta_Y_est) <- paste0("row",1:J)
colnames(teta_Y_est) <- paste0("column",1:J)
compare_mat <- teta_Y_tr
compare_mat[upper.tri(compare_mat)] <- teta_Y_est[upper.tri(teta_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Lam * t(Lam) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


# true covariance and estimated covariance
SCov_Y_tr <- true_samcov1 # we can compare with true_samcov2/true_samcov3
Cov_Y_est <- matrix(0, J, J)
for(s in burn_lef){
  Cov_Y_est <- Cov_Y_est + (1/(1-rho_est[s]^2))*t(lambda_j[,,s])%*%lambda_j[,,s] + sigma_t[s]*diag(J)
}
Cov_Y_est = Cov_Y_est/length(burn_lef)
rownames(SCov_Y_tr) <- paste0("row",1:J)
colnames(SCov_Y_tr) <- paste0("column",1:J)
rownames(Cov_Y_est) <- paste0("row",1:J)
colnames(Cov_Y_est) <- paste0("column",1:J)
compare_mat <- SCov_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true sample cov", column_title = gt_render("Cov(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate


Y_mu <- array(0, dim = c(n,J,T))
r_it_array <- array(0, dim = c(n,J,T,n_sim))

for(s in 1:n_sim){
r_it_array[,,1,s] <- matrix(r_it[,1,s], nrow = n, ncol = J)
r_it_array[,,2,s] <- matrix(r_it[,2,s], nrow = n, ncol = J)
r_it_array[,,3,s] <- matrix(r_it[,3,s], nrow = n, ncol = J)
}
r_it_array

for(i in 1:n){
  for(j in 1:J){
    for(t in 1:T){
      for(s in burn_lef){
        Y_mu <- Y_mu + post_latent[,,,s] - r_it_array[,,,s] - alpha_ijt[,,,s]
      }
    }
  }
}


Y_mu <- post_latent[,,,n_sim] - r_it_array[,,,n_sim] - alpha_ijt[,,,n_sim]

SCov_Y_est <- cov(Y_mu[,,1])
rownames(SCov_Y_est) <- paste0("row",1:J)
colnames(SCov_Y_est) <- paste0("column",1:J)
rownames(SCov_Y_tr) <- paste0("row",1:J)
colnames(SCov_Y_tr) <- paste0("column",1:J)
compare_mat <- SCov_Y_tr
compare_mat[upper.tri(compare_mat)] <- SCov_Y_est[upper.tri(SCov_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true sample cov", column_title = gt_render("upper : estimated sample cov"), name = "mat") 


# true sample covariance vs true cov
SCov_Y_tr <- true_samcov1 # we can compare with true_samcov2/true_samcov3
Cov_Y_tr <- (1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J)
rownames(SCov_Y_tr) <- paste0("row",1:J)
colnames(SCov_Y_tr) <- paste0("column",1:J)
rownames(Cov_Y_tr) <- paste0("row",1:J)
colnames(Cov_Y_tr) <- paste0("column",1:J)
compare_mat <- SCov_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cov_Y_est[upper.tri(Cov_Y_tr)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true sample covariance", column_title = gt_render("Cov(Y_it) <br> upper : true"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate

