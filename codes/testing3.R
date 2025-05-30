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
c <- 1
L <- 30
sigma_base <- 1
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
K <- 3
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
Lambda <- rmvnorm(J, mean = rep(0, K), sigma = diag(1/3, nrow = K))


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


################################################################################
# Gibbs sampling ###############################################################
# parameter setting
n_sim <- 10000
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
  for(i in 1:n){
    for(j in 1:J){
      # m, C, a, R update - no save for whole i, j
      # m, C는 initialization value가 1번째 row에 들어가서 index가 +1됨
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



# phi check ####################################################################
phi1[n_sim] ; phi_g1
phi2[n_sim] ; phi_g2
ts.plot(phi1)
ts.plot(phi2)


# lambda check #################################################################
lambda_j[,,n_sim]
t(Lambda)
ts.plot(lambda_j[1,1,])


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




# rho check ####################################################################
rho
rho_est[n_sim]

ts.plot(rho_est) 
mean(rho_est)


# shuangjie's method ###########################################################
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


Cor_Y_tr <- cov2cor((1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J))
Cor_Y_est <- matrix(0, J, J)
for(s in burn_lef){
  Cor_Y_est <- Cor_Y_est + cov2cor((1/(1-rho_est[s]^2))*t(lambda_j[,,s])%*%lambda_j[,,s] + sigma_t[s]*diag(J))
}
Cor_Y_est = Cor_Y_est/length(burn_lef)
rownames(Cov_Y_tr) <- paste0("row",1:J)
colnames(Cov_Y_tr) <- paste0("column",1:J)
rownames(Cor_Y_est) <- paste0("row",1:J)
colnames(Cor_Y_est) <- paste0("column",1:J)
compare_mat <- Cor_Y_tr
compare_mat[upper.tri(compare_mat)] <- Cor_Y_est[upper.tri(Cor_Y_est)]
compare_mat
rownames(compare_mat) <- paste0("row", 1:J)
colnames(compare_mat) <- paste0("column", 1:J)
Heatmap(compare_mat, column_order = colnames(compare_mat), row_order = rownames(compare_mat),
        row_title = "lower : true", column_title = gt_render("Cor(Y_it) <br> upper : estimated"), name = "mat") # diagonal part is true value, lower - true, upper - posterior estimate
