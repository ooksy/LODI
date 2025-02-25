################################################################################
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
m0 <- 0.5 ; C0 <- 0.8
phi_g1 <- 0.9 ; phi_g2 <- -0.3 
sig_w_t <- 1 # just give 1 to all time points

alpha0 <- matrix(rnorm(n*J, m0, C0), nrow = n, ncol = J)
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

library(mvtnorm)
eta <- array(0, dim=c(K, n, T))

eta0 <- rmvnorm(n, mean = rep(0, K), sigma = diag(K))
eta0 <- t(eta0)

eta[,,1] <- rho*eta0 + t(rmvnorm(n, mean = rep(0, K), sigma = diag(K)))


for(t in 2:T){
  eta[,,t] <- rho*eta[,,t-1] + t(rmvnorm(n, mean = rep(0, K), sigma = diag(K)))
}


# generating factor loading lambda
# Lambda는 관찰값, 시간별로 모두 같다고 가정
Lambda <- rmvnorm(J, mean = rep(0, K), sigma = diag(K))


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
y_exp

# exp(-10) ; exp(-5) ; exp(0) ; exp(3) ; exp(5) ; exp(10) ; exp(11)
# exp(0) 밑은 다 0으로 집계됨
# exp(10)만 되어도 굉장히 큼
# 실제 데이터는 4만이 거의 최대

y_count <- floor(y_exp) # observed count data
y_count 
# some values are too big
# 실제 데이터는 OTU별로 값이 비슷비슷함
# 실제 데이터는 0이 이렇게 많지 않음

data <- y_count

################################################################################
# Gibbs sampling ###############################################################
# 1. sampling the posterior latent variable ####################################
# sampling from truncated normal
# 0인 경우, normal distribution의 왼편에서 sampling
# post_latent0 <- array(0, dim=c(n, J, T))
# 
# # Initial value ?
# for(t in 1:T){
#   for(j in 1:J){
#     for(i in 1:n){
#       if(data[i,j,t] == 0){
#         #post_latent0[i,j,t] <- -abs(rnorm(1, mu[i,j,t], sqrt(sigma_sq)))  #change? absolute value of normal
#         u0  <- runif(1, min = pnorm(-5, 0, 1), max = pnorm(0, 0, 1)) # 그냥 pnorm(-5, 0, 1)로 고정
#         post_latent0[i,j,t] <- sqrt(sigma_sq)*qnorm(u0, 0, 1) + mu[i,j,t]
#       } else {
#         u <- runif(1, min = pnorm((log(data[i,j,t])-mu[i,j,t])/sqrt(sigma_sq), 0, 1), max = pnorm((log(data[i,j,t]+1)-mu[i,j,t])/sqrt(sigma_sq), 0, 1))
#         if(u == 1){u <- 1-1e-10} 
#         post_latent0[i,j,t] <- sqrt(sigma_sq)*qnorm(u, 0, 1) + mu[i,j,t]
#       }
#     }
#   }
# }

# log(0.01)
# plot(post_latent0, log(data+0.01))
# which(post_latent0-log(data+0.01)< -1, arr.ind = T) # 이상치 확인
# post_latent0[which(post_latent0-log(data+0.01)< -1, arr.ind = T)]
# log(data+0.01)[which(post_latent0-log(data+0.01)< -1, arr.ind = T)]
# (data+0.01)[which(post_latent0-log(data+0.01)< -1, arr.ind = T)]


# floor(exp(post_latent0))
# data



# gibbs sampling ###############################################################
# parameter setting
n_sim <- 5000
mu_updated <- array(0, dim=c(n, J, T))
post_latent <- array(0, dim=c(n, J, T, n_sim))
#post_latent[,,,1] <- post_latent0
r_it <- array(0, dim=c(n, T, n_sim))
alpha_ijt <- array(0, dim=c(n, J, T, n_sim))
lambda_j <- array(0, dim=c(K, J, n_sim))
eta_it <- array(0, dim=c(K, n, T, n_sim))
sigma_t <- matrix(NA, nrow = n_sim, ncol = 1) # variance. norm function 안에 들어갈때는 sqrt 해줘야함
sigma_t[1] <- 2 # initial value
m <- rep(0, T+1)
C <- rep(0, T+1)
m[1] <- 1 ; C[1] <- 0.8
a <- rep(0, T)
R <- rep(0, T)
K # prespecified?
phi1 <- rep(0, n_sim)
phi2 <- rep(0, n_sim)
phi1[1] <- 1 # initial value
phi2[1] <- 1 # initial value
m_vec <- matrix(0, nrow = T+1, ncol = K)
C_mat <- array(0, dim = c(K, K, T+1))
m_vec[1,] <- rep(0, K)
C_mat[,,1] <- diag(K)
a_vec <- matrix(0, nrow = T, ncol = K)
R_mat <- array(0, dim = c(K, K, T))
rho_est <- rep(0, n_sim)
rho_est[1] <- 1 # initial value
sig_phi <- 1
sig_w <- 1




for(s in 2:n_sim){
  
  # post_latent[,,,s] <- y
  # r_it[,,s] <- r[,1,]
  # alpha_ijt[,,,s] <- alpha
  # phi1[s] <- phi_g1
  # phi2[s] <- phi_g2
  # lambda_j[,,s] <- t(Lambda)
  # eta_it[,,,s] <- eta[,,]
  # rho_est[s] <- rho
  # sigma_t[s] <- sigma_sq  
  
  #sampling posterior latent variable ########################################
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
          if(u == 1){u <- 1-1e-10}
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
  
  
  # rho update (variable name is rho_est) ####################################
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

post_mean <- numeric()
for(i in 1:n){
  post_mean[i] <- mean(post_latent[i,1,1,n_sim]) 
}

cbind(post_mean, post_latent[,1,1,n_sim], log(data[,1,1]+0.01))
cbind("last samples"=post_latent[,1,1,n_sim],"true latent"= log(data[,1,1]+0.01))



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
  
  png(file = paste("C:/Users/SEC/Desktop/research/25-1/week2/r_i",t,".png", sep=""),
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
for(t in 1:T){
  alpha_mean[t] <- mean(alpha_ijt[1,1,t,])  
}
alpha_mean
cbind("post mean"=alpha_mean, "last value" = alpha_ijt[1,1,,n_sim], "true" = alpha[1,1,])


alpha_ijt[,,,n_sim] ; alpha

par(mar = c(1,1,1,1))

for(t in 1:T){
  
  png(file = paste("C:/Users/SEC/Desktop/research/25-1/week2/alpha_ij",t,".png", sep=""),
      width = 6,
      height = 5,
      units = "in",
      res = 1200)
  
  plot(alpha[,,t], alpha_ijt[,,t,n_sim], xlab = "true alpha",  ylab = "last samples", main = paste("normalized abundance at time point", t))
  abline(0, 1, col = "red")
  
  dev.off()
  
}


graphics.off()

ts.plot(alpha_ijt[1,1,1,])



r_it_estimated <- array(0, dim = c(n, J, T))
r_it_estimated[,,1] <- r_it[,1,n_sim]
r_it_estimated[,,2] <- r_it[,2,n_sim]
r_it_estimated[,,3] <- r_it[,3,n_sim]

plot(r + alpha, r_it_estimated[,,] + alpha_ijt[,,,n_sim], xlab = "true (r + alpha)", ylab = "estimated (r + alpha)", main = "r + alpha")
abline(0, 1, col = "red")




for(t in 1:T){
  
  png(file = paste("C:/Users/SEC/Desktop/research/25-1/week2/mu_ij",t,".png", sep=""),
      width = 6,
      height = 5,
      units = "in",
      res = 1200)
  
  plot(r[,,t] + alpha[,,t], r_it_estimated[,,t] + alpha_ijt[,,t,n_sim], xlab = "true (r + alpha)", ylab = "estimated (r + alpha)", main = paste("r + alpha at time point", t))
  abline(0, 1, col = "red")
  
  dev.off()
  
}


# phi check
phi1[n_sim] ; phi_g1
phi2[n_sim] ; phi_g2
ts.plot(phi1)
ts.plot(phi2)


# lambda check #################################################################
lambda_j[,,n_sim]
t(Lambda)
ts.plot(lambda_j[1,1,])

# # if (!require("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# #
# # BiocManager::install("ComplexHeatmap")
# 
# browseVignettes("ComplexHeatmap")
library(ComplexHeatmap)

# Identifiability issue 때문에 곱 자체를 비교해야함.. 그건 일정하기 때문.. why?
# lambda'lambda 비교
estimated <- t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim]
true <- Lambda %*% t(Lambda)

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
estimated_f <- t(lambda_j[,,n_sim])%*%eta_it[,,1,n_sim]
rownames(estimated_f) <- paste0("row",1:30)
colnames(estimated_f) <- paste0("column",1:20)
true_f <- Lambda%*%eta[,,1]
rownames(true_f) <- paste0("row",1:30)
colnames(true_f) <- paste0("column",1:20)
ht1 <- Heatmap(estimated_f, column_order = colnames(estimated_f), row_order = rownames(estimated_f),
               row_title = "OTU", column_title = "Estimated")
ht2 <- Heatmap(true_f, column_order = colnames(true_f), row_order = rownames(true_f),
               row_title = "OTU", column_title = "True")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Lambda*eta")

# eta check ####################################################################
eta_it[,,1,n_sim]
eta[,,1]



estimated_f <- eta_it[,,1,n_sim]
rownames(estimated_f) <- paste0("row",1:3)
colnames(estimated_f) <- paste0("column",1:20)
true_f <- eta[,,1]
rownames(true_f) <- paste0("row",1:3)
colnames(true_f) <- paste0("column",1:20)
ht1 <- Heatmap(estimated_f, column_order = colnames(estimated_f), row_order = rownames(estimated_f),
               row_title = "OTU", column_title = "Estimated")
ht2 <- Heatmap(true_f, column_order = colnames(true_f), row_order = rownames(true_f),
               row_title = "OTU", column_title = "True")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "Eta")


graphics.off()
ts.plot(eta_it[1,1,1,])
ts.plot(eta_it[2,1,1,])


# rho check ####################################################################
rho
rho_est[n_sim]
mean(rho_est[])

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

true0t1 <- sum(data[,,1] == 0)
true0t2 <- sum(data[,,2] == 0)
true0t3 <- sum(data[,,3] == 0)

cbind("time 1" = true0t1 , "time 2" =  true0t2 , "time 3" = true0t3)


num0ratio1 <- sum(floor(exp(r_it_estimated[,,1] + alpha_ijt[,,1,n_sim])) == 0)/length(data[,,1])
num0ratio2 <- sum(floor(exp(r_it_estimated[,,2] + alpha_ijt[,,2,n_sim])) == 0)/length(data[,,2])
num0ratio3 <- sum(floor(exp(r_it_estimated[,,3] + alpha_ijt[,,3,n_sim])) == 0)/length(data[,,3])

cbind("time 1" = num0ratio1 , "time 2" =  num0ratio2 , "time 3" = num0ratio3)


truenum0ratio1 <- sum(data[,,1] == 0)/length(data[,,1])
truenum0ratio2 <- sum(data[,,2] == 0)/length(data[,,2])
truenum0ratio3 <- sum(data[,,3] == 0)/length(data[,,3])

cbind("time 1" = truenum0ratio1 , "time 2" =  truenum0ratio2 , "time 3" = truenum0ratio3)





# covariance between time point check ##########################################
# Cov(Y_it) 
Cov_Y_tr <- (1/(1-rho^2))*Lambda%*%t(Lambda) + sigma_sq*diag(J)
Cov_Y_est <- (1/(1-rho_est[n_sim]^2))*t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim] + sigma_t[n_sim]*diag(J)

mean_cov <- matrix(0, J, J)
for(i in 1:n_sim){
  mean_cov <- mean_cov + t(lambda_j[,,i])%*%lambda_j[,,i]
}
mean_cov/n_sim
Cov_Y_est <- (1/(1-mean(rho_est[])^2))*(mean_cov/n_sim) + mean(sigma_t[])*diag(J)


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
        row_title = "lower : true", column_title = "upper : estimated", name = "Cov(Y_it)") # diagonal part is true value, lower - true, upper - posterior estimate



# time point1 and time point2
Cov_Y_tr <- (rho/(1-rho^2))*Lambda%*%t(Lambda)
Cov_Y_est <- (rho_est[n_sim]/(1-rho_est[n_sim]^2))*t(lambda_j[,,n_sim])%*%lambda_j[,,n_sim]

mean_cov <- matrix(0, J, J)
for(i in 1:n_sim){
  mean_cov <- mean_cov + t(lambda_j[,,i])%*%lambda_j[,,i]
}
mean_cov/n_sim

Cov_Y_est <- (mean(rho_est[])/(1-mean(rho_est[])^2))*(mean_cov/n_sim)



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
        row_title = "lower : true", column_title = "upper : estimated", name = "Cov(Y_it, Y_i(t+1))") # diagonal part is true value, lower - true, upper - posterior estimate

