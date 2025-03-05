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
L <- 30
c <- 1 # DP parameter. 1, 10, 100 다 해보기

# matrices and vectors
alpha0 <- alpha0_grp <- matrix(0, nrow = n, ncol = J)
mu_jl <- matrix(0, nrow = L, ncol = J) # 모두 양수값이어야하는지?

# weight
w <- c()
v <- c()
v <- rbeta(L-1, shape1 = 1, shape2 = c)
w[1] <- v[1]
for(l in 2:(L-1)){
  w[l] <- v[l]*prod(1-v[1:(l-1)])
}
w[L] <- 1-sum(w)

pie(w, main = "weight probability")

# mu_jl
# prior of mu_jl is N(mu0,sigma0)
mu0 <- 4 # mean
sigma0 <- 2 # standard deviation
mu0_j <- rnorm(J, mean = mu0, sd = sigma0)

mu0_j

for(j in 1:J){
  mu_jl[,j] <- rnorm(L, mean = mu0_j[j], sd = sigma0) # row is OTU, column is truncated infinite mean
}

# alpha0_grp
for(j in 1:J){
  alpha0_grp[,j] <- sample(x = 1:L, size = n, prob = w, replace = T)
}

# alpha0
sigma_base = 1
for(i in 1:n){
  for(j in 1:J){
    group_index <- alpha0_grp[i,j]
    alpha0[i,j] <- rnorm(1, mean = mu_jl[group_index,j], sd = sigma_base)
  }
}


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

y_count <- floor(y_exp) # observed count data
y_count 

data <- y_count


alpha0
grp_num <- unique(as.vector(alpha0_grp))



# make array to matrix
alpha1 <- alpha[,,1]
alpha2 <- alpha[,,2]
alpha3 <- alpha[,,3]

par(mfrow = c(1,4))
for(i in 1:length(grp_num)){
  plot(density(alpha0[alpha0_grp == grp_num[i]]), main = paste("grp",i,", time 0"))
  plot(density(alpha1[alpha0_grp == grp_num[i]]), main = paste("grp",i,", time 1"))
  plot(density(alpha2[alpha0_grp == grp_num[i]]), main = paste("grp",i,", time 2"))
  plot(density(alpha3[alpha0_grp == grp_num[i]]), main = paste("grp",i,", time 3"))
}

# group 1 : density plot along the time
par(mfrow = c(1,4))
  plot(density(alpha0[1:10,]), main = paste("grp",1,", time 0"))
  plot(density(alpha1[1:10,]), main = paste("grp",1,", time 1"))
  plot(density(alpha2[1:10,]), main = paste("grp",1,", time 2"))
  plot(density(alpha3[1:10,]), main = paste("grp",1,", time 3"))

# group 2 : density plot along the time
par(mfrow = c(1,4))
  plot(density(alpha0[11:20,]), main = paste("grp",2,", time 0"))
  plot(density(alpha1[11:20,]), main = paste("grp",2,", time 1"))
  plot(density(alpha2[11:20,]), main = paste("grp",2,", time 2"))
  plot(density(alpha3[11:20,]), main = paste("grp",2,", time 3"))
  

