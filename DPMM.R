n <- 20
J <- 30
L <- 30
c <- 1 # DP parameter. 1, 10, 100 다 해보기

# matrices and vectors
alpha0 <- alpha0_grp <- matrix(0, nrow = n, ncol = J)
mu_jl <- matrix(0, nrow = L, ncol = J) # 모두 양수값이어야하는지?

# data generating ##############################################################
set.seed(250219)

# weight 
w <- c()
v <- c()
v <- rbeta(L-1, shape1 = 1, shape2 = c)
w[1] <- v[1]
for(l in 2:(L-1)){
  w[l] <- v[l]*prod(1-v[1:l-1])
}
w[L] <- 1-sum(w)

pie(w, main = "weight probability")

# mu_jl
# prior of mu_jl is N(mu0,sigma0)
mu0 <- 6 # mean
sigma0 <- 3 # standard deviation
mu0_j <- rnorm(J, mean = mu0, sd = sigma0)

mu0_j

for(j in 1:J){
  mu_jl[,j] <- rnorm(L, mean = mu0_j[j], sd = sigma0)
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

hist(alpha0[,6], probability = T) # histogram of OTU 6
hist(alpha0[,10], probability = T) # histogram of OTU 10


# DPMM estimation ##############################################################
# Using finite mixture distribution estimating 
# prefix L = 30
niter <- 1000
w_update <- matrix(0, nrow = niter, ncol = L)
w_update[1,] <- rep(1/L,L) # initial value
alpha0_grp_update <- array(data = 0, dim = c(n, J, niter))
alpha0_update <- array(data = 0, dim = c(n, J, niter))
mu_jl_update <- array(data = 0, dim = c(J, L, niter))

for(s in 2:niter){
  
  # update group membership
  grp_prob <- c()
  for(i in 1:n){
    for(j in 1:J){
      for(l in 1:L){
        grp_prob[l] <- w_update[s-1,l]*dnorm(alpha0_update[i,j,s-1], mean = mu_jl_update[j,l,s-1], sd = sigma_base)
      }
      alpha0_grp_update[i,j,s] <- sample(x=1:L, size = 1, replace = T, prob = grp_prob/sum(grp_prob))
    }
  }
  
  # update alpha0 | group
  for(i in 1:n){
    for(j in 1:J){
      grp_idx <- alpha0_grp_update[i,j,s]
      alpha0_update[i,j,s] <- rnorm(1, mean = mu_jl_update[j,grp_idx,s-1], sd = sigma_base)
    }
  }
  
  # update weight
  v <- c()
  
  for(l in 1:(L-1)){
    v[l] <- rbeta(1, shape1 = 1 + sum(alpha0_grp_update[,,s]==l), shape2 = c + sum(alpha0_grp_update[,,s]>l))
  }
  
  w_update[s,1] <- v[1]
  for(l in 2:(L-1)){
    w_update[s,l] <- v[l]*prod(1-v[1:l-1])
  }
  w_update[s,L] <- 1-sum(w_update[s,])
  

  # update mu_jl
  for(j in 1:J){
    for(l in 1:L){
      mu_sigma <- 1/(1/sigma0^2 + 1/sigma_base^2)
      alpha_current <- alpha0_update[,,s]
      mu_mean <- mu_sigma*(mu0/sigma0^2 + sum(alpha_current[alpha0_grp_update[,,s]==l])/sigma_base^2)
      mu_jl_update[j,l,s] <- rnorm(1, mean = mu_mean, sd = sqrt(mu_sigma))
    }
  }
}


# simulation result ############################################################
# alpha0
library(ComplexHeatmap)
alpha0 ; alpha0_update[,,niter]

rownames(alpha0) <- paste0("row",1:20)
colnames(alpha0) <- paste0("column",1:30)
rownames(alpha0_update[,,niter]) <- paste0("row",1:20)
colnames(alpha0_update[,,niter]) <- paste0("column",1:30)

ht1 <- Heatmap(alpha0, column_order = colnames(alpha0), row_order = rownames(alpha0),
               , column_title = "alpha0")
ht2 <- Heatmap(alpha0_update[,,niter], column_order = colnames(alpha0_update[,,niter]), row_order = rownames(alpha0_update[,,niter]),
               , column_title = "estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "alpha0")


# latent variable (group member)
alpha0_grp ; alpha0_grp_update[,,niter]

rownames(alpha0_grp) <- paste0("row",1:20)
colnames(alpha0_grp) <- paste0("column",1:30)
rownames(alpha0_grp_update[,,niter]) <- paste0("row",1:20)
colnames(alpha0_grp_update[,,niter]) <- paste0("column",1:30)

ht1 <- Heatmap(alpha0_grp, column_order = colnames(alpha0_grp), row_order = rownames(alpha0_grp),
               , column_title = "alpha0_group")
ht2 <- Heatmap(alpha0_grp_update[,,niter], column_order = colnames(alpha0_grp_update[,,niter]), row_order = rownames(alpha0_grp_update[,,niter]),
               , column_title = "estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "alpha0_grp")

# weight
round(w, digits = 2) ; round(w_update[niter,], digits = 2)
pie(round(w, digits = 2), main= "true weight")
pie(round(w_update[niter,], digits = 2), main= "estimated weight")


# mu_jl
mu_jl ; mu_jl_update[,,niter]

rownames(mu_jl) <- paste0("row",1:30)
colnames(mu_jl) <- paste0("column",1:30)
rownames(mu_jl_update[,,niter]) <- paste0("row",1:30)
colnames(mu_jl_update[,,niter]) <- paste0("column",1:30)

ht1 <- Heatmap(mu_jl, column_order = colnames(mu_jl), row_order = rownames(mu_jl),
               , column_title = "mu_jl")
ht2 <- Heatmap(mu_jl_update[,,niter], column_order = colnames(mu_jl_update[,,niter]), row_order = rownames(mu_jl_update[,,niter]),
               , column_title = "estimated")
ht_list <- ht1 + ht2
draw(ht_list, column_title = "mu_jl")


