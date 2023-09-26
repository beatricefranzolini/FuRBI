#libraries
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
#Generate bivariate data 
set.seed(2019)
#Independent case
n <- 100
n1 <- 20
n2 <- n
mu1_true <- 10
mu2_true <- -10
tau_true <- 1
Y <- matrix(0, ncol = 2, nrow = n)
Y[,1] <- rnorm(n,mu1_true,tau_true)
Y[,2] <- rnorm(n,mu2_true,tau_true)
Y1 <- Y[1:20,1]
Y2 <- Y[,2]
grid <- seq(-15,15, by = 0.1) #evaluation grid

theta <- 1
N <- 150
iterations <- 5000 #iterations
burnin <- 2000

#univariate sampler
init <- list()
init$mu0 <- 0
init$tau0 <- 1
init$sigma <- 1
init$theta <- 1
init$Z <- rnorm(N, mean = init$mu0, sd = sqrt(init$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#independent case
res_univ <- sampler(Y1,init, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
# density estimation
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}

#bivariate initializer
N <- 150
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p
iterations <- 5000 #iterations
burnin <- 2000


#bivariate run
res <- sampler2(Y1,Y2, init, iterations, burnin)

#bivariate density estimation
Z <- res$Z1
Z2 <- res$Z2
p <- res$p
post <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post[j] <- mean(lik)
}
plot(res$rho, type = 'l')

#univariate sampler
init <- list()
init$mu0 <- 0
init$tau0 <- 1
init$sigma <- 1
init$theta <- 1
init$Z <- rnorm(N, mean = init$mu0, sd = sqrt(init$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#Exchangeable case
res_univ2 <- sampler(c(Y1,Y2),init, iterations, burnin)
Z_univ2 <- res_univ2$Z
p_univ2 <- res_univ2$p
post_univ2 <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ2[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ2[j] <- mean(lik)
}

#hierarchical
#hierarchical
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p


res <- sampler_hier_trunc(Y1, Y2, init_hier, iterations,burnin)
Z <- res$Z1
p <- res$p1
post_hier <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_hier[j] <- mean(lik)
}

#true density
true <- rep(0, length(grid))
for(i in 1:length(grid)){
  true[i] <- mean(dnorm(grid[i], mean = mu1_true, sd = sqrt(tau_true)))
}
true2 <- rep(0, length(grid))
for(i in 1:length(grid)){
  true2[i] <- mean(dnorm(grid[i], mean = mu2_true, sd = sqrt(tau_true)))
}

#to_use <- data.frame("Grid" = grid,"True" = true, "Exch" = post_univ2,
                     #"Univ." = post_univ, "FuRBI" = post,
                     #"Hier." = post_hier)
#saveRDS(to_use, file = "Density_firstplot.rds")

#plot
plot(grid,true,col="black", xlab = "Support",ylab = "Density",main = "",pch = 17,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
#plot(grid,true2,col="black", xlab = "Support",ylab = "Density",main = "",pch = 17,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.3, lwd = 1,lty = 2,las = 1,cex.lab = 1.4, ann = F, axes = F)
#par(new = TRUE)
plot(grid,post_univ2,col="gray", xlab = "Support",ylab = "Density",main = "",pch = 15,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,post_univ,col="brown", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,post,col="blue", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 1,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,post_hier,col="red", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-15,15),ylim = c(0,0.4),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 1,las = 1,cex.lab = 1.4, ann = F, axes = F)
legend(-9,0.47, legend = c("True", "Exch.", "Ind.","FuRBI", "Hier."), col = c("black", "grey","brown", "blue", "red"), lwd = 4, pch = c(17,15,19,13,16),lty = c(2,9,10,1,1), cex = 0.9, bty = "n")

#2) Extensive simulation to get the plot with different means
#normal vs normal

#GENERATE DATA
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
set.seed(2019)
mu_true = 10
n <- 100
n1 <- 20
n2 <- n
theta <- 1
Y1 <- rnorm(n1, mu_true, 1)
Y2 <- rnorm(n,-mu_true,1)
interval <- seq(-6,26,by = 2)#c(-5, 0, 5, 10, 15, 20, 25)#-6:26
grid <- seq(-25,25, by = 0.1) #evaluation grid
burnin <- 3000
iterations <- burnin+1000 #iterations
N <- 150

#SELECT CONCENTRATION PARAMETER HDP
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

#true
true <- rep(0,length(grid))
for(i in 1:length(grid)){
  true[i] <- dnorm(grid[i],mu_true,1)
}

#INITIALIZER
#independent
#univariate sampler
init_univ <- list()
init_univ$mu0 <- 0
init_univ$tau0 <- 1
init_univ$sigma <- 1
init_univ$theta <- 1
init_univ$Z <- rnorm(N, mean = init_univ$mu0, sd = sqrt(init_univ$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init_univ$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_univ$p <- p
res_univ <- sampler(Y1,init_univ, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}
perf_ind <- rep(sum(abs(true-post_univ)),length(interval))

#bivariate initializer
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#hierarchical initializer
init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p

n_samples <- 50 #number of samples to average
perf_ANRMI_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_exch_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_hier_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_univ_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
for(l in 1:n_samples){
  Y1 <- rnorm(n1, mu_true, 1)
  Y2 <- rnorm(n,-mu_true,1)
  for(k in 1:length(interval)){
  print(k)
  Y2_new <- Y2 + interval[k]
  #FuRBI
  res <- sampler2(Y1,Y2_new, init, iterations, burnin)
  Z <- res$Z1
  p <- res$p
  post <- rep(0, length(grid))
  for(j in 1:length(grid)){
    lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
    lik <- sapply(1:(iterations-burnin), function(i){
      return(sum(p[i,]*lik[i,]))
    })
    #lik <- diag(p%*%t(lik))
    post[j] <- mean(lik)
  }
  perf_ANRMI_matrix[l, k] <- sum(abs(post-true))
  #exchangeable
  res_univ2 <- sampler(c(Y1,Y2_new),init_univ, iterations, burnin)
  Z_univ2 <- res_univ2$Z
  p_univ2 <- res_univ2$p
  post_univ2 <- rep(0, length(grid))
  for(j in 1:length(grid)){
    lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init_univ$sigma))
    lik <- sapply(1:(iterations-burnin), function(i){
      return(sum(p_univ2[i,]*lik[i,]))
    })
    #lik <- diag(p%*%t(lik))
    post_univ2[j] <- mean(lik)
  }

  perf_exch_matrix[l, k] <- sum(abs(post_univ2-true))
  
  #hierarchical
  res <- sampler_hier_trunc(Y1, Y2_new, init_hier, iterations,burnin)
  Z <- res$Z1
  p <- res$p1
  post_hier <- rep(0, length(grid))
  for(j in 1:length(grid)){
    lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
    lik <- sapply(1:(iterations-burnin), function(i){
      return(sum(p[i,]*lik[i,]))
    })
    #lik <- diag(p%*%t(lik))
    post_hier[j] <- mean(lik)
  }
  perf_hier_matrix[l, k] <- sum(abs(post_hier-true))
  
  #univariate
  res_univ <- sampler(c(Y1),init_univ, iterations, burnin)
  Z_univ <- res_univ$Z
  p_univ <- res_univ$p
  post_univ <- rep(0, length(grid))
  for(j in 1:length(grid)){
    lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
    lik <- sapply(1:(iterations-burnin), function(i){
      return(sum(p_univ[i,]*lik[i,]))
    })
    #lik <- diag(p%*%t(lik))
    post_univ[j] <- mean(lik)
  }
  
  perf_univ_matrix[l, k] <- sum(abs(post_univ-true))
  
  to_print <- c(perf_ANRMI_matrix[l,k], perf_exch_matrix[l,k],perf_univ_matrix[l,k],perf_ind[k], perf_hier_matrix[l,k])
  names(to_print) <- c("FuRBI", "Exch.","Ind.", "Ind.", "Hier.")
  print(to_print)
  }
}
perf_ANRMI <- apply(perf_ANRMI_matrix, 2, median)
perf_exch <- apply(perf_exch_matrix, 2, median)
perf_univ <- apply(perf_univ_matrix, 2, median)
perf_hier <- apply(perf_hier_matrix, 2, median)
#extract errors
Errors <- matrix(c(-10+interval, 0.1*perf_exch, 0.1*perf_ind, 0.1*perf_ANRMI, 0.1*perf_hier), ncol = 5, byrow = F )
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(0,2),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="black", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(0,2),cex.axis = 1.6, cex = 1, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="gray", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(0,2),cex.axis = 1.6, cex = 1, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(0,2),cex.axis = 1.6, cex = 1, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)
#saveRDS(to_use, file = "Errors.rds")

#3) Other simulations: t-student vs exponential

#GENERATE DATA
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
set.seed(2019)
mu_true = 10
n <- 100
n1 <- 20
n2 <- n
theta <- 1
Y1 <- mu_true + rt(n1, df = 3)#rnorm(n1, mu_true, 1)
Y2 <- -mu_true+rexp(n2)-1#rnorm(n,-mu_true,1)
interval <- seq(-6,26,by = 2)#c(-5, 0, 5, 10, 15, 20, 25)#-6:26
grid <- seq(-25,25, by = 0.1) #evaluation grid
burnin <- 3000
iterations <- burnin+1000 #iterations
N <- 150

#SELECT CONCENTRATION PARAMETER HDP
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

#true
true <- rep(0,length(grid))
for(i in 1:length(grid)){
  true[i] <- dt(grid[i]-mu_true,df = 3)
}

#INITIALIZER
#independent
#univariate sampler
init_univ <- list()
init_univ$mu0 <- 0
init_univ$tau0 <- 1
init_univ$sigma <- 1
init_univ$theta <- 1
init_univ$Z <- rnorm(N, mean = init_univ$mu0, sd = sqrt(init_univ$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init_univ$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_univ$p <- p
res_univ <- sampler(Y1,init_univ, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}
perf_ind <- rep(sum(abs(true-post_univ)),length(interval))

#bivariate initializer
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#hierarchical initializer
init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p

n_samples <- 50 #number of samples to average
perf_ANRMI_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_exch_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_hier_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_univ_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
for(l in 1:n_samples){
  print(paste("Sample:", l))
  Y1 <- mu_true + rt(n1, df = 3)#rnorm(n1, mu_true, 1)
  Y2 <- -mu_true+rexp(n2)-1#rnorm(n,-mu_true,1)
  #Y1 <- rnorm(n1, mu_true, 1)
  #Y2 <- rnorm(n,-mu_true,1)
  for(k in 1:length(interval)){
    print(paste("Sample", l, "and second mean:", -mu_true+interval[k]))
    #print(k)
    Y2_new <- Y2 + interval[k]
    #FuRBI
    res <- sampler2(Y1,Y2_new, init, iterations, burnin)
    Z <- res$Z1
    p <- res$p
    numbers_Furbi <- res$numbers_K
    post <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post[j] <- mean(lik)
    }
    perf_ANRMI_matrix[l, k] <- sum(abs(post-true))
    #exchangeable
    res_univ2 <- sampler(c(Y1,Y2_new),init_univ, iterations, burnin)
    Z_univ2 <- res_univ2$Z
    p_univ2 <- res_univ2$p
    numbers_univ2 <- res_univ2$numbers_K
    post_univ2 <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ2[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ2[j] <- mean(lik)
    }
    
    perf_exch_matrix[l, k] <- sum(abs(post_univ2-true))
    
    #hierarchical
    res <- sampler_hier_trunc(Y1, Y2_new, init_hier, iterations,burnin)
    Z <- res$Z1
    p <- res$p1
    numbers_hier <- res$numbers_K
    post_hier <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_hier[j] <- mean(lik)
    }
    perf_hier_matrix[l, k] <- sum(abs(post_hier-true))
    
    #univariate
    res_univ <- sampler(c(Y1),init_univ, iterations, burnin)
    Z_univ <- res_univ$Z
    p_univ <- res_univ$p
    numbers_univ <- res_univ$numbers_K
    post_univ <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ[j] <- mean(lik)
    }
    
    perf_univ_matrix[l, k] <- sum(abs(post_univ-true))
    
    to_print <- c(perf_ANRMI_matrix[l,k], perf_exch_matrix[l,k],perf_univ_matrix[l,k], perf_hier_matrix[l,k])
    names(to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    numbers_to_print <- c(numbers_Furbi, numbers_univ2, numbers_univ, numbers_hier)
    names(numbers_to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    print(to_print)
    print(numbers_to_print)
    cat("")
  }
}
perf_ANRMI <- apply(perf_ANRMI_matrix, 2, median)
perf_exch <- apply(perf_exch_matrix, 2, median)
perf_univ <- apply(perf_univ_matrix, 2, median)
perf_hier <- apply(perf_hier_matrix, 2, median)
#extract errors
Errors <- matrix(c(-10+interval, 0.1*perf_exch, 0.1*perf_ind, 0.1*perf_ANRMI, 0.1*perf_hier), ncol = 5, byrow = F )
y_max <- max(Errors[,-1])
y_min <- min(Errors[,-1])
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="black", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="gray", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)
#saveRDS(to_use, file = "Errors_t_vs_exp.rds")

#4) Other simulations: t-student vs t-student

#GENERATE DATA
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
set.seed(2019)
mu_true = 10
n <- 100
n1 <- 20
n2 <- n
theta <- 1
Y1 <- mu_true + rt(n1, df = 3)#rnorm(n1, mu_true, 1)
Y2 <- -mu_true+rt(n2, df = 3)#rnorm(n,-mu_true,1)
interval <- seq(-6,26,by = 2)#c(-5, 0, 5, 10, 15, 20, 25)#-6:26
grid <- seq(-25,25, by = 0.1) #evaluation grid
burnin <- 3000
iterations <- burnin+1000 #iterations
N <- 150

#SELECT CONCENTRATION PARAMETER HDP
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

#true
true <- rep(0,length(grid))
for(i in 1:length(grid)){
  true[i] <- dt(grid[i]-mu_true,df = 3)
}

#INITIALIZER
#independent
#univariate sampler
init_univ <- list()
init_univ$mu0 <- 0
init_univ$tau0 <- 1
init_univ$sigma <- 1
init_univ$theta <- 1
init_univ$Z <- rnorm(N, mean = init_univ$mu0, sd = sqrt(init_univ$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init_univ$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_univ$p <- p
res_univ <- sampler(Y1,init_univ, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}
perf_ind <- rep(sum(abs(true-post_univ)),length(interval))

#bivariate initializer
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#hierarchical initializer
init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p

n_samples <- 20 #number of samples to average
perf_ANRMI_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_exch_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_hier_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_univ_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
for(l in 1:n_samples){
  #print(paste("Sample:", l))
  Y1 <- mu_true + rt(n1, df = 3)#rnorm(n1, mu_true, 1)
  Y2 <- -mu_true+rt(n2, df = 3)#rnorm(n,-mu_true,1)
  #Y1 <- rnorm(n1, mu_true, 1)
  #Y2 <- rnorm(n,-mu_true,1)
  for(k in 1:length(interval)){
    print(paste("Sample", l, "and second mean:", -mu_true+interval[k]))
    #print(k)
    Y2_new <- Y2 + interval[k]
    #FuRBI
    res <- sampler2(Y1,Y2_new, init, iterations, burnin)
    Z <- res$Z1
    p <- res$p
    numbers_Furbi <- res$numbers_K
    post <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post[j] <- mean(lik)
    }
    perf_ANRMI_matrix[l, k] <- sum(abs(post-true))
    #exchangeable
    res_univ2 <- sampler(c(Y1,Y2_new),init_univ, iterations, burnin)
    Z_univ2 <- res_univ2$Z
    p_univ2 <- res_univ2$p
    numbers_univ2 <- res_univ2$numbers_K
    post_univ2 <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ2[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ2[j] <- mean(lik)
    }
    
    perf_exch_matrix[l, k] <- sum(abs(post_univ2-true))
    
    #hierarchical
    res <- sampler_hier_trunc(Y1, Y2_new, init_hier, iterations,burnin)
    Z <- res$Z1
    p <- res$p1
    numbers_hier <- res$numbers_K
    post_hier <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_hier[j] <- mean(lik)
    }
    perf_hier_matrix[l, k] <- sum(abs(post_hier-true))
    
    #univariate
    res_univ <- sampler(c(Y1),init_univ, iterations, burnin)
    Z_univ <- res_univ$Z
    p_univ <- res_univ$p
    numbers_univ <- res_univ$numbers_K
    post_univ <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ[j] <- mean(lik)
    }
    
    perf_univ_matrix[l, k] <- sum(abs(post_univ-true))
    
    to_print <- c(perf_ANRMI_matrix[l,k], perf_exch_matrix[l,k],perf_univ_matrix[l,k], perf_hier_matrix[l,k])
    names(to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    numbers_to_print <- c(numbers_Furbi, numbers_univ2, numbers_univ, numbers_hier)
    names(numbers_to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    print(to_print)
    print(numbers_to_print)
    cat("")
  }
}
perf_ANRMI <- apply(perf_ANRMI_matrix, 2, median)
perf_exch <- apply(perf_exch_matrix, 2, median)
perf_univ <- apply(perf_univ_matrix, 2, median)
perf_hier <- apply(perf_hier_matrix, 2, median)
#extract errors
Errors <- matrix(c(-10+interval, 0.1*perf_exch, 0.1*perf_ind, 0.1*perf_ANRMI, 0.1*perf_hier), ncol = 5, byrow = F )
y_max <- max(Errors[,-1])
y_min <- min(Errors[,-1])
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="black", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="gray", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)
#saveRDS(to_use, file = "Errors_t_vs_t.rds")


#5) Other simulations: mixtures of two Gaussians

#GENERATE DATA
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
set.seed(2019)
mu_true = 10
n <- 100
n1 <- 50
n2 <- n
theta <- 1
generate_bivariate <- function(n, mu1, mu2, sigma, p){
  Y1 <- rnorm(n, mu1, sigma)
  Y2 <- rnorm(n, mu2, sigma)
  Y <- rep(0, n)
  U <- runif(n)
  Y[U <= p] <- Y1[U <= p]
  Y[U > p] <- Y2[U > p]
  return(Y)
}
Y1 <- mu_true + generate_bivariate(n1, -5, 5, 1, 0.5)
Y2 <- -mu_true+generate_bivariate(n2, -5, 5, 1, 0.5)
interval <- seq(-6,26,by = 2)#c(-5, 0, 5, 10, 15, 20, 25)#-6:26
grid <- seq(-25,25, by = 0.1) #evaluation grid
burnin <- 3000
iterations <- burnin+1000 #iterations
N <- 150

#SELECT CONCENTRATION PARAMETER HDP
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

#true
true <- rep(0,length(grid))
for(i in 1:length(grid)){
  true[i] <- 0.5*dnorm(grid[i]-mu_true, -5, 1)+0.5*dnorm(grid[i]-mu_true, 5, 1)
}

#INITIALIZER
#independent
#univariate sampler
init_univ <- list()
init_univ$mu0 <- 0
init_univ$tau0 <- 1
init_univ$sigma <- 1
init_univ$theta <- 1
init_univ$Z <- rnorm(N, mean = init_univ$mu0, sd = sqrt(init_univ$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init_univ$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_univ$p <- p
res_univ <- sampler(Y1,init_univ, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}
perf_ind <- rep(sum(abs(true-post_univ)),length(interval))

#bivariate initializer
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#hierarchical initializer
init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p

n_samples <- 50 #number of samples to average
perf_ANRMI_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_exch_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_hier_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_univ_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
for(l in 1:n_samples){
  #print(paste("Sample:", l))
  Y1 <- mu_true + generate_bivariate(n1, -5, 5, 1, 0.5)
  Y2 <- -mu_true+generate_bivariate(n2, -5, 5, 1, 0.5)
  #Y1 <- rnorm(n1, mu_true, 1)
  #Y2 <- rnorm(n,-mu_true,1)
  for(k in 1:length(interval)){
    print(paste("Sample", l, "and second mean:", -mu_true+interval[k]))
    #print(k)
    Y2_new <- Y2 + interval[k]
    #FuRBI
    res <- sampler2(Y1,Y2_new, init, iterations, burnin)
    Z <- res$Z1
    p <- res$p
    numbers_Furbi <- res$numbers_K
    post <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post[j] <- mean(lik)
    }
    perf_ANRMI_matrix[l, k] <- sum(abs(post-true))
    #exchangeable
    res_univ2 <- sampler(c(Y1,Y2_new),init_univ, iterations, burnin)
    Z_univ2 <- res_univ2$Z
    p_univ2 <- res_univ2$p
    numbers_univ2 <- res_univ2$numbers_K
    post_univ2 <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ2[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ2[j] <- mean(lik)
    }
    
    perf_exch_matrix[l, k] <- sum(abs(post_univ2-true))
    
    #hierarchical
    res <- sampler_hier_trunc(Y1, Y2_new, init_hier, iterations,burnin)
    Z <- res$Z1
    p <- res$p1
    numbers_hier <- res$numbers_K
    post_hier <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_hier[j] <- mean(lik)
    }
    perf_hier_matrix[l, k] <- sum(abs(post_hier-true))
    
    #univariate
    res_univ <- sampler(c(Y1),init_univ, iterations, burnin)
    Z_univ <- res_univ$Z
    p_univ <- res_univ$p
    numbers_univ <- res_univ$numbers_K
    post_univ <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ[j] <- mean(lik)
    }
    
    perf_univ_matrix[l, k] <- sum(abs(post_univ-true))
    
    to_print <- c(perf_ANRMI_matrix[l,k], perf_exch_matrix[l,k],perf_univ_matrix[l,k], perf_hier_matrix[l,k])
    names(to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    numbers_to_print <- c(numbers_Furbi, numbers_univ2, numbers_univ, numbers_hier)
    names(numbers_to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    print(to_print)
    print(numbers_to_print)
    cat("")
  }
}
#aux_ANRMI <- perf_ANRMI_matrix[1:9,]
#aux_exch <- perf_exch_matrix[1:9,]
#aux_univ <- perf_univ_matrix[1:9,]
#aux_hier <- perf_hier_matrix[1:9,]
perf_ANRMI <- apply(perf_ANRMI_matrix, 2, median)
perf_exch <- apply(perf_exch_matrix, 2, median)
perf_univ <- apply(perf_univ_matrix, 2, median)
perf_hier <- apply(perf_hier_matrix, 2, median)
#extract errors
Errors <- matrix(c(-10+interval, 0.1*perf_exch, 0.1*perf_ind, 0.1*perf_ANRMI, 0.1*perf_hier), ncol = 5, byrow = F )
y_max <- max(Errors[,-1])
y_min <- min(Errors[,-1])
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="black", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="gray", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)

#saveRDS(Errors, file = "Errors_biv_vs_biv.rds")

#6) Other simulations: exponential vs exponential

#GENERATE DATA
source("Functions.R")
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)
set.seed(2019)
mu_true = 10
n <- 100
n1 <- 20
n2 <- n
theta <- 1
Y1 <- mu_true + rexp(n1)-1#rnorm(n1, mu_true, 1)
Y2 <- -mu_true+rexp(n2)-1#rnorm(n,-mu_true,1)
interval <- seq(-6,26,by = 2)#c(-5, 0, 5, 10, 15, 20, 25)#-6:26
grid <- seq(-25,25, by = 0.1) #evaluation grid
burnin <- 3000
iterations <- burnin+1000 #iterations
N <- 150

#SELECT CONCENTRATION PARAMETER HDP
simulate_prob <- function(n, theta, n_simul = 100000){
  u <- matrix(runif(n_simul*n), ncol = n, nrow = n_simul)
  bol <- u <= matrix(theta/(0:(n-1)+theta), ncol = n, nrow = n_simul, byrow = T)
  val <- apply(bol, 1, sum)
  to_ret <- rep(0, n)
  res <- table(val)
  nomi <- as.numeric(names(res))
  to_ret[nomi] <- res
  return(to_ret/n_simul)
}
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2

n_cluster_hier <- function(prob1, prob2, theta0){
  m1 <- length(prob1)
  m2 <- length(prob2)
  res <- 0
  for(i in 1:m1){
    for(j in 1:m2){
      m <- i+j-2
      res <- res+ ifelse(m > 0, sum(theta0/(0:(m-1)+theta0))*prob1[i]*prob2[j], prob1[i]*prob2[j])
    }
  }
  return(res)
}
#expected number of clusters for DP
theta <- 1
prob1 <- simulate_prob(n1, theta)
prob2 <- simulate_prob(n2, theta)
mean1 <- sum((1:n1)*prob1)
mean2 <- sum((1:n2)*prob2)
mean_tot <- mean1+mean2
#We mimick with the HDP
theta_obs <- 5
prob1 <- simulate_prob(n1, theta_obs)
prob2 <- simulate_prob(n2, theta_obs)
theta0 = seq(1, 20, length = 300)
means <- c()
diff <- c()
for(val in theta0){
  value_mean <- n_cluster_hier(prob1, prob2, val)
  means <- c(means, value_mean)
  diff <- c(diff, abs(value_mean-mean_tot))
}
plot(theta0, means, type = "b", pch = 19)
abline(h = mean_tot, col = "red")
pos <- which(diff == min(diff))
theta0 <- theta0[pos]

#true
true <- rep(0,length(grid))
for(i in 1:length(grid)){
  true[i] <- dexp(grid[i]-mu_true+1)
}

#INITIALIZER
#independent
#univariate sampler
init_univ <- list()
init_univ$mu0 <- 0
init_univ$tau0 <- 1
init_univ$sigma <- 1
init_univ$theta <- 1
init_univ$Z <- rnorm(N, mean = init_univ$mu0, sd = sqrt(init_univ$tau0))
#p
V <- rbeta(N-1,shape1 = 1,init_univ$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_univ$p <- p
res_univ <- sampler(Y1,init_univ, iterations, burnin)
Z_univ <- res_univ$Z
p_univ <- res_univ$p
post_univ <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p_univ[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post_univ[j] <- mean(lik)
}
perf_ind <- rep(sum(abs(true-post_univ)),length(interval))

#bivariate initializer
init <-list()
init$mu1 = 0
init$mu2 = 0
init$tau1 = 1
init$tau2 = 1
init$sigma = 1
init$theta = 1
init$rho = 0
Sigma <- matrix(0, nrow = 2, ncol = 2)
diag(Sigma) <- c(init$tau1, init$tau2)
Sigma[2,1] <- init$rho*sqrt(init$tau1*init$tau2)
Sigma[1,2] <- Sigma[2,1]
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2), Sigma = Sigma)
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p

#hierarchical initializer
init_hier <- list()
init_hier$mu0 <- 0
init_hier$tau0 <- 1
init_hier$sigma <- 1
init_hier$theta0 <- theta0
init_hier$theta_obs <- theta_obs
init_hier$N <- N
init_hier$Z0 <- rnorm(N, mean = init_hier$mu0, sd = sqrt(init_hier$tau0))
#p0
V <- rbeta(N-1,shape1 = 1,init_hier$theta0)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p0 <- p
#p1
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p1 <- p
#p2
V <- rbeta(N-1,shape1 = 1,init_hier$theta_obs)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init_hier$p2 <- p

n_samples <- 50 #number of samples to average
perf_ANRMI_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_exch_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_hier_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
perf_univ_matrix <- matrix(0, nrow = n_samples, ncol = length(interval))
for(l in 1:n_samples){
  #print(paste("Sample:", l))
  Y1 <- mu_true + rexp(n1)-1#rnorm(n1, mu_true, 1)
  Y2 <- -mu_true+rexp(n2)-1#rnorm(n,-mu_true,1)
  #Y1 <- rnorm(n1, mu_true, 1)
  #Y2 <- rnorm(n,-mu_true,1)
  for(k in 1:length(interval)){
    print(paste("Sample", l, "and second mean:", -mu_true+interval[k]))
    #print(k)
    Y2_new <- Y2 + interval[k]
    #FuRBI
    res <- sampler2(Y1,Y2_new, init, iterations, burnin)
    Z <- res$Z1
    p <- res$p
    numbers_Furbi <- res$numbers_K
    post <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post[j] <- mean(lik)
    }
    perf_ANRMI_matrix[l, k] <- sum(abs(post-true))
    #exchangeable
    res_univ2 <- sampler(c(Y1,Y2_new),init_univ, iterations, burnin)
    Z_univ2 <- res_univ2$Z
    p_univ2 <- res_univ2$p
    numbers_univ2 <- res_univ2$numbers_K
    post_univ2 <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ2, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ2[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ2[j] <- mean(lik)
    }
    
    perf_exch_matrix[l, k] <- sum(abs(post_univ2-true))
    
    #hierarchical
    res <- sampler_hier_trunc(Y1, Y2_new, init_hier, iterations,burnin)
    Z <- res$Z1
    p <- res$p1
    numbers_hier <- res$numbers_K
    post_hier <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z, sd = sqrt(init_hier$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_hier[j] <- mean(lik)
    }
    perf_hier_matrix[l, k] <- sum(abs(post_hier-true))
    
    #univariate
    res_univ <- sampler(c(Y1),init_univ, iterations, burnin)
    Z_univ <- res_univ$Z
    p_univ <- res_univ$p
    numbers_univ <- res_univ$numbers_K
    post_univ <- rep(0, length(grid))
    for(j in 1:length(grid)){
      lik <- dnorm(grid[j], mean = Z_univ, sd = sqrt(init_univ$sigma))
      lik <- sapply(1:(iterations-burnin), function(i){
        return(sum(p_univ[i,]*lik[i,]))
      })
      #lik <- diag(p%*%t(lik))
      post_univ[j] <- mean(lik)
    }
    
    perf_univ_matrix[l, k] <- sum(abs(post_univ-true))
    
    to_print <- c(perf_ANRMI_matrix[l,k], perf_exch_matrix[l,k],perf_univ_matrix[l,k], perf_hier_matrix[l,k])
    names(to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    numbers_to_print <- c(numbers_Furbi, numbers_univ2, numbers_univ, numbers_hier)
    names(numbers_to_print) <- c("FuRBI", "Exch.","Ind.", "Hier.")
    print(to_print)
    print(numbers_to_print)
    cat("")
  }
}
perf_ANRMI <- apply(perf_ANRMI_matrix, 2, median)
perf_exch <- apply(perf_exch_matrix, 2, median)
perf_univ <- apply(perf_univ_matrix, 2, median)
perf_hier <- apply(perf_hier_matrix, 2, median)
#extract errors
Errors <- matrix(c(-10+interval, 0.1*perf_exch, 0.1*perf_ind, 0.1*perf_ANRMI, 0.1*perf_hier), ncol = 5, byrow = F )
y_max <- max(Errors[,-1])
y_min <- min(Errors[,-1])
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="black", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="gray", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)
#saveRDS(to_use, file = "Errors_exp_vs_exp.rds")

plot(Errors[,1],Errors[,2],col="white", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 1, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
legend(-15,1.7, legend = c( "Exchangeable", "Independent","FuRBI", "Hierarchical"), col = c("gray", "brown", "blue", "red"), lwd = 4, pch = c(15,19,13,16),lty = c(9,10,1,1), cex = 0.9, bty = "n")


#making plots
y_max <- max(Errors[,-1])
y_min <- min(Errors[,-1])
#plot
plot(Errors[,1],Errors[,2],col="gray", xlab = "Mean second group",ylab = "Error",main = "",pch = 15,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(Errors[,1],Errors[,3],col="brown", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1], Errors[,4],col="blue", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(Errors[,1],Errors[,5],col="red", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-16,16),ylim = c(y_min, y_max),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 8,las = 1,cex.lab = 1.4, ann = F, axes = F)
