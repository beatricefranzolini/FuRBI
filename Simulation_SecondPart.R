#Inverse Wishart
source("FunctionMVN_rho.R")
library(LaplacesDemon)
library(matrixStats)
library(MASS)
library(bivariate)
library(rlang)

#Generate bivariate data 
set.seed(2019)
#Independent case
n <- 100
mu1_true <- 10
mu2_true <- -10
mu3_true <- 10
tau_true <- 1
Y <- matrix(0, ncol = 3, nrow = n)
Y[,1] <- rnorm(n,mu1_true,tau_true)
Y[,2] <- rnorm(n,mu2_true,tau_true)
Y[,3] <- rnorm(n,mu3_true,tau_true)
Y1 <- Y[,1]
Y2 <- Y[,2]
Y3 <- Y[,3]
grid <- seq(-15,15, by = 0.1) #evaluation grid

theta <- 1

N <- 150
init <-list()
init$mu1 = 0
init$mu2 = 0
init$mu3 = 0
init$sigma = 1
init$theta = theta
init$tau1 <- 1
init$tau2 <- 1
init$tau3 <- 1
init$rho12 <- 0
init$rho13 <- 0
init$rho23 <- 0
init$nu <- 10
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2, init$mu3), Sigma = diag(3))
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]
init$Z3 <- Z[,3]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p
burnin <- 3000
iterations <- burnin+1000

res <- sampler2(Y1,Y2,Y3, init, iterations, burnin)

p <- res$p
Z1 <- res$Z1
post1 <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z1, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post1[j] <- mean(lik)
}
Z2 <- res$Z2
post2 <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z2, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post2[j] <- mean(lik)
}
Z3 <- res$Z3
post3 <- rep(0, length(grid))
for(j in 1:length(grid)){
  lik <- dnorm(grid[j], mean = Z3, sd = sqrt(init$sigma))
  lik <- sapply(1:(iterations-burnin), function(i){
    return(sum(p[i,]*lik[i,]))
  })
  #lik <- diag(p%*%t(lik))
  post3[j] <- mean(lik)
}
true <- rep(0, length(grid))
for(i in 1:length(grid)){
  true[i] <- mean(dnorm(grid[i], mean = mu1_true, sd = sqrt(tau_true)))
}
plot(p[iterations-burnin,], type = "p", pch = 19)
plot(grid, post1, col = "red", type = "l", lwd = 2, ylim = c(0,0.5))
lines(grid, true, lwd = 2)
lines(grid, post3, lwd = 2, col = "green")
lines(grid, post2, lwd = 2, col = "blue")

rho13 <- res$rho12
rho12 <- res$rho13
rho23 <- res$rho23

plot(rho12, type = "l")
plot(rho13, type = "l")
plot(rho23, type = "l")

#extensive simulation
source("FunctionMVN_rho.R")
set.seed(1920)
interval <- -10:10
theta <- 1

N <- 150
init <-list()
init$mu1 = 0
init$mu2 = 0
init$mu3 = 0
init$sigma = 1
init$theta = theta
init$tau1 <- 1
init$tau2 <- 1
init$tau3 <- 1
init$rho12 <- 0
init$rho13 <- 0
init$rho23 <- 0
init$nu <- 10
Z <- mvrnorm(N, mu = c(init$mu1, init$mu2, init$mu3), Sigma = diag(3))
init$Z1 <- Z[,1]
init$Z2 <- Z[,2]
init$Z3 <- Z[,3]

#p
V <- rbeta(N-1,shape1 = 1,init$theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W
init$p <- p
burnin <- 10000
iterations <- burnin+1000

#generate data
n <- 20
mu1_true <- 10
mu2_true <- -10
mu3_true <- 0
tau_true <- 1
Y <- matrix(0, ncol = 3, nrow = n)
Y[,1] <- rnorm(n,mu1_true,tau_true)
Y[,2] <- rnorm(n,mu2_true,tau_true)
Y[,3] <- rnorm(n,mu3_true,tau_true)
Y1 <- Y[,1]
Y2 <- Y[,2]
Y3 <- Y[,3]
grid <- seq(-15,15, by = 0.1) #evaluation grid

n_samples <- 100
Rho12 <- matrix(0,nrow = n_samples, ncol = length(interval))
Rho13 <- matrix(0,nrow = n_samples, ncol = length(interval))
Rho23 <- matrix(0,nrow = n_samples, ncol = length(interval))
for(k in 1:length(interval)){
  print(paste("Mean:", interval[k]))
  Y3_new <- Y3+interval[k]
  for(l in 1:n_samples){
  res <- sampler2(Y1,Y2,Y3_new, init, iterations, burnin)
  rho12 <- res$rho12
  rho13 <- res$rho13
  rho23 <- res$rho23
  Rho12[l,k] <- median(rho12)
  Rho13[l,k] <- median(rho13)
  Rho23[l,k] <- median(rho23)
  }
  #to_print <- c(median_rho12[k], median_rho13[k],median_rho23[k])
  #names(to_print) <- c("Rho12", "Rho13", "Rho23")
  #print(to_print)
}
#saveRDS(Rho12, file = "Rho12.rds")
#saveRDS(Rho13, file = "Rho13.rds")
#saveRDS(Rho23, file = "Rho23.rds")
interval <- -10:10
median_rho12 <- rep(0, length(interval))
median_rho13 <- rep(0, length(interval))
median_rho23 <- rep(0, length(interval))

for(k in 1:length(interval)){
  median_rho12[k] <- median(Rho12[,k])
  median_rho13[k] <- median(Rho13[,k])
  median_rho23[k] <- median(Rho23[,k])
}


plot(interval,median_rho12,col="black", xlab = "Mean of the second component",ylab = "Correlation",main = "",pch = 17,type = "b", xlim = c(-10,10),ylim = c(-1, 1),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(interval,median_rho13,col="red", xlab = "Support",ylab = "Density",main = "",pch = 15,type = "b", xlim = c(-10,10),ylim = c(-1, 1),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(interval, median_rho23,col="green", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-10,10),ylim = c(-1, 1),cex.axis = 1.6, cex = 0.4, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
#legend(-10,0.5, legend = c("Rho12", "Rho13","Rho23"), col = c(rep("black",1),"red", "green"), lwd = 4, pch = c(17,15,19),lty = c(2,9,10), cex = 1.4, bty = "n")
