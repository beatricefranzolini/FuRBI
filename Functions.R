#FuRBI priors
rbivnorm <- function(n,mu1,mu2,tau1,tau2,rho){
  #samples from the bivariate normal with correlation rho
  Sigma <- matrix(0, nrow = 2, ncol = 2)
  diag(Sigma) <- c(tau1, tau2)
  Sigma[2,1] <- rho*sqrt(tau1*tau2)
  Sigma[1,2] <- rho*sqrt(tau1*tau2)
  Mu <-c(mu1, mu2)
  res <- mvrnorm(n, mu = Mu, Sigma = Sigma)
  return(res)
}
#samplers
K_sampler <- function(Y,Z,p){
  #sample labels
  #Y = n data points
  #Z = N atoms
  #p = N weights
  n <- length(Y)
  N <- length(Z)
  K <- rep(0, n)
  for( i in 1:n){
    lik <- dnorm(Y[i],mean = Z, sd = 1,log = T)
    w <- log(p)+lik
    w <- exp(w-logSumExp(w))
    K[i] <- sample(1:N,size = 1, prob = w)
  }
  return(K)
}
p_sampler <- function(K,theta,N){
  #sample weights
  #K = n classifications labels
  #theta = concentration parameter
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  M <- M[-N]
  #sample betas
  V <- rbeta(length(M),shape1 = 1+M,theta+sum(M)-cumsum(M))
  W <- cumprod(1-V)
  V <- c(V,1)
  W <- c(1,W)
  p <- V*W
  return(p)
}
Z_sampler <- function(Y,K,mu0, sigma, tau0, N){
  #sample atoms - univariate case
  #Y = n data points
  #K = n classifications labels
  #theta = concentration parameter
  #mu0 = prior mean of the mean
  #sigma = variance of the data
  #tau0 = prior variance of the mean
  Z <- rep(0, N)
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  Z[M == 0] = rnorm(sum(M == 0),mean = mu0, sd = sqrt(tau0))
  for(i in pos){
    m <- M[i]
    y_mean <- mean(Y[K == i])
    mu <- (mu0/tau0+m*y_mean/sigma)/(1/tau0+m/sigma)
    tau <- 1/(1/tau0+m/sigma)
    Z[i] <- rnorm(1, mean = mu, sd = sqrt(tau))
  }
  return(Z)
}
Z_sampler2 <- function(Y, K, Z2, mu1, tau1, mu2, tau2, sigma, rho){
  #sample atoms bivariate - NOT USED
  #Y = n data points from Z1
  #K = n classifications labels of Y
  #Z2 = N atoms from X2
  #mu1 = prior mean of Z1
  #tau1 = prior variance of Z1
  #mu2 = prior mean of Z2
  #tau2 = prior variance of Z2
  #sigma = variance of the data
  #rho = correlation between Z1 and Z2
  N <- length(Z2)
  Z <- rep(0, N)
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  tau3 <- tau2*(1-rho^2)
  d <- rho*sqrt(tau2/tau1)
  for(i in 1:N){
    m <- M[i]
    y_mean <- ifelse(m > 0,mean(Y[K == i]),0)
    c <- Z2[i]-mu2
    mu <- (mu1*sigma*tau3+m*y_mean*tau1*tau3+d*(c+d)*tau1*sigma)/(sigma*tau3+m*tau1*tau3+d^2*tau1*sigma)
    if(mu < 6 & mu > 4){
      print("Mu:")
      print(mu)
      print("Y")
      print(y_mean)
      print("c")
      print(c)
    }
    tau <- 1/(1/tau1+m/sigma+d^2/tau3)
    Z[i] <- rnorm(1, mean = mu, sd = sqrt(tau))
  }
  return(Z)
}
Znull <- function(m,mu1,mu2,tau1,tau2,rho){
  #sample atoms without connected data - bivariate case
  res <- rbivnorm(m,mu1,mu2,tau1,tau2,rho)
  return(res)
}
Zused <- function(Pos, Y, K, Z2, mu1, tau1, mu2, tau2, sigma, rho){
  #sample data with connected data - bivariate case
  N <- length(Z2)
  Z <- rep(0, N)
  #multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  tau3 <- tau2*(1-rho^2)
  d <- rho*sqrt(tau2/tau1)
  for(i in Pos){
    m <- M[i]
    y_mean <- ifelse(M[i] > 0,mean(Y[K == i]),0)
    c <- Z2[i]-mu2
    mu <- (mu1*sigma*tau3+m*y_mean*tau1*tau3+d*(c+d)*tau1*sigma)/(sigma*tau3+m*tau1*tau3+d^2*tau1*sigma)
    tau <- 1/(1/tau1+m/sigma+d^2/tau3)
    Z[i] <- rnorm(1, mean = mu, sd = sqrt(tau))
  }
  return(Z)
}
rho_sampler <- function(rho, Z1, Z2, mu1, mu2, tau1, tau2){
  #sample rho using MH
  eps  <- 0.05
  rho_new <- runif(1,rho-eps,rho+eps)
  if(rho_new < -1 | rho_new >1){
    return(rho)
  }
  oldval <- -0.5*log(1-rho^2)-0.5/(1-rho^2)*sum((Z1-mu1)^2/tau1+(Z2-mu2)^2/tau2-2*rho*(Z1-mu1)/sqrt(tau1)*(Z2-mu2)/sqrt(tau2))
  newval <- -0.5*log(1-rho_new^2)-0.5/(1-rho_new^2)*sum((Z1-mu1)^2/tau1+(Z2-mu2)^2/tau2-2*rho_new*(Z1-mu1)/sqrt(tau1)*(Z2-mu2)/sqrt(tau2))
  logratio <- newval-oldval
  rho <- ifelse(log(runif(1)) < logratio, rho_new, rho)
  return(rho)
}
sampler <- function(Y, initializer, iterations,burnin, verbose = F){
  #sampler univariate case
  mu0 <- initializer$mu0
  sigma <- initializer$sigma
  tau0 <- initializer$tau0
  theta <- initializer$theta
  Z <- initializer$Z
  p <- initializer$p
  N <- length(Z) #number of components
  
  Zstore <- matrix(0, ncol = N, nrow = iterations-burnin)
  pstore <- matrix(0, ncol = N, nrow = iterations-burnin)
  
  for(i in 1:burnin){
    K <- K_sampler(Y,Z,p)
    p <- p_sampler(K, theta,N)
    Z <- Z_sampler(Y,K, mu0, sigma, tau0, N)
  }
  if(verbose){
    print("End burnin")
  }
  
  numbers_K <- 0
  for(i in 1:(iterations-burnin)){
    if(i %% 500 == 0 & verbose){
      print(paste("Iteration: ", i))
    }
    #sample K
    K <- K_sampler(Y,Z,p)
    numbers_K <- numbers_K + length(unique(K))
    #sample p
    p <- p_sampler(K, theta,N)
    pstore[i,] <- p
    #sample Z
    Z <- Z_sampler(Y,K, mu0, sigma, tau0, N)
    Zstore[i,] <- Z
  }
  #return objects
  results <- list()
  results$p <- pstore
  results$Z <- Zstore
  results$numbers_K <- numbers_K/(iterations-burnin)
  return(results)
}
sampler2 <- function(Y1,Y2, initializer, iterations, burnin, verbose = F){
  #sample bivariate case
  mu1 <- initializer$mu1
  tau1 <- initializer$tau1
  mu2 <- initializer$mu2
  tau2 <- initializer$tau2
  rho <- initializer$rho
  sigma <- initializer$sigma
  theta <- initializer$theta
  Z1 <- initializer$Z1
  Z2 <- initializer$Z2
  p <- initializer$p
  N <- length(Z1) #number of components
  
  Z1store <- matrix(0, ncol = N, nrow = iterations-burnin)
  pstore <- matrix(0, ncol = N, nrow = iterations-burnin)
  Z2store <- matrix(0, ncol = N, nrow = iterations-burnin)
  Rho <- rep(0, iterations-burnin)
  
  
  for(i in 1:burnin){
    #sample K
    K1 <- K_sampler(Y1,Z1,p)
    K2 <- K_sampler(Y2,Z2,p)
    K <- c(K1,K2)
    #Find multiplicities
    res <- table(K) #summary of labels
    pos <- as.numeric(names(res)) #atoms with at least one label
    #sample p
    p <- p_sampler(K, theta,N)
    #sample Z used
    Z1 <- Zused(pos,Y1,K1,Z2,mu1, tau1, mu2, tau2, sigma, rho)
    Z2 <- Zused(pos,Y2,K2,Z1, mu2, tau2, mu1, tau1, sigma, rho)
    #sample rho
    rho <- rho_sampler(rho, Z1[pos], Z2[pos], mu1, mu2, tau1, tau2)
    #rho = -0.9
    #sample Z unused
    aux <- Znull(N-length(pos),mu1,mu2,tau1,tau2,rho)
    Z1[-pos] <- aux[,1]
    Z2[-pos] <- aux[,2]
  }
  if(verbose){
    print("End burnin")
  }
  
  numbers_K <- 0
  for(i in 1:(iterations-burnin)){
    if(i %% 500 == 0 & verbose){
      print(paste("Iteration: ", i))
    }
    #sample K
    K1 <- K_sampler(Y1,Z1,p)
    K2 <- K_sampler(Y2,Z2,p)
    K <- c(K1,K2)
    numbers_K <- numbers_K+length(unique(K))
    #Find multiplicities
    res <- table(K) #summary of labels
    pos <- as.numeric(names(res)) #atoms with at least one label
    #sample p
    p <- p_sampler(K, theta,N)
    pstore[i,] <- p
    #sample Z used
    Z1 <- Zused(pos,Y1,K1,Z2,mu1, tau1, mu2, tau2, sigma, rho)
    Z2 <- Zused(pos,Y2,K2,Z1, mu2, tau2, mu1, tau1, sigma, rho)
    #sample rho
    rho <- rho_sampler(rho, Z1[pos], Z2[pos], mu1, mu2, tau1, tau2)
    #rho = -0.9
    Rho[i] <- rho
    #sample Z unused
    aux <- Znull(N-length(pos),mu1,mu2,tau1,tau2,rho)
    Z1[-pos] <- aux[,1]
    Z2[-pos] <- aux[,2]
    Z1store[i,] <- Z1
    Z2store[i,] <- Z2
  }
  #return objects
  results <- list()
  results$p <- pstore
  results$Z1 <- Z1store
  results$Z2 <- Z2store
  results$rho <- Rho
  results$numbers_K <- numbers_K/(iterations-burnin)
  return(results)
}


sampler_hier_trunc <- function(Y1, Y2, initializer, iterations,burnin, verbose = F){
  #sampler hierarchical case
  mu0 <- initializer$mu0
  sigma <- initializer$sigma
  tau0 <- initializer$tau0
  theta_obs <- initializer$theta_obs
  theta0 <- initializer$theta0
  Z0 <- initializer$Z0 #at level 0
  Z1 <- duplicate(Z0)
  Z2 <- duplicate(Z0)
  p0 <- initializer$p0
  p1 <- initializer$p1
  p2 <- initializer$p2
  N <- initializer$N #number of components
  
  Z1store <- matrix(0, ncol = N, nrow = iterations-burnin)
  Z2store <- matrix(0, ncol = N, nrow = iterations-burnin)
  p1store <- matrix(0, ncol = N, nrow = iterations-burnin)
  p2store <- matrix(0, ncol = N, nrow = iterations-burnin)
  
  for(i in 1:burnin){
    K1 <- K_sampler(Y1,Z1,p1)
    K2 <- K_sampler(Y2,Z2,p2)
    p1 <- p_sampler(K1, theta_obs,N)
    p2 <- p_sampler(K2, theta_obs,N)
    T1 <- T_sampler_hier(Y1,K1, Z0,p0)
    Z1 <- Z0[T1]
    T2 <- T_sampler_hier(Y2,K2, Z0,p0)
    Z2 <- Z0[T2]
    #update root
    res <- table(K1) #summary of labels
    pos1 <- as.numeric(names(res)) #atoms with at least one label
    res <- table(K2) #summary of labels
    pos2 <- as.numeric(names(res)) #atoms with at least one label
    p0 <- p_sampler(c(T1[pos1],T2[pos2]), theta0,N)
    Z0 <- Z_sampler(c(Y1, Y2),c(T1[K1],T2[K2]), mu0, sigma, tau0, N)
  }
  if(verbose){
    print("End burnin")
  }
  
  numbers_K <- 0
  for(i in 1:(iterations-burnin)){
    if(i %% 500 == 0 & verbose){
      print(paste("Iteration: ", i))
    }
    K1 <- K_sampler(Y1,Z1,p1)
    K2 <- K_sampler(Y2,Z2,p2)
    p1 <- p_sampler(K1, theta_obs,N)
    p2 <- p_sampler(K2, theta_obs,N)
    T1 <- T_sampler_hier(Y1,K1, Z0,p1)
    Z1 <- Z0[T1]
    T2 <- T_sampler_hier(Y2,K2, Z0,p2)
    Z2 <- Z0[T2]
    #update root
    res <- table(K1) #summary of labels
    pos1 <- as.numeric(names(res)) #atoms with at least one label
    res <- table(K2) #summary of labels
    pos2 <- as.numeric(names(res)) #atoms with at least one label
    p0 <- p_sampler(c(T1[pos1],T2[pos2]), theta0,N)
    Z0 <- Z_sampler(c(Y1, Y2),c(T1[K1],T2[K2]), mu0, sigma, tau0, N)
    numbers_K <- numbers_K+length(unique(c(T1[K1],T2[K2])))
    #stare values
    Z1store[i,] <- Z1
    Z2store[i,] <- Z2
    p1store[i, ] <- p1
    p2store[i, ] <- p2
  }
  #return objects
  results <- list()
  results$p1 <- p1store
  results$p2 <- p2store
  results$Z1 <- Z1store
  results$Z2 <- Z2store
  results$numbers_K <- numbers_K/(iterations-burnin)
  return(results)
}
T_sampler_hier <- function(Y,K, Z,p){
  #sample labels
  #Y = n data points
  #K = n labels from observations to first level
  #Z = N atoms of the latent layer
  #p = N weights of the latent layer
  N <- length(Z)
  TT <- rep(0, N)
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  TT[-pos] <- sample(1:N,size = length(TT[-pos]), prob = p, replace = T)
  for( i in pos){
    lik <- dnorm(mean(Y[K == i]),mean = Z, sd = 1,log = T)
    #lik <- dnorm(mean(Y[K == i]),mean = Z, sd = 1,log = T)
    w <- log(p)+lik
    w <- exp(w-logSumExp(w))
    TT[i] <- sample(1:N,size = 1, prob = w)
  }
  return(TT)
}