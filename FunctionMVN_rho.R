#samplers
rho_sampler <- function(rho, Z1, Z2, mu1, mu2, tau1, tau2){
  #sample rho
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
K_sampler <- function(Y,Z,p, sigma){
  #sample labels
  #Y = n data points
  #Z = N atoms
  #p = N weights
  n <- length(Y)
  N <- length(Z)
  K <- rep(0, n)
  for( i in 1:n){
    lik <- dnorm(Y[i],mean = Z, sd = sigma,log = T)
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
Znull <- function(m,mu1,mu2,mu3, W){
  #sample atoms without connected data - bivariate case
  res <- mvrnorm(m,c(mu1, mu2, mu3), W)
  return(res)
}
Zused <- function(Pos,Y,K,Z2,Z3, mu, mu1, mu2, tau,Sigma12,Sigma22, sigma, label){
  #sample data with connected data - bivariate case
  N <- length(Z2)
  Z <- rep(0, N)
  #multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  #print(dim(Sigma12))
  #print(dim(Sigma22))
  Sigma_inv <-Sigma12%*%solve(Sigma22)
  #print(dim(Sigma_inv))
  #print(dim(t(Sigma12)))
  tau_bar <- ifelse(tau^2-Sigma_inv%*%t(Sigma12) > 0, tau^2-Sigma_inv%*%t(Sigma12), 1e-6)
  mu_vec <- matrix(c(mu1, mu2), ncol = 1)
  for(i in Pos){
    m <- M[i]
    y_mean <- ifelse(M[i] > 0,mean(Y[K == i]),0)
    Z_vec <- matrix(c(Z2[i], Z3[i]), ncol = 1)
    mu_bar <- mu+Sigma_inv%*%(Z_vec-mu_vec)
    new_mu <- (m*tau_bar*y_mean+mu_bar*sigma^2)/(sigma^2+m*tau_bar)
    new_tau <- sigma^2*tau_bar/(sigma^2+m*tau_bar)
    #print(paste("Label", label))
    #print(paste("Y_mean", y_mean))
    #print(paste("Tau", tau))
    #print(paste("Sigma_inv", Sigma_inv%*%t(Sigma12)))
    #print(paste("Tau_bar", tau_bar))
    #print(paste("New_Mu", new_mu))
    #print(paste("New_Tau", sqrt(new_tau)))
    #print(" ")
    Z[i] <- rnorm(1, mean = new_mu, sd = sqrt(new_tau))
  }
  #print(paste("Label", label))
  #print(summary(Z[pos]))
  return(Z)
}
sampler2 <- function(Y1,Y2,Y3, initializer, iterations, burnin){
  #sample bivariate case
  mu1 <- initializer$mu1
  mu2 <- initializer$mu2
  mu3 <- initializer$mu3
  tau1 <- initializer$tau1
  tau2 <- initializer$tau2
  tau3 <- initializer$tau3
  sigma <- initializer$sigma
  theta <- initializer$theta
  rho12 <- initializer$rho12
  rho13 <- initializer$rho13
  rho23 <- initializer$rho23
  Z1 <- initializer$Z1
  Z2 <- initializer$Z2
  Z3 <- initializer$Z3
  p <- initializer$p
  N <- length(Z1) #number of components
  
  Z1store <- matrix(0, ncol = N, nrow = iterations-burnin)
  Z2store <- matrix(0, ncol = N, nrow = iterations-burnin)
  Z3store <- matrix(0, ncol = N, nrow = iterations-burnin)
  pstore <- matrix(0, ncol = N, nrow = iterations-burnin)
  rho12store <- rep(0, iterations-burnin)
  rho13store <- rep(0, iterations-burnin)
  rho23store <- rep(0, iterations-burnin)
  
  for(i in 1:burnin){
    if(i %% 500 == 0){
      print(paste("Iteration: ", i))
    }
    #sample K
    K1 <- K_sampler(Y1,Z1,p, sigma)
    K2 <- K_sampler(Y2,Z2,p, sigma)
    K3 <- K_sampler(Y3,Z3,p, sigma)
    K <- c(K1,K2, K3)
    #Find multiplicities
    res <- table(K) #summary of labels
    pos <- as.numeric(names(res)) #atoms with at least one label
    #sample p
    p <- p_sampler(K, theta,N)
    #sample Z used
    W <- matrix(c(tau1^2,rho12*tau1*tau2, rho13*tau1*tau3,rho12*tau1*tau2,tau2^2, rho23*tau2*tau3,rho13*tau1*tau3, rho23*tau2*tau3, tau3^2), nrow = 3, byrow = T)
    V <- eigen(W)$vectors
    d <- eigen(W)$values
    d[d <= 0] <- 1e-6
    W <- V%*%diag(d)%*%t(V)
    Z1 <- Zused(pos,Y1,K1,Z2,Z3, mu1, mu2, mu3,tau1^2, t(as.matrix(W[1,c(2,3)])),W[2:3, 2:3], sigma, 1)
    Z2 <- Zused(pos,Y2,K2,Z1,Z3, mu2, mu1, mu3, tau2^2, t(as.matrix(W[2,c(1,3)])),W[c(1,3), c(1,3)], sigma, 2)
    Z3 <- Zused(pos,Y3,K3,Z1,Z2, mu3, mu1, mu2, tau3^3, t(as.matrix(W[3,c(1,2)])),W[c(1,2), c(1,2)], sigma, 3)
    #sample rho
    rho12 <- rho_sampler(rho12, Z1[pos], Z2[pos], mu1, mu2, tau1, tau2)
    rho13 <- rho_sampler(rho13, Z1[pos], Z3[pos], mu1, mu3, tau1, tau3)
    rho23 <- rho_sampler(rho23, Z2[pos], Z3[pos], mu2, mu3, tau2, tau3)
    #sample Z unused
    aux <- Znull(N-length(pos),mu1,mu2,mu3, W)
    Z1[-pos] <- aux[,1]
    Z2[-pos] <- aux[,2]
    Z3[-pos] <- aux[,3]
  }
  print("End burnin")
  
  for(i in 1:(iterations-burnin)){
    if(i %% 500 == 0){
      print(paste("Iteration: ", i))
    }
    #sample K
    K1 <- K_sampler(Y1,Z1,p, sigma)
    K2 <- K_sampler(Y2,Z2,p, sigma)
    K3 <- K_sampler(Y3,Z3,p, sigma)
    K <- c(K1,K2, K3)
    #Find multiplicities
    res <- table(K) #summary of labels
    pos <- as.numeric(names(res)) #atoms with at least one label
    #sample p
    p <- p_sampler(K, theta,N)
    pstore[i,] <- p
    #sample Z used
    W <- matrix(c(tau1^2,rho12*tau1*tau2, rho13*tau1*tau3,rho12*tau1*tau2,tau2^2, rho23*tau2*tau3,rho13*tau1*tau3, rho23*tau2*tau3, tau3^2), nrow = 3, byrow = T)
    V <- eigen(W)$vectors
    d <- eigen(W)$values
    d[d <= 0] <- 1e-6
    W <- V%*%diag(d)%*%t(V)
    Z1 <- Zused(pos,Y1,K1,Z2,Z3, mu1, mu2, mu3,tau1^2, t(as.matrix(W[1,c(2,3)])),W[2:3, 2:3], sigma, 1)
    Z2 <- Zused(pos,Y2,K2,Z1,Z3, mu2, mu1, mu3, tau2^2, t(as.matrix(W[2,c(1,3)])),W[c(1,3), c(1,3)], sigma, 2)
    Z3 <- Zused(pos,Y3,K3,Z1,Z2, mu3, mu1, mu2, tau3^3, t(as.matrix(W[3,c(1,2)])),W[c(1,2), c(1,2)], sigma, 3)
    #sample rho
    rho12 <- rho_sampler(rho12, Z1[pos], Z2[pos], mu1, mu2, tau1, tau2)
    rho13 <- rho_sampler(rho13, Z1[pos], Z3[pos], mu1, mu3, tau1, tau3)
    rho23 <- rho_sampler(rho23, Z2[pos], Z3[pos], mu2, mu3, tau2, tau3)
    rho12store[i] <- rho12
    rho13store[i] <- rho13
    rho23store[i] <- rho23
    #sample Z unused
    aux <- Znull(N-length(pos),mu1,mu2,mu3, W)
    Z1[-pos] <- aux[,1]
    Z2[-pos] <- aux[,2]
    Z3[-pos] <- aux[,3]
    Z1store[i,] <- Z1
    Z2store[i,] <- Z2
    Z3store[i,] <- Z3
  }
  #return objects
  results <- list()
  results$p <- pstore
  results$Z1 <- Z1store
  results$Z2 <- Z2store
  results$Z3 <- Z3store
  results$rho12 <- rho12store
  results$rho13 <- rho13store
  results$rho23 <- rho23store
  return(results)
}