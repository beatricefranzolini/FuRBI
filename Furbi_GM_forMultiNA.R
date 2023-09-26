#Furbi_GM_forMultiNA
#missing data
library(progress)#to draw the progress bar
library(hypergeo) #for hypergeometric function
library(mvtnorm) #for multivariate normal
library(fossil) #to compute rand indexes
library(mice) #to impute missing data
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(matrixcalc) #check if positive def
library(visdat) #to plot missing
library(LaplacesDemon) #to sample invwishart
library(Matrix)#to ensure pos def of covariance
library(naniar)#to plot missing
library(ggplot2)
library(coda)

sim_data <- function(scenario = 1, seed = 0){
  set.seed(seed)
  p = 3 #number of variables
  #correlation matrix 
  rho = diag(p) 
  off_diag = c(-0.9, 0.9, -0.9) 
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  rho[col<row] = off_diag
  rho = rho + t(rho)
  
  K_true = 4 #true number of clusters
  theta_true = rmvnorm(K_true, rep(0,3), rho)
  
  nc_true = rep(250,K_true) #true clusters numerosities
  n = sum(nc_true)
  data = NULL; c_true = NULL
  for (c in 1:K_true){
    data = rbind(data,rmvnorm(nc_true[c], theta_true[c,], diag(p)))
    c_true = c(c_true, rep(c, nc_true[c]))
  }
  data = as.data.frame(data)
  data_true = scale(data, scale = FALSE)
  
  if(scenario == 1){
    set.seed(seed)
    data = data_true
    miss1 = floor(runif(floor(n/5))*n); miss2 = floor(runif(floor(n/5))*n);
    miss3 = floor(runif(floor(n/5))*n); miss3 = setdiff(setdiff(miss3,miss1), miss2)
    data[miss1, 1] = NA
    data[miss2, 2] = NA
    data[miss3, 3] = NA
    
  }else if(scenario == 2){
    set.seed(seed)
    data = data_true
    data[c_true == 1, 1] = NA
    data[c_true == 3, 3] = NA
    
    
  }else if(scenario == 3){
    set.seed(seed)
    data = data_true
    data[floor(runif(n/3)*n), 1] = NA
    data[floor(runif(n/3)*n), 2] = NA
    data[seq(1,n)[rowSums(is.na(data)) == 0], 3] = NA
    
  }else if(scenario == 4){
    set.seed(seed)
    data = data_true
    data[sample(which(c_true == 1), floor(n/5) ), 1] = NA
    data[sample(which(c_true == 2), floor(n/5) ), 2] = NA
    data[sample(which(c_true %in% c(3,4) ), floor(n/10) ), 1] = NA
    data[sample(which(c_true %in% c(3,4) ), floor(n/10) ), 2] = NA
    data[seq(1,n)[rowSums(is.na(data)) == 0], 3] = NA
    
  }else if(scenario == 5){
    set.seed(seed)
    data = data_true
    data[sample(which(c_true == 1), floor(n/6) ), 1] = NA
    data[sample(which(c_true == 3), floor(n/6) ), 3] = NA
    data[sample(which(c_true %in% c(2,3,4)), floor(n/5) ), 2] = NA
    
  }else if(scenario == 6){
    set.seed(seed)
    data = data_true
    data[(which(c_true == 1) ), 1] = NA
    data[(which(c_true == 2) ), 2] = NA
    data[sample(which(c_true %in% c(3,4) ), floor(n/20) ), 1] = NA
    data[sample(which(c_true %in% c(3,4) ), floor(n/20) ), 2] = NA
    data[seq(1,n)[rowSums(is.na(data)) == 0], 3] = NA
  }
  return(list("data" = data, "data_true" = data_true, "c_true" = c_true, 
              "theta_true" = theta_true))
}

avg_sil <- function(k) {
  km.res <- kmeans(df, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}


#compute the marginal likelihood of a norm-norm model (marg lik of a cluster) univariate obs
margnn <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  sigma = sqrt(sigma2); s = sqrt(s2)
  n = length(y); Y2 = sum(y**2); Y = sum(y)
  p = - n / 2 * log(2 * pi) - (n - 1) * log(sigma) - 1/2 * log(n * s2 + sigma2) - 
    Y2 / (2 * sigma2) - m**2 / (2 * s2) + (s * Y / sigma + sigma * m / s ) ** 2 / (2 * (n * s2 + sigma2))
  return(p)
} 


sample_v <- function(group, group_of_cc, nc, nbar, mbar, alpha, z){
  p = rep(0,2)
  n0 = sum(group == group_of_cc); n1 = sum(group != group_of_cc)
  p[1] = lgamma(alpha + n0 - nbar ) - 
    lgamma(2 * alpha + n0 + n1 - mbar - nbar ) 
  log( genhypergeo(c(alpha * (1 + z), alpha + n1 - nbar , 
                     alpha + n1 - mbar),
                   c(alpha * (1 + z) + n1, 2 * alpha + n0 + n1 - nbar - mbar),
                   1) )
  
  p[2] = lgamma(alpha + n0 - nbar - nc) - 
    lgamma(2 * alpha + n0 + n1 - mbar - nbar - nc ) +
    log(z) - log(1-z) +
    log( genhypergeo(c(alpha * (1 + z), alpha + n1 - nbar - nc , alpha + n1 - mbar),
                     c(alpha * (1 + z) + n1, 2 * alpha + n0 + n1 - nbar - mbar  - nc),
                     1) )
  p = p - max(p)
  p = exp(p) / sum( exp(p) )
  if(sum(is.na(p)) == 2){p = c(0.5, 0.5); print("check convergence of 3F2")}
  return(sample(c(0,1), 1, prob = p))
}

sample_v_new <- function(group, group_of_cc, nc, nbar, mbar, alpha, z){
  p = rep(0,2)
  n0 = sum(group == group_of_cc); n1 = sum(group == group_of_cc)
  p[1] = lgamma(alpha + n0 - nbar ) - 
    lgamma(alpha + n0 + n1 - nbar ) +
    log(1-z) +
    log( genhypergeo(c(alpha + n1 - mbar - alpha * (z) + 
                         n0 - nbar, 
                       n0, n1),
                     c(alpha + n0 + n1 - mbar, 
                       alpha + n0 + n1 - nbar),
                     1) )
  
  p[2] = lgamma(alpha + n0 - nbar - nc ) - 
    lgamma(alpha + n0 + n1 - nbar - nc ) +
    log(z) +
    log( genhypergeo(c(alpha + n1 - mbar - alpha * (z) + 
                         n0 - nbar - nc, 
                       n0, n1),
                     c(alpha + n0 + n1 - mbar, 
                       alpha + n0 + n1 - nbar - nc),
                     1) )
  p = p - max(p)
  p = exp(p) / sum( exp(p) )
  if(sum(is.na(p)) == 2){p = c(0.5, 0.5); print("check convergence of 3F2")}
  return(sample(c(0,1), 1, prob = p))
}

sample_v_new2 <- function(group, group_of_cc, nc, nbar, mbar, alpha, z){
  p = rep(0,2)
  n0 = sum(group == group_of_cc); n1 = sum(group != group_of_cc)
  p[1] = lgamma(alpha + n0 - nbar ) - 
    lgamma(alpha + n0 + n1 - nbar ) +
    log(1-z) +
    log( genhypergeo_contfrac(c(alpha + n1 - mbar - alpha * (z) + 
                                  n0 - nbar, 
                                n0, n1),
                              c(alpha + n0 + n1 - mbar, 
                                alpha + n0 + n1 - nbar),
                              1))
  
  p[2] = lgamma(alpha + n0 - nbar - nc ) - 
    lgamma(alpha + n0 + n1 - nbar - nc ) +
    log(z) +
    log( genhypergeo_contfrac(c(alpha + n1 - mbar - alpha * (z) + 
                                  n0 - nbar - nc, 
                                n0, n1),
                              c(alpha + n0 + n1 - mbar, 
                                alpha + n0 + n1 - nbar - nc),
                              1) )
  p = p - max(p)
  p = exp(p) / sum( exp(p) )
  if(sum(is.na(p)) == 2){
    p = c(log(1-z),log(z)) 
    p = p - max(p)
    p = exp(p) / sum( exp(p) )
  }
  return(sample(c(0,1), 1, prob = p))
}


Furbi_GM_forMultiNA <- function(data, hyper = c(0.1, 0.5, 0.1), c_true = NA, totiter = 1000, verbose = FALSE){
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = totiter, clear = FALSE, width= 100)
  
  #data = scale(data, scale = FALSE)
  
  alpha = hyper[1]
  z = hyper[2]
  eps = hyper[3]
  
  
  n = dim(data)[1]
  p = dim(data)[2]
  miss_obs = seq(1,n)[rowSums(is.na(data)) > 0]
  
  group = c(rep(0,n))
  if(length(miss_obs)>0){
    for (i in miss_obs){
      group[i] = as.numeric(paste(which(is.na(data[i, ])),collapse=""))
    }
  }
  
  #initialization
  theta = matrix(rep(rnorm(p), n), ncol = p, byrow = TRUE)
  c = group
  v = rep(1, n)
  common = floor(length(unique(group))/2) + 1
  v[c%in%unique(group)[1:common] ] = 0
  
  if(z == 0){v = rep(0,n); c = sample(c(1,2,3),n, replace = TRUE) }
  if(z == 1){v = rep(1,n)}
  
  data_complete = data 
  for (pp in 1:p){
    data_complete[is.na(data_complete[,pp]),pp] = 0
    for(cc in unique(c)){
      theta[c == cc,pp] = mean(data_complete[c == cc,pp])
    }
  }
  
  u_theta = unique(theta)
  
  if(dim(u_theta)[1] == 1){ 
    cov = matrix(0, nrow = p, ncol= p)
  }else{
    cov = cov(u_theta)  
  }
  #Sigma0 = rinvwishart(length(unique(c))+ 3, diag(p)*3 + cov*length(unique(c)) )
  Sigma0 = rinvwishart(p + 3, diag(p)*p + cov*length(unique(c)) )
  
  
  for (iter in (1:totiter)){
    pb$tick()
    for (i in (1:n)){
      #identify from which DP it comes
      j = group[i]
      vi = v[i]
      if(vi == 1){
        #specific
        restjc = c[group == j & v == 1 & seq(1,n)!=i]
        theta_minus_i = theta[ group == j & v == 1 & seq(1,n)!=i,,drop = FALSE]
        z_temp = z
      }else{
        #common
        restjc = c[v==0 & seq(1,n)!=i]
        theta_minus_i = theta[v==0 & seq(1,n)!=i,,drop = FALSE]
        z_temp =  1 - z
      }
      if(length(restjc)>0){ #if it is not alone in its DP, sample from the Polya
        nc = table(restjc)
        c_unique = unique(restjc)
        prob = NULL; h = 1
        for (cc in c_unique){
          prob[h] = dmvnorm((data_complete[i,]), theta_minus_i[restjc==cc,,drop=FALSE][1,], diag(p), log = TRUE ) +
            log(nc[as.character(cc)])
          h = h + 1
        }
        prob[h] = 0 
        for (pp in 1:p){
          prob[h] = prob[h] + margnn(data_complete[i,pp], s2 = Sigma0[pp,pp])
        }
        prob[h] = prob[h] + log(alpha) + log(z_temp)
        prob = prob - max(prob)
        prob = exp(prob) / sum( exp(prob) )
        c[i] = sample(c(c_unique, setdiff(1:n,c[-i])[1]), 1, prob = prob)
        if(c[i] %in% c_unique){
          theta[i,] = theta_minus_i[restjc==c[i],,drop =FALSE][1,]
        }else{
          invS = solve(Sigma0 + diag(1, p))
          Sigman = Sigma0%*%invS%*%diag(1, p)
          theta[i,] = rmvnorm(1,Sigma0%*%invS%*%data_complete[i,],Sigman)
        }
      }else{
        invS = solve(Sigma0 + diag(1, p))
        Sigman = Sigma0%*%invS%*%diag(1, p)
        theta[i,] = rmvnorm(1,Sigma0%*%invS%*%data_complete[i,],Sigman)
      }
    }
    #sample v###################################################################
    if( z != 1 & z != 0 ){
      #sample those clusters observed only in one group
      v_old = v
      temp = (rowSums(table(c, group)>0)>1)
      c_corresp_to_v_to_sample = as.numeric(labels(temp)[temp==FALSE])
      for(cc in c_corresp_to_v_to_sample){
        group_of_cc = group[c==cc][1]
        nc = sum(c == cc)
        nbar = sum(v[group == group_of_cc] == 1)
        mbar = sum(v[group != group_of_cc] == 1)
        v[c == cc] = sample_v_new2(group, group_of_cc, nc, nbar, mbar, alpha, z)
      }
    }
    u_theta = unique(theta)
    #Sigma0 = cov(u_theta)
    if(dim(u_theta)[1] == 1){ 
      cov = matrix(0, nrow = p, ncol= p)
    }else{
      cov = cov(u_theta)  
    }
    #Sigma0 = rinvwishart(length(unique(c))+ 3, diag(p)*3 + cov(u_theta)*length(unique(c)) )
    Sigma0 = rinvwishart(p + 3, diag(p)*p + cov(u_theta)*length(unique(c)) )
    #sample theta, i.e, clusters parameters
    for(cc in unique(c)){
      nc = sum(c == cc)
      invS = solve(Sigma0 + diag(1/nc, p))
      Sigman = Sigma0%*%invS%*%diag(1/nc, p)
      mun = Sigma0%*%invS%*%colMeans(data_complete[c == cc,,drop=FALSE])
      theta[c == cc,] = matrix( rep(rmvnorm(1,mun,Sigman),nc),
                                ncol = p, byrow = TRUE)
      data_complete[c == cc,] = rmvnorm(nc,theta[c==cc,,drop=FALSE][1,],diag(p))*
        (is.na(data[c==cc,])) + 
        data_complete[c==cc,]*(!is.na(data[c==cc,])) 
    }
    print(c(iter, rand.index(c,c_true), length(unique(c)), sum(table(c) > n*0.1)))
    
    #if (floor(iter / 100) == (iter/100)){
    #       plot(as.data.frame(data_complete), col = c , pch = group)
    #}
    if(!is.na(c_true[1])){
      write.table(rand.index(c,c_true), file = "ri_saved.csv", append = TRUE, row.names = FALSE,
                  col.names= FALSE)
    }
    
    write.table(t(c), file = "c_saved.csv", append = TRUE, row.names = FALSE,
                col.names= FALSE)
    write.table(t(v), file = "v_saved.csv", append = TRUE, row.names = FALSE,
                col.names= FALSE)
    write.table(t(as.vector(Sigma0)), file = "Sigma0_saved.csv", append = TRUE, row.names = FALSE,
                col.names= FALSE)
    write.table(t(as.vector(unique(theta))), file = "theta_saved.csv", append = TRUE, row.names = FALSE,
                col.names= FALSE)
    write.table(data_complete, file = "data_complete_saved.csv", append = TRUE, row.names = FALSE,
                col.names= FALSE)
  }
}
