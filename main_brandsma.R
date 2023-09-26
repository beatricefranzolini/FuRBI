library(scatterplot3d)
#missing data
library(progress)#to draw the progress bar
library(hypergeo) #for hypergeometric function
library(mvtnorm) #for multivariate normal
library(fossil) #to compute rand indexes
library(mice) #to impute missing data
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(VIM)
library(visdat)
library(mice)
library(mcclust.ext) #to minimize VI
library(T4cluster) #to compute psm with psm() function
data = brandsma

data = data.frame(brandsma$iqv, brandsma$iqp, brandsma$lpr, brandsma$apr)
temp = apply(data, 1, FUN=is.na)
sum(colSums(temp)>0)
sum(is.na(data))
vis_miss(data)

set.seed(0)
n = dim(data)[1]
group = c(rep(0,n))
miss_obs = seq(1,n)[rowSums(is.na(data)) > 0]
for (i in miss_obs){
  group[i] = as.numeric(paste(which(is.na(data[i, ])),collapse=""))
}
table(group)


margnn_multi_cor <- function(y, sigma , rho, log = TRUE){
  n = dim(y)[1]
  ybarn = colSums(y)
  y2barn = sum(y**2)
  Sigma0 = sigma*rho*t(sigma)
  Sigmanew =  solve(n + solve(Sigma0)) 
  
  p = - n / 2 * log(2 * pi) - 1/2 * log(2 * pi * det(Sigma0)) +
    1 / 2 * log(2 * pi * det(Sigmanew) ) +
    1 / 2 * (t(ybarn) %*% Sigmanew %*% ybarn - y2barn)
  return(p)
} 

sample_v <- function(group, group_of_cc, nc, nbar, mbar, alpha, z){
  p = rep(0,2)
  n0 = sum(group == group_of_cc); n1 = sum(group != group_of_cc)
  p[1] = lgamma(alpha + n0 - nbar - nc ) - 
    lgamma(2 * alpha + n0 + n1 - mbar - nbar - nc ) +
    log(z) - log(1-z) +
    log( genhypergeo(c(alpha * (1 + z), alpha + n1 - nbar - nc , alpha + n1 - mbar),
                     c(alpha * (1 + z) + n1, 2 * alpha + n0 + n1 - nbar - mbar - nc),
                     1) )
  
  p[2] = lgamma(alpha + n0 - nbar) - 
    lgamma(2 * alpha + n0 + n1 - mbar - nbar) +
    log( genhypergeo(c(alpha * (1 + z), alpha + n1 - nbar , alpha + n1 - mbar),
                     c(alpha * (1 + z) + n1, 2 * alpha + n0 + n1 - nbar - mbar),
                     1) )
  p = p - max(p)
  p = exp(p) / sum( exp(p) )
  if(sum(is.na(p)) == 2){p = c(0.5, 0.5); print("check convergence of 3F2")}
  return(sample(c(0,1), 1, prob = p))
}

sample_rho <-function(theta, p, eps, rho, sigma){
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE) 
  rho_pro = diag(p) * 0.5
  rho_pro[col>row] = runif(p * (p - 1) / 2, rho[col>row] - eps, rho[col>row] + eps)
  rho_pro = rho_pro + t(rho_pro)
  sigma_pro = matrix(rep(rgamma(p, 10*sigma[1,], rep(10,p)), p), ncol = p, byrow = TRUE)
  
  if(prod(rho_pro[col>row]> - 0.95) == 1 & prod(rho_pro[col>row]< 0.95) == 1){
    ratio = (sum(dmvnorm(unique(theta), rep(0,p), sigma_pro*rho_pro*t(sigma_pro), log = TRUE)) - 
               sum(dmvnorm(unique(theta), rep(0,p),  sigma*rho*t(sigma), log = TRUE)))
    if(!is.na(ratio)){
      if(runif(1)<exp(ratio) ){
        return(list("rho" = rho_pro, "sigma" = sigma_pro))
      }
    }
  }
  return(list("rho" = rho, "sigma" = sigma))
}

totiter = 25000

data = scale(data)
alpha = 0.1
sigma2 = 1
z = 0.5
eps = 0.1

#data = scale(data)
n = length(group)
p = dim(data)[2]
miss_obs = seq(1,n)[rowSums(is.na(data)) > 0]


rho = diag(p)
sigma = matrix(rep(rep(1, p), p), ncol = p, byrow = TRUE) 

u_group = unique(group)
#initialization
if (z == 1){
  c = group
  v = rep(0, n)
}else{
  c = sample(c(0,1), n, replace = TRUE)
  v =  rep(1, n)
  if (z!=0){
    for (j in unique(group)){
      ngroup = sum(group==j)
      if(ngroup>1){
        not_shared = sample(which(group==j), floor(ngroup/2))
        c[not_shared] = j + 2
        v[not_shared] = 0
      }
    }
  }
}

theta = matrix(rep(rnorm(p),n), ncol = p, byrow = TRUE)
data_complete = data
for (i in miss_obs){
  p_miss = which(is.na(data_complete[i,]))
  p_nonmiss = which(!is.na(data_complete[i,]))
  S = diag(p)
  S11 = S[p_miss, p_miss, drop = FALSE]
  S12 = S[p_miss,, drop = FALSE][, p_nonmiss,drop = FALSE]
  S22 = S[p_nonmiss, p_nonmiss, drop = FALSE]
  non_miss = data_complete[i,p_nonmiss]
  mu_incomplete = theta[c == c[i],,drop=FALSE][1, p_miss ]
  mu_complete = theta[c == c[i],,drop=FALSE][1, p_nonmiss ]
  data_complete[i,is.na(data_complete[i,])] = 
    rmvnorm(1, mu_incomplete + S12%*%solve(S22)%*%(non_miss - mu_complete),
            S11 - S12%*%solve(S22)%*%t(S12))
}

rho = cor(data_complete)
print(sum(is.na(data_complete)))

for (iter in (1:totiter)){
  print("mcmc running")
  #sample c###################################################################
  for (i in (1:n)){
    j = group[i]
    vi = v[i]
    if(vi == 1){
      restjc = c[v == 1 & seq(1,n)!=i]
      y_minus_i = data_complete[v == 1 & seq(1,n)!=i,,drop = FALSE]
    }else{
      restjc = c[group == j & v==0 & seq(1,n)!=i]
      y_minus_i = data_complete[group == j & v==0 & seq(1,n)!=i,,drop = FALSE]
    }
    if(length(restjc)>0){
      nc = table(restjc)
      c_unique = unique(restjc)
      prob = NULL; h = 1
      for (cc in c_unique){
        prob[h] = margnn_multi_cor(rbind(y_minus_i[t(restjc == cc),], data_complete[i,]), sigma , rho, log = TRUE ) - 
          margnn_multi_cor(y_minus_i[t(restjc == cc),,drop = FALSE], sigma , rho, log = TRUE ) +
          log(nc[as.character(cc)])
        h = h + 1
      }
      prob[h] = margnn_multi_cor(data_complete[i,,drop = FALSE], sigma , rho, log = TRUE) + log(alpha)
      prob = prob - max(prob)
      prob = exp(prob) / sum( exp(prob) )
      c[i] = sample(c(c_unique, setdiff(1:n,c[-i])[1]), 1, prob = prob)
    }
  }
  #sample clusters' parameters###############################################################
  tab = aggregate(list(data_complete), by = list(c), FUN = function(x) { my.mean = mean(x, na.rm = TRUE) } )
  for (cc in unique(c)){
    sample_means_by_cluster = tab[tab[,1]==cc,2:(p+1)]
    S = sigma*rho*t(sigma)
    invM = solve( diag(p)*sum(c==cc) + solve(S) )
    b = colSums( data_complete[c==cc,,drop = FALSE])
    theta[c==cc,] = matrix(rep(rmvnorm(1, invM%*%b, invM),
                               sum(c == cc)), byrow = TRUE , nrow = sum(c == cc))
  }
  #sample v###################################################################
  if( z != 1 & z!=0){
    temp = (rowSums(table(c, group)>0)>1)
    c_corresp_to_v_to_sample = as.numeric(labels(temp)[temp==FALSE])
    for(cc in c_corresp_to_v_to_sample){
      group_of_cc = group[c==cc][1]
      nc = sum(c == cc)
      nbar = sum(v[group == group_of_cc] == 0)
      mbar = sum(v[group != group_of_cc] == 0)
      v[c == cc]=sample_v(group, group_of_cc, nc, nbar, mbar, alpha, z) 
    }
  }
  #sample rho#################################################################
  out = sample_rho(theta, p, eps, rho, sigma)
  rho = out$rho; sigma = out$sigma
  #sample missing data #######################################################
  data_complete = data
  for (i in miss_obs){
    p_miss = which(is.na(data_complete[i,]))
    p_nonmiss = which(!is.na(data_complete[i,]))
    S = diag(p)
    S11 = S[p_miss, p_miss, drop = FALSE]
    S12 = S[p_miss,, drop = FALSE][, p_nonmiss,drop = FALSE]
    S22 = S[p_nonmiss, p_nonmiss, drop = FALSE]
    non_miss = data_complete[i,!is.na(data_complete[i,])]
    mu_incomplete = theta[c == c[i],,drop = FALSE][1, p_miss ]
    mu_complete = theta[c == c[i],,drop = FALSE][1, p_nonmiss ]
    data_complete[i,is.na(data_complete[i,])] = 
      rmvnorm(1, mu_incomplete + S12%*%solve(S22)%*%(non_miss - mu_complete),
              S11 - S12%*%solve(S22)%*%t(S12))
    
  }
  print(c(iter, length(unique(c))))
  print(table(c))
  #if (floor(iter / 100) == (iter/100)){
  #  plot(as.data.frame(data_complete), col = c , pch = group)
  #}
  # write.table(t(c), file = "c_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  # write.table(t(v), file = "v_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  # write.table(rho, file = "rho_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  # write.table(sigma, file = "sigma_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  # write.table(unique(theta), file = "theta_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  # write.table(data_complete, file = "data_complete_saved.csv", append = TRUE, row.names = FALSE,
  #             col.names= FALSE)
  end_time <- Sys.time()
  print(end_time - start_time)
}

c_saved = read.table("c_saved.csv", header = FALSE)
c_saved = c_saved + 1

dim(c_saved)

library(dplyr)
data_complete %>%
  group_by(c) %>%
  summarise_at(c(1,2,3,4), list(name = mean))
data_complete %>%
  group_by(c) %>%
  summarise_at(c(1,2,3,4), list(name = median))

plot(data_complete[1:4], col = data_complete$c)

data_complete = data_complete %>% 
  rename(
    IQV = x,
    IQP = x.1,
    LPR = x.2,
    APR = x.3
  )

data = data.frame(brandsma$iqv, brandsma$iqp, brandsma$lpr, brandsma$apr)
prova = scale(data)
data_comp_notscaled = t(t(data_complete[,1:4]) * attr(prova, "scaled:scale") + attr(prova, "scaled:center"))
data_comp_notscaled = as.data.frame(data_comp_notscaled)
data_comp_notscaled$c = data_complete$c
plot(data_comp_notscaled[1:4], col = data_comp_notscaled$c)


c_saved = as.matrix(c_saved)
c_saved <- mapply(c_saved, FUN=as.integer)
class(c_saved)<- "numeric"

psm = comp.psm(c_saved[10001:50000,])
VImin = minVI(psm, c_saved[10001:50000,], method = "draws")


results = read.table("results.csv", header = FALSE)

colnames(results) = c("ID", "IQV", "IQP", "LPR", "APR", "CLUSTER")
prop.table(table(results$CLUSTER))
results %>%
  group_by(CLUSTER) %>%
  summarise_at(vars("IQV", "IQP", "LPR", "APR"), list(name = mean))

results %>%
  group_by(CLUSTER) %>%
  summarise_at(vars("IQV", "IQP", "LPR", "APR"), list(name = sd))

plot(results[,2:5], col = results$CLUSTER, pch = results$CLUSTER, cex = 1.5 )

