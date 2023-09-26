devtools::install_github("tommasorigon/LSBP")
library(LSBP)
library(ggplot2)

set.seed(2019)
mu_true = 10
Y1 <- rnorm(30, mu_true, 1) 
Y2 <- rnorm(100,-mu_true,1)
interval <- 0:30
grid <- seq(-15,15, by = 0.1) #evaluation grid
R <- 10000 #iterations
burnin <- 2000
#N <- 150 

#kernel and weights#################################################
data = data.frame( X = c(rep(1,30), rep(2,100)), Y = c(Y1, Y2))

fit_gibbs <- LSBP_Gibbs(Y ~ X|X, data=data, H=3, 
          control=control_Gibbs(R=R,burn_in=burnin,method_init="random") )
plot(data)
lines(data$X,colMeans(predict(fit_gibbs))) # Posterior mean

X2 = cbind(1,(round(quantile(data$X ,c(0.1,0.9)),2)))  
X1 = X2

pred_Gibbs <- array(0,c(R,length(grid),2)) 

for(r in 1:R){      # Cycle over the iterations of the MCMC chain
  for(i in 1:length(grid)){
    pred_Gibbs[r,i,] = LSBP_density(grid[i], X1, X2,
             fit_gibbs$param$beta_mixing[r,,],
             fit_gibbs$param$beta_kernel[r,,],
             fit_gibbs$param$tau[r,])
  }
}

estimate_Gibbs <- apply(pred_Gibbs,c(2,3),mean)
lower_Gibbs    <- apply(pred_Gibbs,c(2,3),function(x) quantile(x,0.05))
upper_Gibbs    <- apply(pred_Gibbs,c(2,3),function(x) quantile(x,0.95))



data.plot <- data.frame(
  density = c(c(estimate_Gibbs)),
  lower = c(c(lower_Gibbs)),
  upper = c(c(upper_Gibbs)),
  grid = rep(grid,2),
  x = (rep(c(1,2),each=301))
)
ggplot(data=data.plot) + geom_line(aes(x=grid,y=density), color = "blue",
                                   linetype = "dashed") + 
  facet_grid( 1~x,scales="free_y")+ 
  geom_ribbon(alpha=0.3,aes(x=grid,ymin=lower,ymax=upper), fill = "blue")+ 
  theme_classic() +
  xlab("")

#only weights##########################################################
data = data.frame( X = c(rep(1,30), rep(2,100)), Y = c(Y1, Y2))

fit_gibbs <- LSBP_Gibbs(Y ~ 1|X, data=data, H=3, 
                        control=control_Gibbs(R=R,burn_in=burnin,method_init="random") )
plot(data)
lines(data$X,colMeans(predict(fit_gibbs))) # Posterior mean

X2 = cbind(1,(round(quantile(data$X ,c(0.1,0.9)),2)))  
X1 = as.matrix(c(1,1))

pred_Gibbs <- array(0,c(R,length(grid),2)) 

for(r in 1:R){      # Cycle over the iterations of the MCMC chain
  for(i in 1:length(grid)){
    pred_Gibbs[r,i,] = LSBP_density(grid[i], X1, X2,
                        fit_gibbs$param$beta_mixing[r,,],
                        matrix(fit_gibbs$param$beta_kernel[r,]),
                        fit_gibbs$param$tau[r,])
  }
}

estimate_Gibbs <- apply(pred_Gibbs,c(2,3),mean)
lower_Gibbs    <- apply(pred_Gibbs,c(2,3),function(x) quantile(x,0.05))
upper_Gibbs    <- apply(pred_Gibbs,c(2,3),function(x) quantile(x,0.95))

data.plot <- data.frame(
  density = c(c(estimate_Gibbs)),
  lower = c(c(lower_Gibbs)),
  upper = c(c(upper_Gibbs)),
  grid = rep(grid,2),
  x = (rep(c(1,2),each=301))
)
ggplot(data=data.plot) + geom_line(aes(x=grid,y=density), color = "blue",
                                   linetype = "dashed") + 
  facet_grid( 1~x,scales="free_y")+ 
  geom_ribbon(alpha=0.3,aes(x=grid,ymin=lower,ymax=upper), fill = "blue")+ 
  theme_classic() +
  xlab("")