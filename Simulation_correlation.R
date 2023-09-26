library(VGAM)
corr_FURBI <- function(x, rho){
  res <- pbinorm(x, x, cov12 = rho)-pnorm(x)^2
  res <- res/(pnorm(x)*(1-pnorm(x)))
  return(res)
  
}
len <- 100
grid <- seq(-5,5, length = len)
rho <- c(-0.99, -0.5, 0, 0.5, 0.99)
corrs <- matrix(0, nrow = 5, ncol = len)
for(i in 1:len){
  for(j in 1:5){
    corrs[j,i] <- corr_FURBI(grid[i], rho[j])
  }
}

#plot
plot(grid,corrs[1,],col="orange", xlab = "Support",ylab = "Correlation",main = "",pch = 17,type = "b", xlim = c(-5,5),ylim = c(-1,1),cex.axis = 1.6, cex = 0.8, lwd = 4,lty = 2,las = 1,cex.lab = 1.4)
par(new = TRUE)
plot(grid,corrs[2,],col="red", xlab = "Support",ylab = "Density",main = "",pch = 15,type = "b", xlim = c(-5,5),ylim = c(-1,1),cex.axis = 1.6, cex = 0.5, lwd = 4,lty = 9,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,corrs[3,],col="brown", xlab = "Support",ylab = "Density",main = "",pch = 19,type = "b", xlim = c(-5,5),ylim = c(-1, 1),cex.axis = 1.6, cex = 0.2, lwd = 4,lty = 10,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,corrs[4,],col="blue", xlab = "Support",ylab = "Density",main = "",pch = 13,type = "b", xlim = c(-5,5),ylim = c(-1,1),cex.axis = 1.6, cex = 0.5, lwd = 4,lty = 1,las = 1,cex.lab = 1.4, ann = F, axes = F)
par(new = TRUE)
plot(grid,corrs[5,],col="black", xlab = "Support",ylab = "Density",main = "",pch = 16,type = "b", xlim = c(-5,5),ylim = c(-1, 1),cex.axis = 1.6, cex = 0.5, lwd = 4,lty = 1,las = 1,cex.lab = 1.4, ann = F, axes = F)
#legend(-5,0, legend = c("Rho = -0.99", "Rho = -0.5", "Rho = 0","Rho = 0.5", "Rho = 0.99"), col = c(rep("black",5)), lwd = 4, pch = c(17,15,19,13,16),lty = c(2,9,10,1,1), cex = 1, bty = "n")
