rm(list=ls())

source("Furbi_GM_forMultiNA.R")
out = sim_data(1) #1,2,3,4
data = out$data; c_true = out$c_true; n=dim(data)[1]
group = c(rep(0,n))
miss_obs = seq(1,n)[rowSums(is.na(data)) > 0]
if(length(miss_obs)>0){
  for (i in miss_obs){
    group[i] = as.numeric(paste(which(is.na(data[i, ])),collapse=""))
  }
}
plot(as.data.frame(out$data_true), col = c_true , pch = group)
vis_miss(as.data.frame(data))

missingplot = data.frame("cluster" = rep(c(1,2,3,4),3),
                         "variable" = c(rep("V1", 4), rep("V2", 4),rep("V3", 4)),
                         "miss" = c(0, 12))
i = 1
for (pp in 1:3){
  for (cc in 1:4){
    missingplot$miss[i] = sum(is.na(data[c_true==cc, pp])) / length(data[c_true==cc, pp]) 
    i = i+1
  }
}

ggplot(missingplot, aes(x = cluster, y = variable, fill = miss)) +
  theme_minimal()+
  geom_tile(color = "black") +
  scale_fill_continuous(breaks=c(0,0.25,0.5,0.75,1),
                        labels=c(0,0.25,0.5,0.75,1),
                        limits=c(0,1)) +
  coord_fixed()

#mice + k-means##################################
set.seed(0)
tempData <- mice(data)
completedData <- complete(tempData,1)

set.seed(0)
df = as.matrix(completedData)
# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     main = "Average Silhouettes - kmeans sim n.1")

set.seed(2)
k_est = which(avg_sil_values == max(avg_sil_values)) + 1
k_est
c = kmeans(df, centers = k_est)$cluster
rand.index(c, c_true)
#nFurbi#################################
set.seed(0)
Furbi_GM_forMultiNA(data, hyper = c(0.1, 0.5, 0.1), c_true=c_true, totiter = 25000)
set.seed(0) #Dirichlet
Furbi_GM_forMultiNA(df, hyper = c(0.1, 0, 0.1), c_true=c_true, totiter = 25000)
set.seed(1)
Furbi_GM_forMultiNA(data, hyper = c(0.1, 0.2, 0.1), c_true=c_true, totiter = 25000)
set.seed(2)
Furbi_GM_forMultiNA(data, hyper = c(0.1, 0.8, 0.1), c_true=c_true, totiter = 25000)

c_saved = read.table("c_saved.csv", header = FALSE)
K = apply(c_saved, 1, function(x) length(unique(x)))
mean(K[12500:25000])
ri_saved = read.table("ri_saved.CSV", header = FALSE)
mean(ri_saved[12500:25000,])

plot(K, type = "l", xlab="MCMC iterations")
title(main="Number of clusters - Trace plot",
      cex.lab=0.75)
abline(v = 12500, col="red", lwd=3, lty=2)

plot(K[12501:25000], type = "l",  xlab="MCMC iterations - after burn-in")
title(main="Number of clusters - Trace plot",
      cex.lab=0.75)

plot(ri_saved[,], type = "l", xlab="MCMC iterations")
title(main="Rand Index - Trace plot",
      cex.lab=0.75)
abline(v = 12500, col="red", lwd=3, lty=2)
plot(ri_saved[12501:25000,], type = "l", xlab="MCMC iterations - after burn-in")
title(main="Rand Index - Trace plot",
      cex.lab=0.75)
ESS(K[12501:25000])/ 12500
ESS(ri_saved[12501:25000,])/ 12500

ri_adj = NULL
for (iter in 12501:25000){
  ri_adj = c(ri_adj, adj.rand.index(as.numeric(c_saved[iter,]), c_true))
}
mean(ri_adj)
