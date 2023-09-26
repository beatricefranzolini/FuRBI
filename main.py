# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:13:42 2021
"""
import os
from pathlib import Path
import time 

import numpy as np
import pandas as pd

import acrmgm as agm

import matplotlib.pyplot as plt


#set the working directory to source file location
mypath = Path().absolute()
os.chdir(mypath)

np.random.seed(0) #set seed for replicability 

#read data
data = pd.read_excel('data.xlsx', sheet_name = 'data2021M1', header = 0)
X1 = pd.DataFrame.to_numpy(data['X1'])[0:49]
X2 = pd.DataFrame.to_numpy(data['X2'])

X = np.append (X1, X2)
X1 = ( X1 - np.mean(X) ) / np.std(X)
X2 = ( X2 - np.mean(X) ) / np.std(X)

n = np.array ([len(X1), len(X2)])

#set model hyperparameters
theta = 1
alphatheta = 3
betatheta = 3

z = 0.5
alphaz = 2
betaz = 2

alpha1 = 2
beta1 = 4
lambda1 = 1 
#mu1 = 5.8591
mu1 =  ( 10 - np.mean(X) ) / np.std(X)

alpha2 = 2
beta2 = 4
lambda2 = 1 
#mu2 = 3.9731
mu2 =  ( -10 - np.mean(X) ) / np.std(X)

print(X1)
print(X2)
rho = 0

#set algorithm parameter
tot_iter = 10000
burnin = 5000
M = 10000
eps = 0.1
verbose = 'T' #if set to T,every 500 iterations a message is displayed

temp = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan

#initialize values
[c1, c2, V1, V2, theta1, theta2, tau1, tau2] = agm.initialize(n)

#create arrays to save simul 
theta1save = np.ones((n[0], tot_iter))
theta2save = np.ones((n[1], tot_iter))
tau1save = np.ones((n[0], tot_iter))
tau2save = np.ones((n[1], tot_iter))
c1save = np.ones((n[0], tot_iter))
c2save = np.ones((n[1], tot_iter))
V1save = np.ones((n[0], tot_iter))
V2save = np.ones((n[1], tot_iter))
rhosave = np.ones( tot_iter )
thetasave = np.ones( tot_iter )
zsave = np.ones( tot_iter )
V1new = np.ones( (2, tot_iter) )
V2new = np.ones( (2, tot_iter) )


grid = np.arange(np.min(np.append(X1,X2))-0.5, np.max(np.append(X1,X2))+0.5, \
                 step = (np.max(np.append(X1,X2))-np.min(np.append(X1,X2))+1) / 200)
estdens1 = np.zeros((len(grid),tot_iter))
estdens2 = np.zeros((len(grid),tot_iter))

simtau1p = np.random.gamma(alpha1, 1/ beta1, size = M)
simtau2p = np.random.gamma(alpha2, 1/ beta2, size = M)

#MCMC
start = time.time()

for _iter in range(tot_iter): 
    if (verbose=='T' and _iter > 0 and ((_iter)/100==int((_iter)/100))):
        print('simulation number', _iter)
        end = time.time()
        print('time in seconds for 100 iter:', end - start)
        print(np.unique(c2))
        print(sum(rhosave[0:(_iter-1)]>0)/_iter)
        start = time.time()
    
    #1.sample unique atoms 
    for j in np.intersect1d(c1, c2):
        
        nc1 = np.count_nonzero(c1 == j)
        nc2 = np.count_nonzero(c2 == j)
        
        mean1 = np.mean(X1[c1 == j])
        mean2 = np.mean(X2[c2 == j])
        
        var1 = np.var(X1[c1 == j])
        var2 = np.var(X2[c2 == j])
        
        [newtheta1, newtau1, newtheta2, newtau2] = \
            agm.samplepseudotie(nc1, nc2, mean1, mean2, var1, var2, \
                    alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,\
                        rho, M)
                
        theta1[c1 == j] = newtheta1
        tau1[c1 == j] = newtau1
        
        theta2[c2 == j] = newtheta2
        tau2[c2 == j] = newtau2   
        
    for j in np.setdiff1d(c1, c2):
        
        nc1 = np.count_nonzero(c1 == j)
        mean1 = np.mean(X1[c1 == j])
        var1 = np.var(X1[c1 == j])
        
        [newtheta1, newtau1] = \
            agm.sampletie(nc1, mean1, var1, alpha1, beta1, mu1, lambda1)
        
        theta1[c1 == j] = newtheta1
        tau1[c1 == j] = newtau1
        
    for j in np.setdiff1d(c2, c1):
        
        nc2 = np.count_nonzero(c2 == j)
        mean2 = np.mean(X2[c2 == j])
        var2 = np.var(X2[c2 == j])
        
        [newtheta2, newtau2] = \
            agm.sampletie(nc2, mean2, var2, alpha2, beta2, mu2, lambda2)
        
        theta2[c2 == j] = newtheta2
        tau2[c2 == j] = newtau2
        
    #2.sample v
    for j in np.setdiff1d(c1, c2):
        nc = np.count_nonzero( c1 == j )
        
        nbar = np.count_nonzero( V1[c1 != j] == 0 )
        mbar = np.count_nonzero( V2 == 0 )
        
        newv1 = agm.sampleVstar(n, nbar, nc, mbar, theta, z)
        V1[c1 == j] = newv1
        
    for j in np.setdiff1d(c2, c1):
        nc = np.count_nonzero( c2 == j )
        
        nbar = np.count_nonzero( V2[c2 != j]==0 )
        mbar = np.count_nonzero( V1 == 0 )
        
        newv2 = agm.sampleVstar([n[1], n[0]], nbar, nc, mbar, theta, z)
        V2[c2 == j] = newv1
    
    #3.sample c1
    for i in range(n[0]):
        #3.a sample label
        newc = agm.samplelabel1(i, V1, c1, V2, c2, theta1, theta2, tau1, tau2, \
                 alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,
                 rho, z, theta, simtau1p, X1)
        
        #3.b update atom value
        dropc = np.delete(c1, i)
        if np.count_nonzero( dropc == newc ) > 0: #ties
            
            c1[i] = newc
            theta1[i] = np.delete(theta1, i)[dropc == newc][0]
            tau1[i] = np.delete(tau1, i)[dropc == newc][0]
        
        elif newc>0: #pseudotie
            
            c1[i] = newc
            
            theta2temp = theta2[c2 == newc][0]
            tau2temp = tau2[c2 == newc][0]
            
            [newtheta, newtau] = \
            agm.samplepsuedotie_conditional(1, X1[i], 0, theta2temp, tau2temp,\
                     alpha1, beta1, mu1, mu2, lambda1, lambda2,\
                         rho, M)
            theta1[i] = newtheta
            tau1[i] = newtau
        
        else: #new atom ("tie" not yet observed)
            #search a free label
            c1[i] = np.setdiff1d( range(np.sum(n)), np.union1d(dropc, c2) )[0]
            
            [newtheta, newtau] = agm.sampletie(1, X1[i], 0, alpha1, beta1, mu1, lambda1)
            theta1[i] = newtheta
            tau1[i] = newtau
            
    
    #4.sample c2
    for i in range(n[1]):
        #3.a sample label
        newc = agm.samplelabel2(i, V1, c1, V2, c2, theta1, theta2, tau1, tau2, \
                 alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,
                 rho, z, theta, simtau2p, X2)
        
        #3.b update atom value
        dropc = np.delete(c2, i)
        if np.count_nonzero( dropc == newc ) > 0: 
            
            c2[i] = newc
            theta2[i] = np.delete(theta2, i)[dropc == newc][0]
            tau2[i] = np.delete(tau2, i)[dropc == newc][0]
        
        elif newc>0:
            
            c2[i] = newc
            
            theta1temp = theta1[c1 == newc][0]
            tau1temp = tau1[c1 == newc][0]
            
            [newtheta, newtau] = \
            agm.samplepsuedotie_conditional(1, X2[i], 0, theta1temp, tau1temp,\
                     alpha2, beta2, mu2, mu1, lambda2, lambda1,\
                         rho, M)
            theta2[i] = newtheta
            tau2[i] = newtau
        
        else:
            #search a free label
            c2[i] = np.setdiff1d( range(np.sum(n)), np.union1d(dropc, c1) )[0]
            
            [newtheta, newtau] = agm.sampletie(1, X2[i], 0, alpha2, beta2, mu2, lambda2)
            theta2[i] = newtheta
            tau2[i] = newtau
        
    #5.sample z
    z = agm.samplez(z, n, c1, c2, V1, V2, theta, alphaz, betaz, eps)
    
    #6.sample theta
    theta = agm.sampletheta(theta, n, c1, c2, V1, V2, z, alphatheta, betatheta)
    
        
    #7.sample rho
    rho = agm.samplerho(rho, c1, c2, theta1, theta2, tau1, tau2,\
                       lambda1, lambda2, mu1, mu2, eps)
    
    #save simul
    theta1save[:,_iter] = theta1
    theta2save[:,_iter] = theta2
    tau1save[:,_iter] = tau1
    tau2save[:,_iter] = tau2
    V1save[:,_iter] = V1
    V2save[:,_iter] = V2
    c1save[:,_iter] = c1
    c2save[:,_iter] = c2
    rhosave[_iter] = rho
    thetasave[_iter] = theta
    zsave[_iter] = z
    
    
#compute predictive
print("mcmc ended, density estimation has started")
start = time.time()
for _iter in range(burnin,tot_iter): 
    if (verbose=='T' and ((_iter)/100==int((_iter)/100))):
        print('simulation number', _iter)
        end = time.time()
        print('time in seconds for 100 iter:', end - start)
        start = time.time()
    c1 = c1save[:,_iter]
    c2 = c2save[:,_iter]
    V1 = V1save[:,_iter]
    V2 = V2save[:,_iter]
    theta1 = theta1save[:,_iter]
    theta2 = theta2save[:,_iter]
    tau1 = tau1save[:,_iter]
    tau2 = tau2save[:,_iter]
    #theta = thetasave[_iter]
    theta = 20
    z = zsave[_iter]
    rho = rhosave[_iter]
    
    [v1new, v2new] = agm.predweights(V1, V2, theta, z, M)
    [dens0, dens1] = agm.predictive(grid, c1, c2, V1, V2, theta1, theta2, tau1, tau2, \
               alpha1, alpha2, beta1, beta2, lambda1, lambda2, mu1, mu2,
               theta, z, rho, simtau1p)
    estdens1[:, _iter] = v1new[0] * dens0 + v1new[1] * dens1
    estdens1[:, _iter] = dens1
    [dens0, dens1] = agm.predictive(grid, c2, c1, V2, V1, theta2, theta1, tau2, tau1, \
               alpha2, alpha1, beta2, beta1, lambda2, lambda1, mu2, mu1,
               theta, z, rho, simtau2p)
    estdens2[:, _iter] = v2new[0] * dens0 + v2new[1] * dens1
    estdens2[:, _iter] = dens1

estdens1 = estdens1 / np.std(X)
estdens2 = estdens2 / np.std(X)

grid = grid * np.std(X) + np.mean(X)
density1 = np.mean(estdens1[:,burnin:tot_iter], axis=1)
density1_1 = np.quantile(estdens1[:,burnin:tot_iter], q = 0.01, axis=1)
density1_5 = np.quantile(estdens1[:,burnin:tot_iter], q = 0.05, axis=1)
density1_95 = np.quantile(estdens1[:,burnin:tot_iter], q = 0.95, axis=1)
density1_99 = np.quantile(estdens1[:,burnin:tot_iter], q = 0.99, axis=1) 
plt.hist(X1 * np.std(X) + np.mean(X), density = True, stacked=True, color='grey', bins=10, alpha=0.3, label = "Empirical distribution") 
plt.fill_between(grid, density1_5, \
                     density1_95, alpha=0.7, color ='#bbddff', label = "95% interval" )
plt.fill_between(grid, density1_1, \
                     density1_5, alpha=0.7, color ='#ffe4e1', label = "99% interval" )
plt.fill_between(grid, density1_95, \
                     density1_99, alpha=0.7, color ='#ffe4e1')
plt.plot(grid, density1, color= 'blue', linewidth = 2, linestyle='--',\
     label="Estimated density")
plt.legend()
plt.ylim((0,0.16))
plt.xlim((-10,20))
plt.savefig('stocksGMacrm', dpi=300)
plt.show()



density2 = np.mean(estdens2[:,burnin:tot_iter], axis=1)
density2_1 = np.quantile(estdens2[:,burnin:tot_iter], q = 0.01, axis=1)
density2_5 = np.quantile(estdens2[:,burnin:tot_iter], q = 0.05, axis=1)
density2_95 = np.quantile(estdens2[:,burnin:tot_iter], q = 0.95, axis=1)
density2_99 = np.quantile(estdens2[:,burnin:tot_iter], q = 0.99, axis=1) 
plt.hist(X2 * np.std(X) + np.mean(X), density = True, stacked=True, color='grey', bins=10, alpha=0.3, label = "Empirical distribution") 
plt.fill_between(grid, density2_5, \
                     density2_95, alpha=0.7, color ='#bbddff', label = "95% interval" )
plt.fill_between(grid, density2_1, \
                     density2_5, alpha=0.7, color ='#ffe4e1', label = "99% interval" )
plt.fill_between(grid, density2_95, \
                     density2_99, alpha=0.7, color ='#ffe4e1')
plt.plot(grid, density2, color= 'blue', linewidth = 2, linestyle='--',\
     label="Estimated density")
plt.legend()
plt.ylim((0,0.16))
plt.xlim((-10,20))
plt.savefig('commoditiesGMacrm', dpi=300)
plt.show()

np.savez("GMACRMrealstandard", theta1save = theta1save,\
         theta2save = theta2save,\
    tau1save = tau1save,\
    tau2save = tau2save,\
    V1save= V1save,\
    V2save = V2save,\
    c1save = c1save,\
    c2save= c2save,\
    rhosave = rhosave,\
    thetasave = thetasave,\
    zsave = zsave,\
        X1=X1,\
            X2=X2,\
       estdens1 = estdens1,\
       estdens2 = estdens2) 
    

results = np.load("GMACRMrealstandard.npz")
V1save = results['V1save']
V2save = results['V2save']
theta1save = results['theta1save']
theta2save  = results['theta2save']
tau1save = results['tau1save']
tau2save = results['tau2save']
c1save = results['c1save']
c2save = results['c2save']
rhosave = results['rhosave']
thetasave = results['thetasave']
zsave = results['zsave']
estdens1 = results['estdens1']
estdens2 = results['estdens2']
