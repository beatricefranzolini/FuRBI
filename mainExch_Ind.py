# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:24:39 2021
"""

import os
from pathlib import Path
import time 

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import t as student


which = 1 #### 0 for independence, 1 for exchangeable
#set the working directory to source file location
mypath = Path().absolute()
os.chdir(mypath)

np.random.seed(0) #set seed for replicability 

#I. DATA########################################################################
data = pd.read_excel('data.xlsx', sheet_name = 'data2021M1', header = 0)
X1 = pd.DataFrame.to_numpy(data['X1'])[0:49]
X2 = pd.DataFrame.to_numpy(data['X2'])
if  (which == 0) :
    #X = X2
    #mu = 7.8851
    X = X1
    X = (X - np.mean(np.append(X1,X2)) ) / np.std(np.append(X1,X2))
    #mu = (3.9731 - np.mean(np.append(X1,X2)) ) / np.std(np.append(X1,X2))
    mu = (5.8591 - np.mean(np.append(X1,X2)) ) / np.std(np.append(X1,X2))
else:
    X = np.append(X1,X2)
    X = (X - np.mean(X) ) / np.std(X)
    mu = (4.8617 - np.mean(np.append(X1,X2)) ) / np.std(np.append(X1,X2))
n = len(X)

#II. HYPERPARAMETER##############################################################
#1.base dist
lam = 1
alpha = 2
beta = 4

#2.hyperparameters for the concentration
#omega:
a = 3
b = 3 

#set algorithm parameter
tot_iter = 10000
burnin = 5000
M = 10000
eps = 0.1
verbose = 'T' #if set to T,every 100 iterations a message appears

#III. INITIALIZATION############################################################
c = np.arange(n)
K = np.max(c)+1
label = np.arange(n)
nk = np.ones(n)
theta = np.random.normal(mu, 1, n)
sigma = 1 / np.random.gamma(alpha, 1/beta, n)
omega = 1

#IV. FUNCTIONS##################################################################

def samplec(x, K, label, nk, theta, sigma, omega, mu, lam, alpha, beta):
    _temp = len(label)
    p = np.empty((2,_temp+1))
    p[0,0:_temp] = label
    p[0,_temp] = K
    for k in range(_temp):
        p[1,k] = nk[k] * norm.pdf(x, theta[k], sigma[k]**(1/2))
    p[1,_temp] = alpha * student.pdf(x, df = 2 * alpha,\
                                    loc = mu, \
                scale = (beta * (lam + 1) / (alpha * lam))**(1/2))
    p[1,:] = p[1,:] / np.sum(p[1,:])
    c = np.random.choice(p[0,:], p = p[1,:])
    return(c)

def sampleomega(omega, a, b, K, n):
    eta = np.random.beta(omega+1, n)
    p = np.array ([a + K - 1, n * (b - np.log(eta))])
    p = p / np.sum(p)
    mix = np.random.choice([0,1],p=p)
    return( (1-mix) * np.random.gamma(a + K, 1 / (b-np.log(eta))) +\
             mix * np.random.gamma(a + K - 1, 1 / (b-np.log(eta))))

#create arrays to save simul 
Ksave = np.ones(tot_iter) * np.nan
csave = np.ones((n,tot_iter)) * np.nan
thetasave = np.ones((n,tot_iter)) * np.nan
labelsave = np.ones((n,tot_iter)) * np.nan
sigmasave = np.ones((n,tot_iter)) * np.nan
nksave = np.ones((n,tot_iter)) * np.nan
omegasave = np.ones(tot_iter) * np.nan

#grid = np.arange(np.min(X)-0.5, np.max(X)+0.5, step = (np.max(X)-np.min(X)+1) / 200)
grid = np.arange(-2.1267283449534635-0.5, 2.8250032481699052+0.5, step = 5.951731593123369 / 200)
estdens = np.zeros((len(grid),tot_iter))

start = time.time()
for sim in range(tot_iter):
    if (verbose=='T' and sim > 0 and ((sim)/100==int((sim)/100))):
        print('simulation number', sim)
        end = time.time()
        print('time in seconds for 100 iter:', end - start)
        start = time.time()
    for i in range(n):
        
        index = np.where(label==c[i]) 
        nk[index] = nk[index] - 1
        c[i] = samplec(X[i], K, label, nk, theta, sigma, omega, mu, lam, alpha, beta)
        if (c[i] < K):
            index = np.where(label==c[i])
            nk[index] = nk[index] + 1
        else:
            c[i] = np.setdiff1d(range(2*n), label)[0]
            label = np.append(label, c[i])
            
            nk = np.append(nk, 1)
            
            a = alpha + 1 / 2
            b = beta + lam / (2 * (lam + 1)) * (X[i] - mu)
            s = 1 / np.random.gamma(a, 1/b)
            m = 1 / (1  + lam ) * X[i] + \
                lam  / (1  + lam ) * mu
            t = 1 * 1/s + lam * 1/s
            mm =  np.random.normal(m, t ** (-1/2))
            theta = np.append(theta, mm)
            sigma = np.append(sigma, s)
        
    
    K = np.count_nonzero(nk)
    empty = np.where(nk==0)
    nk = np.delete(nk, empty)
    label = np.delete(label, empty)
    theta = np.delete(theta, empty)
    sigma = np.delete(sigma, empty)

    for k in range(K):
        nc = np.count_nonzero(c==label[k])
        obs = X[np.where(c==label[k])]
        a = alpha + nc / 2
        b = beta + 1/2 * nc * np.var(obs) +\
            nc * lam / (2 * (lam + nc)) * (np.mean(obs) - mu) ** 2
        sigma[k]= 1 / np.random.gamma(a, 1/b)
        
        m = nc / (nc  + lam ) * np.mean(obs) + \
            lam  / (nc  + lam ) * mu
        t = nc * 1/sigma[k] + lam * 1 / sigma[k]
        theta[k] = np.random.normal(m, t ** (-1/2))
    
    omega = sampleomega(omega, a, b, K, n)
    
    Ksave[sim] = K
    csave[:,sim] = c
    thetasave[0:K,sim] = theta
    sigmasave[0:K,sim] = sigma
    nksave[0:K,sim] = nk
    omegasave[sim] = omega
    labelsave[0:K,sim] = label

#compute predictive
print("mcmc ended, density estimation has started")
start = time.time()

for sim in range(burnin,tot_iter):
    if (verbose=='T' and ((sim)/100==int((sim)/100))):
        print('simulation number', sim)
        end = time.time()
        print('time in seconds for 100 iter:', end - start)
        start = time.time()
    for i in range(len(grid)):
        K = int(Ksave[sim])
        theta = thetasave[0:K,sim]
        sigma = sigmasave[0:K,sim]
        nk = nksave[0:K,sim]
        omega = omegasave[sim]
        estdens[i,sim] = (np.sum(nk * norm.pdf(grid[i], loc = theta, \
                                          scale = sigma**(1/2))) +
                 omega * student.pdf(grid[i], df = 2 * alpha,\
                                                 loc = mu, \
                    scale = (beta * (lam + 1) / (alpha * lam))**(1/2)))\
                     /(n + omega) 

estdens = estdens / np.std(np.append(X1,X2))
#estdens = estdens / np.std(X1)
grid = grid * np.std(np.append(X1,X2)) + np.mean(np.append(X1,X2))
#grid = grid * np.std(X1) + np.mean(X1)
#plots
if  (which == 0) :
    np.savez("INDrealX2standard", thetasave = thetasave,\
         Ksave = Ksave,\
             csave = csave,\
                 sigmasave = sigmasave,\
                     nksave = nksave,\
                         omegasave = omegasave,\
                             estdens = estdens,\
                              labelsave = labelsave)
        
    density = np.mean(estdens[:,burnin:tot_iter], axis=1)
    density_1 = np.quantile(estdens[:,burnin:tot_iter], q = 0.01, axis=1)
    density_5 = np.quantile(estdens[:,burnin:tot_iter], q = 0.05, axis=1)
    density_95 = np.quantile(estdens[:,burnin:tot_iter], q = 0.95, axis=1)
    density_99 = np.quantile(estdens[:,burnin:tot_iter], q = 0.99, axis=1) 
    plt.hist(X* np.std(np.append(X1,X2)) + np.mean(np.append(X1,X2)), density = True, stacked=True, color='grey', bins=10, alpha=0.3,\
             label = "Empirical distribution") 
    plt.fill_between(grid, density_5, \
                         density_95, alpha=0.7, color ='#bbddff', label = "95% interval" )
    plt.fill_between(grid, density_1, \
                         density_5, alpha=0.7, color ='#ffe4e1', label = "99% interval" )
    plt.fill_between(grid, density_95, \
                         density_99, alpha=0.7, color ='#ffe4e1')
    plt.plot(grid, density, color= 'blue', linewidth = 2, linestyle='--',\
         label="Estimated density Ind.")
    plt.ylim((0,0.16))
    plt.xlim((-10,20))
    plt.legend()
    plt.savefig('comm_Ind', dpi=300)
    plt.show()
else:
    np.savez("EXCHrealstandard", thetasave = thetasave,\
         Ksave = Ksave,\
             csave = csave,\
                 sigmasave = sigmasave,\
                     nksave = nksave,\
                         omegasave = omegasave,\
                             estdens = estdens,\
                              labelsave = labelsave)
        
    density = np.mean(estdens[:,burnin:tot_iter], axis=1)
    density_1 = np.quantile(estdens[:,burnin:tot_iter], q = 0.01, axis=1)
    density_5 = np.quantile(estdens[:,burnin:tot_iter], q = 0.05, axis=1)
    density_95 = np.quantile(estdens[:,burnin:tot_iter], q = 0.95, axis=1)
    density_99 = np.quantile(estdens[:,burnin:tot_iter], q = 0.99, axis=1) 
    plt.hist(X1, density = True, stacked=True, color='grey', bins=10, alpha=0.3,\
             label = "Empirical distribution") 
    plt.fill_between(grid, density_5, \
                         density_95, alpha=0.7, color ='#bbddff', label = "95% interval" )
    plt.fill_between(grid, density_1, \
                         density_5, alpha=0.7, color ='#ffe4e1', label = "99% interval" )
    plt.fill_between(grid, density_95, \
                         density_99, alpha=0.7, color ='#ffe4e1')
    plt.plot(grid, density, color= 'blue', linewidth = 2, linestyle='--',\
         label="Estimated density Exch.")
    plt.ylim((0,0.16))
    plt.xlim((-10,20))
    plt.legend()
    plt.savefig('stocksExch', dpi=300)
    plt.show()
    
    plt.hist(X2, density = True, stacked=True, color='grey', bins=10, alpha=0.3,\
             label = "Empirical distribution") 
    plt.fill_between(grid, density_5, \
                         density_95, alpha=0.7, color ='#bbddff', label = "95% interval" )
    plt.fill_between(grid, density_1, \
                         density_5, alpha=0.7, color ='#ffe4e1', label = "99% interval" )
    plt.fill_between(grid, density_95, \
                         density_99, alpha=0.7, color ='#ffe4e1')
    plt.plot(grid, density, color= 'blue', linewidth = 2, linestyle='--',\
         label="Estimated density Exch.")
    plt.ylim((0,0.16))
    plt.xlim((-10,20))
    plt.legend()
    plt.savefig('commoditiesExch', dpi=300)
    plt.show()