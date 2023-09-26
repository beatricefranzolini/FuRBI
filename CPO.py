# -*- coding: utf-8 -*-
"""
Created on Thu May 13 13:47:29 2021
"""
import pandas as pd
import numpy as np
from scipy.stats import norm
import seaborn as sns
import matplotlib.pyplot as plt

#read data
data = pd.read_excel('data.xlsx', sheet_name = 'data2021M1', header = 0)
X1 = pd.DataFrame.to_numpy(data['X1'])[0:49]
X2 = pd.DataFrame.to_numpy(data['X2'])

n = np.array ([len(X1), len(X2)])
X = np.append (X1, X2)
X1 = ( X1 - np.mean(X) ) / np.std(X)
X2 = ( X2 - np.mean(X) ) / np.std(X)
X = (X - np.mean(X) ) / np.std(X)

#Furbi priors full borrowing
results = np.load("Furbifull.npz")
#results = np.load("GMACRMrealstandard.npz")
theta1save = results['theta1save']
theta2save  = results['theta2save']
tau1save = results['tau1save']
tau2save = results['tau2save']

tot_iter = np.shape(theta1save)[1]
#burnin = int(np.shape(theta1save)[1] / 2)
burnin = 1
temp = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan

for _iter in range(burnin, tot_iter):
    temp[_iter-burnin,0:n[0]] = norm.pdf(X1, loc = theta1save[:,_iter], scale = tau1save[:,_iter]**(-1/2))
    temp[_iter-burnin,n[0]:n[0]+n[1]] = norm.pdf(X2, loc = theta2save[:,_iter], \
                                            scale = tau2save[:,_iter]**(-1/2))

#temp = temp[0::8,:]
prova = temp ** (-1)
CPOFurbi = np.mean(prova, axis=0) ** (-1)
ALCPOFurbi = np.mean(np.log(CPOFurbi))
MLCPOFurbi = np.median(np.log(CPOFurbi))
print('ALCPO FuRBI:', ALCPOFurbi, 'MLCPO Furbi', MLCPOFurbi )

#Furbi priors POS
results = np.load("Furbipos.npz")
theta1save = results['theta1save']
theta2save  = results['theta2save']
tau1save = results['tau1save']
tau2save = results['tau2save']

tot_iter = np.shape(theta1save)[1]
burnin = int(np.shape(theta1save)[1] / 2)
temp = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan

for _iter in range(burnin, tot_iter):
    temp[_iter-burnin,0:n[0]] = norm.pdf(X1, loc = theta1save[:,_iter], scale = tau1save[:,_iter]**(-1/2))
    temp[_iter-burnin,n[0]:n[0]+n[1]] = norm.pdf(X2, loc = theta2save[:,_iter], \
                                            scale = tau2save[:,_iter]**(-1/2))

#temp = temp[0::8,:]
prova = temp ** (-1)
CPOFurbipos = np.mean(prova, axis=0) ** (-1)
ALCPOFurbipos = np.mean(np.log(CPOFurbipos))
MLCPOFurbipos = np.median(np.log(CPOFurbipos))
print('ALCPO FuRBI:', ALCPOFurbipos, 'MLCPO Furbi', MLCPOFurbipos )

#Furbi priors NEG
results = np.load("furbineg.npz")
theta1save = results['theta1save']
theta2save  = results['theta2save']
tau1save = results['tau1save']
tau2save = results['tau2save']

tot_iter = np.shape(theta1save)[1]
burnin = int(np.shape(theta1save)[1] / 2)
temp = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan

for _iter in range(burnin, tot_iter):
    temp[_iter-burnin,0:n[0]] = norm.pdf(X1, loc = theta1save[:,_iter], scale = tau1save[:,_iter]**(-1/2))
    temp[_iter-burnin,n[0]:n[0]+n[1]] = norm.pdf(X2, loc = theta2save[:,_iter], \
                                            scale = tau2save[:,_iter]**(-1/2))

#temp = temp[0::8,:]
prova = temp ** (-1)
CPOFurbineg = np.mean(prova, axis=0) ** (-1)
ALCPOFurbineg = np.mean(np.log(CPOFurbineg))
MLCPOFurbineg = np.median(np.log(CPOFurbineg))
print('ALCPO FuRBI:', ALCPOFurbineg, 'MLCPO Furbi', MLCPOFurbineg )

#GM-dependent 
results = np.load("GMdep.npz")
theta1save = results['theta1save']
theta2save  = results['theta2save']
tau1save = results['tau1save']
tau2save = results['tau2save']

tot_iter = np.shape(theta1save)[1]
burnin = int(np.shape(theta1save)[1] / 2)
temp = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan

for _iter in range(burnin, tot_iter):
    temp[_iter-burnin,0:n[0]] = norm.pdf(X1, loc = theta1save[:,_iter], scale = tau1save[:,_iter]**(-1/2))
    temp[_iter-burnin,n[0]:n[0]+n[1]] = norm.pdf(X2, loc = theta2save[:,_iter], \
                                            scale = tau2save[:,_iter]**(-1/2))

#temp = temp[0::8,:]
prova = temp ** (-1)
CPOFurbiGM = np.mean(prova, axis=0) ** (-1)
ALCPOFurbiGM = np.mean(np.log(CPOFurbiGM))
MLCPOFurbiGM = np.median(np.log(CPOFurbiGM))
print('ALCPO FuRBI:', ALCPOFurbiGM, 'MLCPO Furbi', MLCPOFurbiGM)



#Exch

tot_iter = 10000
burnin = 5000
tempex = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan
results = np.load("EXCHrealstandard.npz")
thetasave = results['thetasave']
sigmasave = results['sigmasave']
csave = results['csave']
nksave = results['nksave']
labelsave = results['labelsave']

temp2 = np.ones((int(n[0]+n[1]),10000)) * np.nan
for i in range(n[0]+ n[1]):
    prova = (np.where(labelsave == csave[i,:]))[0]
    provaiter = (np.where(labelsave == csave[i,:]))[1]
    temp2[i, provaiter] = prova

#theta12save = np.ones((104,10000))
#for _iter in range(0, tot_iter):
 #   theta12save[:,_iter] = thetasave[temp2[:,_iter].astype(int),_iter]

for _iter in range(burnin, tot_iter):
    tempex[_iter-burnin,:] = norm.pdf(X, \
                loc = thetasave[temp2[:,_iter].astype(int),_iter],\
                scale = sigmasave[temp2[:,_iter].astype(int),_iter]**(1/2))

CPOexch = np.mean(tempex ** (-1), axis=0) ** (-1)
ALCPOexch = np.mean(np.log(CPOexch))
MLCPOexch = np.median(np.log(CPOexch))
print('ALCPO exch:', ALCPOexch, 'MLCPO exch:', MLCPOexch)


#Ind

tot_iter = 10000
burnin = 5000
tempex2 = np.ones( (tot_iter-burnin, int(np.sum(n)) ) ) * np.nan
results = np.load("INDrealX1standard.npz")
thetasave = results['thetasave']
sigmasave = results['sigmasave']
csave = results['csave']
nksave = results['nksave']
labelsave = results['labelsave']

temp2 = np.ones((int(n[0]),10000)) * np.nan
for i in range(n[0]):
    prova = (np.where(labelsave == csave[i,:]))[0]
    provaiter = (np.where(labelsave == csave[i,:]))[1]
    temp2[i, provaiter] = prova
    
for _iter in range(burnin, tot_iter):
    tempex2[_iter-burnin,0:n[0]] = norm.pdf(X1, \
                loc = thetasave[temp2[:,_iter].astype(int),_iter],\
                scale = sigmasave[temp2[:,_iter].astype(int),_iter]**(1/2))
        
results = np.load("INDrealX2standard.npz")
thetasave = results['thetasave']
sigmasave = results['sigmasave']
csave = results['csave']
nksave = results['nksave']
labelsave = results['labelsave']

temp2 = np.ones((int(n[1]),10000)) * np.nan
for i in range(n[1]):
    prova = (np.where(labelsave == csave[i,:]))[0]
    provaiter = (np.where(labelsave == csave[i,:]))[1]
    temp2[i, provaiter] = prova
    
for _iter in range(burnin, tot_iter):
    tempex2[_iter-burnin,n[0]:n[0]+n[1]] = norm.pdf(X2, \
                loc = thetasave[temp2[:,_iter].astype(int),_iter],\
                scale = sigmasave[temp2[:,_iter].astype(int),_iter]**(1/2))
        
CPOind = np.mean(tempex2 ** (-1), axis=0) ** (-1)
ALCPOind = np.mean(np.log(CPOind))
MLCPOind = np.median(np.log(CPOind))
print('ALCPO ind:', ALCPOind, 'MLCPO ind:', MLCPOind)

model = [] 
for i in range(n[0]+n[1]):
    model.append("FuRBI \n full")
for i in range(n[0]+n[1]):
    model.append("FuRBI \n -0.95 ")
for i in range(n[0]+n[1]):
    model.append("FuRBI \n 0.95 ")
for i in range(n[0]+n[1]):
    model.append("Exch")    
for i in range(n[0]+n[1]):
    model.append("GM-dep")    
for i in range(n[0]+n[1]):
    model.append("Ind")

data = {'CPO': np.concatenate(((CPOFurbi),(CPOFurbineg), (CPOFurbipos),\
                          (CPOexch), (CPOFurbiGM), (CPOind))),
        'model':model}
CPO = pd.DataFrame(data)

#data.plot.box(grid='True')

sns.set(rc = {'figure.figsize':(5,7)})
sns.set_theme(style="whitegrid")
sns.set_theme(style="ticks", palette="pastel")

ax = sns.boxplot(x = "model", y = "CPO", palette = sns.color_palette("Set2"), data = CPO)
ax.set(xlabel=None)

ax.savefig('boxplot', dpi=300)
ax.savefig('boxplot', dpi=300)
