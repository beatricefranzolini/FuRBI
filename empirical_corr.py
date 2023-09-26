# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 10:33:35 2021
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

#set the working directory to source file location
mypath = Path().absolute()
os.chdir(mypath)

np.random.seed(0) #set seed for replicability 

#read data
df = pd.read_excel('data.xlsx', sheet_name = 'avg_returns', header=0, \
                   index_col=0, parse_dates=True, squeeze=True)
X1 = pd.DataFrame.to_numpy(data['stocks'])
X2 = pd.DataFrame.to_numpy(data['comm'])

T = len(X1)
corr = np.ones(T) * np.nan
for t in range(6, T-6):
    corr[t] = np.corrcoef(X1[t-6:t+6],X2[t-6:t+6])[0,1]

df['corr'] = corr

df['corr'].plot(color= 'blue', linewidth = 2, label="Empirical correlation")
plt.ylim((-1,1))
plt.legend()
plt.savefig('empcorr', dpi=300)
plt.show()

