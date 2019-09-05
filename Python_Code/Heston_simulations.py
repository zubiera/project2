#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is to simulate the heston model
"""
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import random

#np.random.seed(0)
		
S0    =  241025.0   # 100      
V0    =  1281.0734463283723
 #4        
theta =  0.0015#50 #0.02  
alpha =  12134.44  #4    
eta   =  6.8968#0.5    
rho   =  0.0938 #0.7   

dt    =     1 #0.005

V = np.zeros(1381-800)
S = np.zeros(1381-800)
V[0] = V0
S[0] = S0

# data fram to store all the simulations
S_df = pd.DataFrame()
V_df = pd.DataFrame()

for j in range(100):
    for i in range(1,len(V)):
        
        Zv   =  np.random.normal(0, 1)
        Zs   =  rho*Zv   + np.sqrt(1-rho*rho)*np.random.normal(0, 1)
        #V_plus   = max(0,V[i-1])
        V[i] = V[i-1] + theta*(alpha - max(0,V[i-1]))*dt + eta*np.sqrt(max(0,V[i-1]))*Zv*np.sqrt(dt)
        S[i] = S[i-1] + np.sqrt(max(0,V[i-1]))*Zs*np.sqrt(dt)
    S_df[j]  = S 
    V_df[j]  = V 
        
plt.figure()
plt.plot(V_df,linewidth=0.1)
plt.plot(np.array(RVask60[800:]),color='k')
plt.title("Heston Model Variance vs Rolling Variance")
plt.xlabel("time")
plt.ylabel("variance")
plt.tight_layout()

plt.figure()
plt.plot(S_df,linewidth=0.1)
plt.plot(np.array(S_a[800:]),color='k')
plt.title("Heston Model Price vs Real Price")
plt.xlabel("time")
plt.ylabel("price")
plt.tight_layout()


#************************************* K-S Test **************************************
from scipy.stats import kstest
from scipy.stats import ks_2samp
tester = np.random.normal(0,1,1000)
kstest(tester, 'norm')


Real_data = S_a[800:]
for i in range (S_df.shape[1]):
    sim_data  = S_df[i]
    print(ks_2samp(S_df[i],S_df[i]))



#****************************** BackOut Initial Parameters from Simulations ****************

from heston_calibration import  Heston_cali

parameter_df = pd.DataFrame(columns=['theta', 'xi', 'alpha','rho'])
#parameter_df.append(pd.Series([1,2,3,4],index = parameter_df.columns),ignore_index=True)

for i in range(np.shape(V_df)[1]):
    ask = S_df[i]
    var = V_df[i]
    theta, xi, alpha, rho = Heston_cali(ask,var,dt)
    parameter_df= parameter_df.append(pd.Series([theta,xi,alpha,rho],index = parameter_df.columns),ignore_index=True)

import seaborn as sns
import pylab 
import scipy.stats as stats





#******************************  Dist plots of the Parameters    ************************************
plt.subplot(2,2,1)
plt.title('KDE Plot Of Estimated theta')
plt.grid()
sns.distplot(parameter_df['theta'])

plt.subplot(2,2,2)
plt.title('KDE Plot Of Estimated xi ')
plt.xlabel('xi')
plt.grid()
sns.distplot(parameter_df['xi'])

plt.subplot(2,2,3)
plt.title('KDE Plot Of Estimated alpha')
plt.grid()
sns.distplot(parameter_df['alpha'])

plt.subplot(2,2,4)
plt.title('KDE Plot Of Estimated  rho ')
plt.grid()
sns.distplot(parameter_df['rho'])

plt.tight_layout()



#******************************  QQ plots of the Parameters    ***************************
fig = plt.figure()

ax = fig.add_subplot(221)
r = stats.probplot(parameter_df['theta']*0.005/2, dist='lognorm',sparams=(0.5968081), plot=plt,rvalue=True)
ax.get_lines()[0].set_marker('.')
ax.get_lines()[0].set_markeredgecolor('k')
ax.get_lines()[0].set_markerfacecolor('w')
ax.get_lines()[0].set_markersize(5.0)
ax.get_lines()[1].set_linewidth(3.0)
ax.get_lines()[1].set_color('purple')
plt.title("QQ plot for theta")


ax = fig.add_subplot(222)
r = stats.probplot(parameter_df['xi'], dist='norm', plot=plt,rvalue=True)
ax.get_lines()[0].set_marker('.')
ax.get_lines()[0].set_markeredgecolor('k')
ax.get_lines()[0].set_markerfacecolor('w')
ax.get_lines()[0].set_markersize(5.0)
ax.get_lines()[1].set_linewidth(3.0)
ax.get_lines()[1].set_color('purple')
plt.title("QQ plot for xi")

ax = fig.add_subplot(223)
r = stats.probplot(parameter_df['alpha'], dist='norm', plot=plt,rvalue=True)
ax.get_lines()[0].set_marker('.')
ax.get_lines()[0].set_markeredgecolor('k')
ax.get_lines()[0].set_markerfacecolor('w')
ax.get_lines()[0].set_markersize(5.0)
ax.get_lines()[1].set_linewidth(3.0)
ax.get_lines()[1].set_color('purple')
plt.title("QQ plot for alpha")

ax = fig.add_subplot(224)
r = stats.probplot(parameter_df['rho'], dist='norm', plot=plt, rvalue=True)
ax.get_lines()[0].set_marker('.')
ax.get_lines()[0].set_markeredgecolor('k')
ax.get_lines()[0].set_markerfacecolor('w')
ax.get_lines()[0].set_markersize(5.0)
ax.get_lines()[1].set_linewidth(3.0)
ax.get_lines()[1].set_color('purple')
plt.title("QQ plot for rho")
plt.show()
plt.tight_layout()
