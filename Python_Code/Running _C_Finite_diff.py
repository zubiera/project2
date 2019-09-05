#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file runs the c++ file that solves the PDe for the stochastic volatility variateion
Here we aim to jus run the c+++ file and read the output.

"""

import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Axes3D import has side effects, it enables using projection='3d' in add_subplot
import matplotlib.pyplot as plt
import random
import pandas as pd


#Have to find the location of executable file and run that
theta = "0.0006"  #"0.02"
alpha = "6948.75" #"0.09"
xi    = "4.2365"  #"0.01"
rho   = "0.0037"  #"0.7"
gamma = ["0.01","0.01", "0.0001", "0.0001"]  #"0.01"   #0.05
A     = ["1.1956", "0.218", "1.1956", "0.218"]  #0.1 
k     = ["0.0435",  "0.03533", "0.0435",  "0.03533"] #0.3
T      = "1380"
nu_max = "20000"
label  = ["gamma01_opti", "gamma01_pessi", "gamma0001_opti", "gamma0001_pessi"]


for i in range(len(A)):
    file   = "/Users/amitarfan1/Documents/Computational\ Finance/Xcode\ Output/Build/Products/Debug/finite_difference_implicit "
    params = theta + " " + alpha + " " + xi  + " " + rho + " " + gamma[i] + " " + A[i] + " " + k[i] + " " + T  + " " + nu_max  + " "+ label[i] 
    
    #test line, delete when done
    os.system(file + params)




f_b = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_implicit_Qchange.txt", "r")
if f_b.mode == 'r':
    delt_b =f_b.read()

from data_load import create_3d_df #Converst pde soln str to df
T = 1382      #Time points 
N = 501      #Variance points
Q = 21       #inventory
nu_max = 20000 #Maximum nu used in finite diff grid
q_max  = 10
values = create_3d_df(delt_b,T,N,Q)


#Run the C++ file  and read thre cout line as a big string
delta = os.popen(file+params).readlines()
delta = delta[0].split(",")            #Split the string by comma's
delta = list(map(float,delta[:-1]))    #turns the strings into floats and removes the last space
delta = np.array(delta)



T = 1382
N = 501
Q = 21

#set up a 3d array to store the valeus for each combo
values  = np.zeros((T,N,Q))  


counter  = 0
for t in range(T):
    for j in range(N):
        for q in range(Q):
            values[t][j][q] =  delta[counter]
            counter+=1




#plot the delta vs time, make  S constant because it doesnt matter in geants paper
# this is becuse out numercal method solves backwards in time, so we need to rever plot it to see \
#It in the real order
delta_vs_time = list(reversed(values[:,0,:])) 
delta_vs_time = pd.DataFrame(delta_vs_time)  #make into df for ease of plotting
time  = np.arange(T)
fig = plt.figure()
ax = fig.add_subplot(111)
delta_vs_time = list(reversed(values[:,200,:])) 
delta_vs_time = pd.DataFrame(delta_vs_time) 
for i in range(0,20):
    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='blue')
delta_vs_time = list(reversed(values[:,100,:])) 
delta_vs_time = pd.DataFrame(delta_vs_time) 
for i in range(0,20):
    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='red')
delta_vs_time = list(reversed(values[:,50,:])) 
delta_vs_time = pd.DataFrame(delta_vs_time) 
for i in range(0,20):
    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='black')    


import matplotlib.patches as mpatches
blue = mpatches.Patch(color='blue', label='nu=1')  
red = mpatches.Patch(color='red', label='nu=0.5')
black = mpatches.Patch(color='black', label='nu=0.25')
ax.legend(handles=[blue,red, black])
ax.set_xlabel('Time')
ax.set_ylabel('$\delta^b$')



#***************************************Plot 3d Grapgs*********************************
import matplotlib.patches as mpatches
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(np.arange(N),np.arange(T))
zs1 = np.array(list(reversed(values[:,:,1])))
#beware, something funny happens with axis, therefore the mesgrid has to be done wrongway around
Z1 = zs1.reshape(X.shape)

for i in range(4,17,2):
    zs = np.array(list(reversed(values[:,:,i])))
    Z  = zs.reshape(X.shape)
    ax.plot_surface(X, Y, Z)
   
#blue = mpatches.Patch(color='blue', label='q=1')
#orange = mpatches.Patch(color='orange', label='q=2')
#green = mpatches.Patch(color='green', label='q=3')
#red = mpatches.Patch(color='red', label='q=4')
#purp = mpatches.Patch(color='purple', label='q=5')
#brown = mpatches.Patch(color='brown', label='q=6')
#ax.legend(handles=[blue,orange, green, red, purp, brown])

ax.set_xlabel('j')
ax.set_ylabel('i')
ax.set_zlabel('$\delta^b$')
plt.show()







# This is to plot delta vs nu for stoctastic vol as we change theta and xi

theta = ["0.02","0.05","0.1","1.0","2.0","10.0"]
alpha = "0.09"
xi    = ["0.01", "0.001", "0.0001","0.000001","0.0000001","0.0000000001"]
file   = "/Users/amitarfan1/Documents/Computational\ Finance/Xcode\ Output/Build/Products/Debug/Finite_difference "
delta_store_2 = pd.DataFrame()
for i in range(len(theta)):
    params = theta[i] + " " + alpha + " " + xi[i]  
    #Run the C++ file  and read thre cout line as a big string
    delta = os.popen(file+params).readlines()
    delta = delta[0].split(",")            #Split the string by comma's
    delta = list(map(float,delta[:-1]))    #turns the strings into floats and removes the last space
    delta = np.array(delta)
    
    T = 301
    N = 201
    Q = 21
    
    #set up a 3d array to store the valeus for each combo
    values  = np.zeros((T,N,Q))  
    
    
    counter  = 0
    for t in range(T):
        for j in range(N):
            for q in range(Q):
                values[t][j][q] =  delta[counter]
                counter+=1
    
    delta_vs_nu = list((values[:,:,5]))
    delta_vs_nu = pd.DataFrame(delta_vs_nu)
    delta_vs_nu = delta_vs_nu.T
    plt.plot(np.arange(N),delta_vs_nu[150])    
    delta_store_2 = pd.concat([delta_store_2,delta_vs_nu[150]],axis=1)

plt.xlabel('j')
plt.ylabel('$\delta^b$')

