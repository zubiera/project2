#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 14:20:15 2019

@author: amitarfan1
"""

import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Axes3D import has side effects, it enables using projection='3d' in add_subplot
import matplotlib.pyplot as plt
import random
import pandas as pd

from data_load import create_3d_df #Converst pde soln str to df
from data_load import create_df



#Have to find the location of executable file and run that
theta = "0.0006"  #"0.02"
alpha = "6948.75" #"0.09"
xi    = "4.2365"  #"0.01"
rho   = "0.5"#"0.0037"  #"0.7"
gamma = "0.01"   #0.05
A     = "1.1956"  #0.1 
k     = "0.0435"  #0.3
T      = "1380"
nu_max = "12000"


#     ********** Evening paramters avaraged for 4-7 july 
theta = "0.0015"
alpha = "11459"
xi    = "6.6095"
rho   = "-0.0172"
gamma = "0.01"   #0.05
A     = "1.1956"  #0.1 
k     = "0.0435"  #0.3
T      = "1380"
nu_max = "20000"


#####Ideallist parameters ###
#theta = "0.02"
#alpha = "0.09"
#xi    = "0.01"
#rho   = "0.7"
#gamma = "0.05"
#A     = "0.1" 
#k     = "0.3"
#T      = "300"
#nu_max = "100"

file   = "/Users/amitarfan1/Documents/Computational\ Finance/Xcode\ Output/Build/Products/Debug/Finite_difference "
params = theta + " " + alpha + " " + xi  + " " + rho + " " + gamma + " " + A + " " + k + " " + T  + " " + nu_max   

#test line, delete when done
os.system(file + params)


f_b = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_str_cali.txt", "r")
if f_b.mode == 'r':
    delt_b =f_b.read()

#delta = delt_b.split(",")            #Split the string by comma's
#delta = list(map(float,delta[:-1]))    #turns the strings into floats and removes the last space
#delta = np.array(delta)

T_steps = 1382 #1381      #Time points 
N       = 501      #Variance points
Q       = 21       #inventory


values = create_3d_df(delt_b,T_steps,N,Q)




#************************************* Plots    ******************************************
import matplotlib.patches as mpatches
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(np.arange(0,12001,12000/500),np.arange(T_steps))


for i in range(4,17,2):
    zs = np.array(list(reversed(values[:,:,i])))
    Z  = zs.reshape(X.shape)
    ax.plot_surface(X, Y, Z)
    
ax.set_xlabel('Variance')
ax.set_ylabel('Time')
ax.set_zlabel('$\delta^b$')
ax.set_title("evening")
plt.show()




##Delete the belwo when done.
#delta_vs_time = list(reversed(values[:,100,:])) 
#delta_vs_time = pd.DataFrame(delta_vs_time)  #make into df for ease of plotting
#time  = np.arange(T_steps)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#for i in range(0,20,2):
#    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='blue')

#delta_vs_time = list(reversed(values[:,10,:])) 
#delta_vs_time = pd.DataFrame(delta_vs_time)  #make into df for ease of plotting
#time  = np.arange(T_steps)
#
#for i in range(0,20,2):
#    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='green')
#    
#
#delta_vs_time = list(reversed(values[:,100,:])) 
#delta_vs_time = pd.DataFrame(delta_vs_time)  #make into df for ease of plotting
#time  = np.arange(T_steps)
#
#for i in range(0,20,2):
#    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='red')
#    
#delta_vs_time = list(reversed(values[:,200,:])) 
#delta_vs_time = pd.DataFrame(delta_vs_time)  #make into df for ease of plotting
#time  = np.arange(T_steps)
#
#for i in range(0,20,2):
#    plt.plot(time,delta_vs_time[i],linewidth = 0.5,color='black')