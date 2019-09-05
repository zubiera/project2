#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:55:03 2019

@author: amitarfan1
"""

import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Axes3D import has side effects, it enables using projection='3d' in add_subplot
import matplotlib.pyplot as plt
import random
import pandas as pd


#Have to find the location of executable file and run that
theta = "0.02"
alpha = "0.09"
xi    = "0.01"
file   = "/Users/amitarfan1/Documents/Computational\ Finance/Xcode\ Output/Build/Products/Debug/Finite_difference "
params = theta + " " + alpha + " " + xi  
#Run the C++ file  and read thre cout line as a big string
deltas = os.popen(file+params).readlines()
os.system(file + params)

#******THe belos is not needed right now because I have saved the output from the C++ file, no need to read lines*******
#f_b = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_str.txt","w+")
#f_b.write(deltas[0])
#f_b.close()



#**************verify it works


"""
THe below is jsut to see if our storage of the strings is not being corrupted and is being read correctly
"""
#
#f=open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_str.txt", "r")
#if f.mode == 'r':
#    delta_b =f.read()
#    
#delta_b = delta_b.split(",")            #Split the string by comma's
#delta_b = list(map(float,delta_b[:-1]))    #turns the strings into floats and removes the last space
#delta_b = np.array(delta_b)
#
#
#T = 301
#N = 201
#Q = 21
#
##set up a 3d array to store the valeus for each combo
#values  = np.zeros((T,N,Q))  
#
#
#counter  = 0
#for t in range(T):
#    for j in range(N):
#        for q in range(Q):
#            values[t][j][q] =  delta_b[counter]
#            counter+=1
#
##***************************************Plot 3d Grapgs*********************************
#import matplotlib.patches as mpatches
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#X, Y = np.meshgrid(np.arange(N),np.arange(T))
#zs1 = np.array(list(reversed(values[:,:,1])))
##beware, something funny happens with axis, therefore the mesgrid has to be done wrongway around
#Z1 = zs1.reshape(X.shape)
#
#for i in range(4,17,2):
#    zs = np.array(list(reversed(values[:,:,i])))
#    Z  = zs.reshape(X.shape)
#    ax.plot_surface(X, Y, Z)