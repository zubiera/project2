#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 15:15:18 2019

@author: amitarfan1
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

from data_load import create_df               #funtion i made to convert dataframe
from heston_calibration import  Heston_cali   #function i made to calibrate

days = ['03','04','05','06','07','10','11','12','13','14','17','18', '19','20','21','24','25','26','27','28','31']
days = ['18']

Volumes   = pd.DataFrame()
Variances = pd.DataFrame()

for i in range(len(days)):
    
    df = create_df('data/1707'+ days[i] +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
                   
    
                   
    # the below is to get the data in minutley format
    time      = np.array(df["time_from_open"])
    ask_price = np.array(df["a1"])
    ask_vol   = df["av1"] + df["av2"] + df["av3"] + df["av4"] + df["av5"]
    bid_vol   = df["bv1"] + df["bv2"] + df["bv3"] + df["bv4"] + df["bv5"]
    total_vol = ask_vol + bid_vol
    dt=6e10   
    S_a  = np.array(ask_price[0])  
    vol   = np.array(total_vol[0])
    minutes = np.array([0])   
    
    
    i=0
    while dt*i <  time[-1]:
        
        S_a =   np.append (S_a, [ask_price[np.searchsorted(time,dt*i)]]) 
        vol  =   np.append (vol,  [ask_vol[np.searchsorted(time,dt*i)]])
        minutes =   np.append (minutes, i)
        i+=1
    S_a  = np.delete(S_a, 0)  #Remove first dummy eleemnt
    vol  = np.delete(vol,0)    #Remove first dummy element
    minutes = np.delete(minutes, 0)  #Remove first dummy element
    
    
    S_a = pd.Series(S_a)    #make into a panda time series
    vol  = pd.Series(vol)
    minutes = pd.Series(minutes)
    
    #Pandas series to use
    ask        = S_a
    RVask60    = S_a.rolling(60).var()
    #Vk = RVask60.apply(lambda x: np.where(x > RVask60[:].quantile(0.9),RVask60[:].quantile(0.9),x))
    Vk = Vk[:]
    volume = vol
    Volumes = pd.concat([Volumes, volume], axis=1)
    
    Variances = pd.concat([Variances,RVask60],axis=1)
    
#Variances.columns = ['03','04','05','06','07','10','11','12','13','14','17','18', '19','20','21','24','25','26','27','28','31']
#Variances['03'].plot()
#Variances['04'].plot()
#Variances['05'].plot()
#Variances['06'].plot()
#Variances['07'].plot()

    