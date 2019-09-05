#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is to calibrate the parameters for the heston model using daily data for 
ten diffent trading days.
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

from data_load import create_df               #funtion i made to convert dataframe
from heston_calibration import  Heston_cali   #function i made to calibrate


days = ['03','04','05','06','07','10','11','12','13','14','17' ,'18', '19','20','21','24','25','26','27','28','31']
days = ['07']
Parameters = pd.DataFrame(columns=['theta', 'xi', 'alpha','rho'])

for i in range(len(days)):
    
    df = create_df('data/1707'+ days[i] +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
                   
    
                   
    # the below is to get the data in minutley format
    time      = np.array(df["time_from_open"])
    ask_price = np.array(df["a1"])
    ask_vol   = df["av1"] + df["av2"] + df["av3"] + df["av4"] + df["av5"]
    dt=6e10   
    S_a  = np.array(ask_price[0])  
    av   = np.array(ask_vol[0])
    minutes = np.array([0])   
    i=0
    while dt*i <  time[-1]:
        
        S_a =   np.append (S_a, [ask_price[np.searchsorted(time,dt*i)]]) 
        av  =   np.append (av,  [ask_vol[np.searchsorted(time,dt*i)]])
        minutes =   np.append (minutes, i)
        i+=1
    S_a = np.delete(S_a, 0)  #Remove first dummy eleemnt
    av  = np.delete(av,0)    #Remove first dummy element
    minutes = np.delete(minutes, 0)  #Remove first dummy element
    
    S_a = pd.Series(S_a)    #make into a panda time series
    av  = pd.Series(av)
    minutes = pd.Series(minutes)
    
    #Pandas series to use
    ask        = S_a[:]
    ask        = ask[:]
    RVask60    = S_a[:].rolling(60).var()
    RVask60    = RVask60[:]
    
    theta, xi, alpha, rho = Heston_cali(ask,RVask60,1)
    Parameters = Parameters.append(pd.Series([theta,xi,alpha,rho],index = Parameters.columns),ignore_index=True)
    
    