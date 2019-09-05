#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:13:43 2019

This File will use the files produced from solving the PDE. and run a  pessi back test

"""
import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Axes3D import has side effects, it enables using projection='3d' in add_subplot
import matplotlib.pyplot as plt
import random
import pandas as pd

from data_load import create_3d_df #Converst pde soln str to df
from data_load import create_df    #Makes a df out of the trading data




#********************** Open The Files Needed and convert into a 3d dataframe ***********************
f_b = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_implicit_Qchange_gamma001_pessi.txt", "r")
if f_b.mode == 'r':
    df_b =f_b.read()
    
f_a = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_a_implicit_Qchange_gamma001_pessi.txt", "r")
if f_a.mode == 'r':
    df_a =f_a.read()


T = 1382      #Time points 
N = 501      #Variance points
Q = 21       #inventory
nu_max = 200 #Maximum nu used in finite diff grid
q_max  = 10
df_b = create_3d_df(df_b,T,N,Q)
df_a = create_3d_df(df_a,T,N,Q)


# Load the PAst data to back test on
data_file_list = ['17','18', '19','20','21','24','25','26','27','28','31']

Max_inv=0    #P0aramteer to check max_inv 

PnL  = np.array([['No_Trades','Cash',"maxinv", "inv_at_T"]])


nu_max = 20000
q_max  = 10
dt = 60e9
#Loop through each backtest
for  file in  data_file_list:
    #**************************Process of the raw data**************************
    df = create_df('data/1707'+str(file)+'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
    #df = df[:100000]
    
                
    time      = np.array(df["time_from_open"])  
    ask_price = np.array(df["a1"])
    ask_vol   = np.array(df["av1"])
    ask_vol2   = np.array(df["av2"])
    
    bid_price = np.array(df["b1"])
    bid_vol   = np.array(df["bv1"])
    bid_vol2   = np.array(df["bv2"])
    
    #Getting the price per minute, and then calcing rolling vriance
    S_a = ask_price[0]
    S_b = bid_price[0]
    k=0
    while dt*k <  time[-1]:
        
        S_a =   np.append (S_a, [ask_price[np.searchsorted(time,dt*k)]])
        S_b =   np.append (S_b, [bid_price[np.searchsorted(time,dt*k)]])
        k+=1
    S_a = np.delete(S_a, 0)
    S_b = np.delete(S_b, 0) 
    RV_ask_60 = np.array(pd.DataFrame(S_a).rolling(60).var())
    RV_bid_60 = np.array(pd.DataFrame(S_b).rolling(60).var())
    
    RV_ask_60[:59]=RV_ask_60[59]
    RV_bid_60[:59]=RV_bid_60[59]
    
    

   #overwriting with const volatility 
    #RV_ask_60 = np.ones(1382)*50
    #RV_bid_60 = np.ones(1382)*50
    
    
    #*******************function to round the quotes to nearest tick***************
    
    
    # This code rounds to nearest25
    df_b = np.around(df_b/25, decimals=0)*25
    df_a = np.around(df_a/25, decimals=0)*25
    
    #add something here to set qutes realy far for q_max and -q_max
    df_a[:,:,0]   = 500  #sts ask soo low at -q_max so we dont sell more.
    df_b[:,:,20]  = 500  #sets bid so high that we cant hold more than q_max

    #***************************** Run Strategy **************************************
    
      #The algo posts quotes every dt.  (60billion nano-seconds = 1 minute)
    
    #initlisise paramters
    q   = 0      #Initial inventory 
    x   = 0      #Initial cashflow
    S_a = [0]    #vector tracking quoted ask price at each dt
    S_b = [0]
    idx = [0]    #vector tracking index locations
    event_tracker  = np.array([['minute','time_idx','event','inventory','cash','ex_price','delta', 'nu']])
    quote_tracker  = np.array([['time_idx','ask_quote','bid_quote']])
    
    
    i=0
    while np.searchsorted(time,dt*(i+1))< len(time):
        
        # can manually overwride this to determin fixed strategy vs DataFrame strategy
        nu = RV_ask_60[i]
        #nu = RV_ask_60[np.searchsorted(time,dt*(i+1))]
        nu_idx =min(int( nu/nu_max * N),N-1)   #index of nu in 3d data fram. #min to cap very large nu going out of array.
        delta_a =  df_a[1381-i,nu_idx,q_max + q]
        
        nu = RV_bid_60[i]
        #nu = RV_bid_60[np.searchsorted(time,dt*(i+1))]
        nu_idx =min(int( nu/nu_max * N),N-1)  #index of nu in 3d data fram.
        delta_b =  df_b[1381-i,nu_idx,q_max + q]
        
        
                
        #setting the Limit order quotes
        LO_a = ask_price[np.searchsorted(time,dt*i)] + delta_a
        LO_b = bid_price[np.searchsorted(time,dt*i)] - delta_b
        

        #keeping track of what is quoted at what time.    
        quote_tracker = np.concatenate((quote_tracker, np.array([[np.searchsorted(time,dt*i), LO_a, LO_b]])))
        
        is_b_hit = False
        is_a_hit = False
      
                
        
        #iterate within each timestep
        for j in range( np.searchsorted(time,dt*i), np.searchsorted(time,dt*(i+1)) ):
            
        #************************below checks for ask side executions****************************
            #check if delta implies a market order  (if so, insntant execution)
            if (delta_a < 0) and (is_a_hit == False):
                x = x + bid_price[j]  #execute as market order at best bid
                q = q - 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'MO_sell', int(q), float(x), bid_price[j],delta_a,int(nu ) ]])))
                is_a_hit = True
    

            if LO_a <  ask_price[j+1] and (is_a_hit == False):
                x = x + LO_a
                q = q - 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'a at MO walking over our price', int(q), float(x), LO_a, delta_a, int(nu ) ]])))
                is_a_hit = True
                    
                    

                
        #*************Below checks for  bid sides executions****************************************
            #check if delta implies a market order
            if (delta_b < 0) and (is_b_hit == False):
                x = x - ask_price[j]  #execute as market order at best bid
                q = q + 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'MO_buy', int(q), float(x), ask_price[j], delta_b, int(nu ) ]])))
                is_b_hit = True
                
            if LO_b >  bid_price[j+1] and (is_b_hit == False):
                x = x - LO_b
                q = q + 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'b at MO walking over our price', int(q), float(x), LO_b, delta_b, int(nu ) ]])))
                is_b_hit = True
                    

            #break the loop if both ask and bid have been executed in this dt
            if (is_b_hit == True) and (is_a_hit == True):
                break
        i+=1
    
    
    
    #*************** Liquidate any inventory at time T *****************
    if float(event_tracker[-1][3]) < 0:
        cash = float(event_tracker[-1][4]) - abs(float(event_tracker[-1][3]))*ask_price[-1]
    
    if float(event_tracker[-1][3]) > 0:
        cash = float(event_tracker[-1][4]) + abs(float(event_tracker[-1][3]))*bid_price[-1]
    
    if float(event_tracker[-1][3]) == 0:
        cash = float(event_tracker[-1][4])
        
    inv_at_T = event_tracker[-1][3]

    #convnert event_tracker to a dataframe
    event_tracker = pd.DataFrame(event_tracker)
    event_tracker.columns = event_tracker.iloc[0]  #make row zero column headings
    event_tracker = event_tracker.reindex(event_tracker.index.drop(0))  #drop row zero
    event_tracker['inventory'] = event_tracker['inventory'].astype('float64') 
    Max_inv = max(Max_inv, max(abs(min(event_tracker["inventory"])), abs(max(event_tracker["inventory"])) ) )
    
    PnL = np.concatenate((PnL, np.array([[len(event_tracker),cash, Max_inv,inv_at_T]])))







