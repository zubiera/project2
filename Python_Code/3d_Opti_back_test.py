#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:13:43 2019

This File will use the files produced from solving the PDE. and run a back test

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
f_b = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_implicit_Qchange_gamma001_opti.txt", "r")
if f_b.mode == 'r':
    df_b =f_b.read()
    
f_a = open("/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_a_implicit_Qchange_gamma001_opti.txt", "r")
if f_a.mode == 'r':
    df_a =f_a.read()


T = 1382      #Time points 
N = 501      #Variance points
Q = 21       #inventory
nu_max = 200 #Maximum nu used in finite diff grid
q_max  = 10
df_b = create_3d_df(df_b,T,N,Q)
df_a = create_3d_df(df_a,T,N,Q)
##df_a = df_a*0+25
##df_b = df_b*0+25


##Creating a dummy dataframe to test the backtest
#
#folder       = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/vary_vol/gamm001/"
#
#
#ask_deltas = "vol_sigma_10d_ave_Opti_delta_a.csv"
#bid_deltas = "vol_sigma_10d_ave_Opti_delta_b.csv"
#
##turn the data sheets into dataframes
#dfa    =   pd.read_csv(folder + ask_deltas ) #,header=None
#dfb    =   pd.read_csv(folder + bid_deltas ) #,header=Non
#
#dfa= dfa.drop(["Time_Step"], axis=1)
#dfb= dfb.drop(["Time_Step"], axis=1)
#dfa = dfa.append(dfa.loc[[len(dfa)-1]],ignore_index=True)
#dfb = dfb.append(dfb.loc[[len(dfb)-1]],ignore_index=True)
#
#dfa = np.array(dfa)
#dfb = np.array(dfb)
#dfa = dfa[::-1]
#dfb = dfb[::-1]
#
#
#T=1382
#N=501
#Q=21
#
#df_a  = np.zeros((T,N,Q)) #3d array set
#df_b  = np.zeros((T,N,Q)) #3d array set
#
##Just making a 3d array with the 2d repeated
#for i in range(0,N):
#    df_a[:,i,:] = dfa
#    df_b[:,i,:] = dfb


#df_a  = np.zeros((T,N,Q))
#counter  = 0
#for t in range(T):
#    for j in range(N):
#        for q in range(Q):
#            df_a[t][j][q] =  24.0 #purposely 24, to check it rounds off later to 25.
#            counter+=1 
#df_b  = np.zeros((T,N,Q)) #3d array set
#
#counter  = 0
#for t in range(T):
#    for j in range(N):
#        for q in range(Q):
#            df_b[t][j][q] =  24.0 #purposely 24, to check it rounds off later to 25.
#            counter+=1 
#

# Load the PAst data to back test on
data_file_list = ['17','18', '19','20','21','24','25','26','27','28','31']

Max_inv=0    #P0aramteer to check max_inv 

PnL  = np.array([['No_Trades','Cash',"maxinv","inv_at_T"]])


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
    event_tracker  = np.array([['minute','time_idx','event','inventory','cash','ex_price','delta','nu']])
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
        
        
        #*************************  Iterate Within each Timestep    **********************
        for j in range( np.searchsorted(time,dt*i), np.searchsorted(time,dt*(i+1)) ):
            
        #************************below checks for ask side executions****************************
            #check if delta implies a market order  (if so, insntant execution)
            if (delta_a < 0) and (is_a_hit == False):
                x = x + bid_price[j]  #execute as market order at best bid
                q = q - 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'MO_sell', int(q), float(x), bid_price[j],delta_a,int(nu )]])))
                is_a_hit = True
            
            #considers the corner case where Quote is priced at best ask
            elif delta_a == 0 and (is_a_hit == False):
                if ( (ask_price[j] == ask_price[j+1] == LO_a) and (ask_vol[j]>ask_vol[j+1]) ) \
                or (ask_price[j] >  LO_a):
                    x = x + LO_a             # Update cashflow by adding sale revenue               
                    q = q - 1                 # Update inventory Level 
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'a at 0', int(q), float(x), LO_a, delta_a , int(nu )]])))
                    is_a_hit = True          
            
            #This is the case where ref price becomes equal to the quote
            elif LO_a == ask_price[j+1] and (is_a_hit == False):
                #subcase #The case an MO walks 1 level trhough the book, eating some 
                #volume where our quote is
                if (ask_price[j+1]==LO_a) and (ask_price[j] < LO_a) and (ask_vol[j+1]<ask_vol2[j]):
                    x = x + LO_a
                    q = q - 1
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'a at MO eating 1+ level', int(q), float(x), LO_a,delta_a, int(nu ) ]])))
                    is_a_hit = True
                #subcase where our quote becomes equal to best ask, but no volume
                #has been been consumed just yet.
                elif ask_price[j] == ask_price[j+1] and ask_vol[j] > ask_vol[j+1]:
                    x = x + LO_a
                    q = q - 1
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'a at MO eating our quote', int(q), float(x), LO_a, delta_a, int(nu ) ]])))
                    is_a_hit = True
            
            #The case where a large MO results in the ref price becoming larger
            #than the quote (i.e it consumes our quote while walking the book
            elif LO_a <  ask_price[j+1] and (is_a_hit == False):
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
                
            #considers the corner case where Quote is priced at best bid
            elif delta_b == 0 and (is_b_hit == False):
                if ( (bid_price[j] == bid_price[j+1]) and (bid_vol[j]>bid_vol[j+1]) ) \
                or (bid_price[j] <  LO_b):
                    x = x - LO_b             
                    q = q + 1
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'b at 0', int(q), float(x), LO_b, delta_b, int(nu )]])))
                    is_b_hit = True
                    
            #This is the case where ref price becomes equal to the quote
            elif LO_b == bid_price[j+1] and (is_b_hit == False):
                    
                #subcase #The case an MO walks atleast 1 level through the book, eating some 
                #volume where our quote is
                if (bid_price[j+1]==LO_b) and (bid_price[j] > LO_b) and (bid_vol[j+1]<bid_vol2[j]):
                    x = x - LO_b
                    q = q + 1
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'b at MO eating 1+ level', int(q), float(x), LO_b, delta_b, int(nu ) ]])))
                    is_b_hit = True
                #subcase where our quote becomes equal to best bid, but no volume
                #has been been consumed just yet.
                elif bid_price[j] == bid_price[j+1] and bid_vol[j] > bid_vol[j+1]:
                    x = x - LO_b
                    q = q + 1
                    event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'b at MO eating our quote', int(q), float(x), LO_b, delta_b, int(nu ) ]])))
                    is_b_hit = True
            
            #The case where a large MO results in the  Best-Bid Price becoming smaller
            #than the quote (i.e it consumes our quote while walking the book
            elif LO_b >  bid_price[j+1] and (is_b_hit == False):
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
    Max_inv =  max(abs(min(event_tracker["inventory"])), abs(max(event_tracker["inventory"])) ) 
    
    PnL = np.concatenate((PnL, np.array([[len(event_tracker),cash,Max_inv, inv_at_T]])))
#print(bid_vol[-1]," ", ask_vol[-1])




        
        





 

    
    