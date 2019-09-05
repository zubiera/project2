import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


from data_load import create_df

## *********** Ths code here is to be able to calibrate the parameters for the heston model
#
#
#
#days = ['03']
#
#df = create_df('data/1707'+ '06' +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
#               
#
#               
## the below is to get the data in minutley format
#time      = np.array(df["time_from_open"])
#ask_price = np.array(df["a1"])
#ask_vol   = df["av1"] + df["av2"] + df["av3"] + df["av4"] + df["av5"]
#dt=6e10   
#S_a  = np.array(ask_price[0])  
#av   = np.array(ask_vol[0])
#minutes = np.array([0])   
#i=0
#while dt*i <  time[-1]:
#    
#    S_a =   np.append (S_a, [ask_price[np.searchsorted(time,dt*i)]]) 
#    av  =   np.append (av,  [ask_vol[np.searchsorted(time,dt*i)]])
#    minutes =   np.append (minutes, i)
#    i+=1
#S_a = np.delete(S_a, 0)  #Remove first dummy eleemnt
#av  = np.delete(av,0)    #Remove first dummy element
#minutes = np.delete(minutes, 0)  #Remove first dummy element
#
#
#S_a = pd.Series(S_a)    #make into a panda time series
#av  = pd.Series(av)
#minutes = pd.Series(minutes)
#
##Pandas series to use
#ask        = S_a
#ask        = ask[:800]
#RVask30    = S_a.rolling(60).var()
#RVask30    = RVask30[:800]
#ask_volume = av
#Vk         = RVask30.dropna()      #rolling volatility used as instantaneous vol
#Vk = Vk.apply(lambda x: np.where(x < 20,20,x))
##Vk         = Vk[Vk > 20]           #********** Dropping anomalies*********
#Vk_minus_1  = Vk.shift(1)[1:]     #previous time point volatility here we also drop the final emntry
#Vk         = Vk[1:]               #drop final entry to align with line above
##aligh th arrays for stock rpice Sk
#Sk         =  ask[len(ask)-len(Vk)-1:] 
#Sk_minus_1 =  Sk.shift(1)[1:]
#Sk         = Sk[1:]
#
#
##************************** Heston  Model Estimation   ************************
#
#delta = 1  #1 minute this is time step size
#n     =  len(Vk) #This is the number of time points in minutes
#
#V_prod    =  Vk*Vk_minus_1   
#V_divide  =  Vk/Vk_minus_1 
#
##Theorem one in paper (parameters are labelled cosistent with out notation, not the paper's)
#P_numerator   = (1/n)*(np.sum(np.sqrt(V_prod))) - (1/(n*n))*np.sum(np.sqrt(V_divide))*np.sum(Vk_minus_1)  
#P_denominator = delta/2 - (delta/2)*(1/(n*n))*np.sum(1/Vk)*np.sum(Vk_minus_1)
#P             = P_numerator / P_denominator
#
#theta         = (2/delta)*( 1 + (P*delta/2)*(1/n)*np.sum(1/Vk_minus_1) - (1/n)*np.sum(np.sqrt(V_divide)) )
#eta           = np.sqrt( (4/delta)*(1/n)*np.sum( (np.sqrt(Vk) - np.sqrt(Vk_minus_1) - \
#                          (delta/(2*np.sqrt(Vk_minus_1)))*(P-theta*Vk_minus_1))**2 )  ) 
#alpha         = (P + eta*eta/4)/theta
#
#
##dW_1       = (np.log( Sk ) - np.log(Sk_minus_1) ) / np.sqrt(Vk_minus_1)
#dW_1       = ( Sk  - Sk_minus_1 ) / np.sqrt(Vk_minus_1)
#
#dW_2       = (Vk - Vk_minus_1 - theta*(alpha-Vk_minus_1)*delta)/(eta*np.sqrt(Vk_minus_1))
#rho        = (1/(n*delta))*np.sum(dW_1*dW_2)


def Heston_cali(price,var,delta):
    """
    This function takes input of
    price : array of stock prices
    var   : corresponding array of instanteous variance 
    delta : this is the time step between each point in time. i.e delta = 1mimute    
    
    len(price) must equal len(var)
    """
    
    Vk          = var.dropna()                      #dropping any nans
    Vk_90       = Vk.quantile(0.9)                 # 75th  quantile for upper cap.
    Vk = Vk.apply(lambda x: np.where(x < 20,20,x))  #removing near-zero values
    #Vk = Vk.apply(lambda x: np.where(x > Vk_90, Vk_90, x))   #removing numbers above
    Vk_minus_1  = Vk.shift(1)[1:]                  #previous time point volatility here we also drop the final emntry
    Vk          = Vk[1:]                           #drop final entry to align with line above
    n           = len(Vk)
    #aligh th arrays for stock price and Sk
    Sk         =  price[len(price)-len(Vk)-1:] 
    Sk_minus_1 =  Sk.shift(1)[1:]
    Sk         =  Sk[1:]
        
    V_prod    =  Vk*Vk_minus_1   
    V_divide  =  Vk/Vk_minus_1 
    
    #Theorem one in paper (parameters are labelled cosistent with out notation, not the paper's)
    P_numerator   = (1/n)*(np.sum(np.sqrt(V_prod))) - (1/(n*n))*np.sum(np.sqrt(V_divide))*np.sum(Vk_minus_1)  
    P_denominator = (delta/2)    - (delta/2)*(1/(n*n))*np.sum(1/Vk)*np.sum(Vk_minus_1)
    P             = P_numerator / P_denominator
    
    theta         = (2/delta)*( 1 + (P*delta/2)*(1/n)*np.sum(1/Vk_minus_1) - (1/n)*np.sum(np.sqrt(V_divide)) )
    eta           = np.sqrt( (4/delta)*(1/n)*np.sum( (np.sqrt(Vk) - np.sqrt(Vk_minus_1) - \
                              (delta/(2*np.sqrt(Vk_minus_1)))*(P-theta*Vk_minus_1))**2 )  ) 
    alpha         = (P + eta*eta/4)/theta
    
    
    #dW_1       = (np.log( Sk ) - np.log(Sk_minus_1) ) / np.sqrt(Vk_minus_1)
    dW_1       = (( Sk ) - (Sk_minus_1) ) / np.sqrt(Vk_minus_1)
    dW_2       = (Vk - Vk_minus_1 - theta*(alpha-Vk_minus_1)*delta)/(eta*np.sqrt(Vk_minus_1))
    rho        = (1/(n*delta))*np.sum(dW_1*dW_2)
    
    return theta, eta, alpha, rho 


#a,b,c,d = Heston_cali(ask,RVask30,1)


