import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from statsmodels.tsa.arima_model import ARIMA

from data_load import create_df

#************replicate data frame to have it in 1 minute slots ***********

days = ['03','04','05','06','07','10','11','12','13','14','17',\
                  '18', '19','20','21','24','25','26','27','28','31']

days = ['05','06','07','10','11','12','13','14']
time_mins = {}
ask     = pd.DataFrame()           #minutely ask
bid     = {}           #minutely bid
diffa_1 = {}           #differecnes ask price 1 min
diffa_2 = {}           #differecnes ask price 2 min
diffa_5 = {}           #differecnes ask price 5 min
double_diffa = {}
RVask10    = pd.DataFrame()           #rolling variance 10 min
RVask30    = pd.DataFrame()           #rolling variance 30 min
RVask60    = pd.DataFrame()           #rolling variance 30 min
ask_volume = pd.DataFrame()
 
for j in range (len(days)):
    
    df = create_df('data/1707'+ days[j] +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
                   
    time      = np.array(df["time_from_open"])
    ask_price = np.array(df["a1"])
    bid_price = np.array(df["b1"]) 
    bid_vol   = df["bv1"] + df["bv2"] + df["bv3"] + df["bv4"] + df["bv5"]   
    ask_vol   = df["av1"] + df["av2"] + df["av3"] + df["av4"] + df["av5"]
    dt=6e10   
    S_a  = np.array(ask_price[0])  
    S_b  = np.array(bid_price[0])
    av   = np.array(ask_vol[0])
    bv   = np.array(bid_vol[0])
    minutes = np.array([0])   
    i=0
    while dt*i <  time[-1]:
        
        S_a =   np.append (S_a, [ask_price[np.searchsorted(time,dt*i)]]) 
        S_b =   np.append (S_b, [bid_price[np.searchsorted(time,dt*i)]])
        bv  =   np.append (bv,  [bid_vol[np.searchsorted(time,dt*i)]])
        av  =   np.append (av,  [ask_vol[np.searchsorted(time,dt*i)]])
        minutes =   np.append (minutes, i)
        i+=1
    
    S_a = np.delete(S_a, 0)  #Remove first dummy eleemnt
    S_b = np.delete(S_b, 0)  #Remove first dummy element
    av  = np.delete(av,0)    #Remove first dummy element
    bv  = np.delete(bv,0)    #Remove first dummy element
    minutes = np.delete(minutes, 0)  #Remove first dummy element
    
    
    S_a = pd.Series(S_a)    #make into a panda time series
    #S_a = S_a.pct_change()  #caclulate minutly retursn
    S_b = pd.Series(S_b)    #make into a panda time series
    #S_b = S_b.pct_change()  #caclulate minutly retursn
    bv  = pd.Series(bv)
    av  = pd.Series(av)
    minutes = pd.Series(minutes)
    
    #all the ask-side series
    time_mins["{0}".format(days[j])] = minutes
    ask[days[j]]                     = S_a
    
    diffa_1["{0}".format(days[j])]    = S_a.diff(1)
    diffa_2["{0}".format(days[j])]    = S_a.diff(2)
    diffa_5["{0}".format(days[j])]    = S_a.diff(5)
    double_diffa["{0}".format(days[j])] = S_a.diff(1).diff(1)
    RVask10[days[j]]                    = S_a.rolling(10).var()
    RVask30[days[j]]                    = S_a.rolling(30).var()
    RVask60[days[j]]                    = S_a.rolling(60).var()
    ask_volume[days[j]]                 = av
    
    fig = plt.figure()
    
    plt.subplot(3, 2, 1)
    plt.plot(time_mins[days[j]], ask_volume[days[j]], linewidth=0.3, color ='red' )
    plt.title('Ask Volume')
    plt.ylabel('Volume')
    plt.xlabel('time (mins)')
    
    plt.grid()
    
    plt.subplot(3, 2, 2)
    plt.plot(time_mins[days[j]], diffa_1[days[j]], linewidth=0.3, color ='red')
    plt.title('differenced t=1')
    plt.ylabel('$S_t-S_{t-1}$')
    plt.xlabel('time (mins)')
    plt.grid()
    
    plt.subplot(3, 2, 3)
    plt.plot(time_mins[days[j]], RVask10[days[j]], linewidth=0.3, color ='red')
    plt.title('Rolling 10min Variance')
    plt.ylabel('Variance')
    plt.xlabel('time (mins)')
    plt.grid()
    
    plt.subplot(3, 2, 4)
    plt.plot(time_mins[days[j]],  diffa_2[days[j]], linewidth=0.3, color ='red')
    plt.title('differenced t=2')
    plt.ylabel('$S_t-S_{t-2}$')
    plt.xlabel('time (mins)')
    plt.grid()
    
    plt.subplot(3, 2, 5)
    plt.plot(time_mins[days[j]],  RVask30[days[j]], linewidth=0.3, color ='red')
    plt.title('Rolling 30min Variance')
    plt.ylabel('Variance')
    plt.xlabel('time (mins)')
    plt.grid()
    
    plt.subplot(3, 2, 6)
    plt.plot(time_mins[days[j]],  diffa_5[days[j]], linewidth=0.3, color ='red')
    plt.title('differenced t=5')
    plt.ylabel('$S_t-S_{t-5}$')
    plt.xlabel('time (mins)')
    plt.grid()
    
    plt.tight_layout()
    plt.show()



#***********Plotting the volatility average
RV10_average = RVask10.mean(axis=1) 
RV30_average = RVask30.mean(axis=1) 
RV60_average = RVask60.mean(axis=1) 
RV10_average.plot()
RV30_average.plot(color = "purple")
RV60_average.plot()
#variance for first segment of day, for each day (divide by N degreefreedom=0)
sigma_test = np.sqrt( ask.var(axis=0,ddof=0)/1380)

sigma1 = np.sqrt( ask[:550].var(axis=0,ddof=0)/550)
sigma2 = np.sqrt(ask[550:850].var(axis=0,ddof=0)/(850-550))
sigma3 = np.sqrt(ask[850:1100].var(axis=0,ddof=0)/(1100-850))
sigma4 = np.sqrt(ask[1100:1380].var(axis=0,ddof=0)/(1380-1100))


sigma1_mean = np.mean(sigma1)        # get average variance
sigma2_mean = np.mean(sigma2)        # get average variance

plt.plot(RVask30["12"])
plt.plot(ask["12"])
plt.plot(RV30_average)
plt.plot(RVask30["07"])
plt.plot(RVask1380["05"])



#***************  Do some time series stuff ••••••••••••••

from pandas.tools.plotting import autocorrelation_plot
plt.plot(S_a)
autocorrelation_plot(S_a)



