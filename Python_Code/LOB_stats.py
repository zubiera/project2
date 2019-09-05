import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


from data_load import create_df

#create df
df = create_df('data/170707#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')

time  = np.array(df["time_from_open"])
ask_price = np.array(df["a1"])
bid_price = np.array(df["b1"])            
total_ask_vol = np.array(df["av1"]) +np.array(df["av2"])+np.array(df["av3"])+np.array(df["av4"]) +np.array(df["av5"])
total_bid_vol = np.array(df["bv1"]) +np.array(df["bv2"])+np.array(df["bv3"])+np.array(df["bv4"]) +np.array(df["bv5"])

plt.plot(time,total_ask_vol,linewidth=0.2, color='red', label='Ask Volume')
plt.plot(time,-total_bid_vol,linewidth=0.2, color = 'blue', label='Bid Volume')
plt.title('Total bid and ask volume Throughout the day.')
plt.ylabel('Volume ($50 Contract Units)')
plt.xlabel("Time from Market Open (nanoseconds)")
plt.legend()
plt.tight_layout()


#****************************** Point Estimation Method ******************************



def mean_variance(day,dt):
    df = create_df('data/1707'+ day +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
                   
    time  = np.array(df["time_from_open"])
    ask_price = np.array(df["a1"])
    bid_price = np.array(df["b1"])  

    
    Sa_sq_sum = 0   #reference price squared
    Sa_sum    = 0   #sum of refernce price
    Sb_sq_sum = 0
    Sb_sum    = 0
    
    i=0
    while dt*i <  time[-1]:
        S_a =   ask_price[np.searchsorted(time,dt*i)] 
        S_b =   bid_price[np.searchsorted(time,dt*i)]
        
        #sum of sq and sum for asks
        Sa_sq_sum += S_a**2
        Sa_sum    += S_a
        
        #sum of sq and sum for bids
        Sb_sq_sum += S_b**2
        Sb_sum    += S_b
        
        
        i+=1
    
    sigma_a_sq = Sa_sq_sum/i - (Sa_sum/i)**2
    sigma_b_sq = Sb_sq_sum/i - (Sb_sum/i)**2
        
    mu_a = (ask_price[-1]-ask_price[0]) / i  
    mu_b = (bid_price[-1]-bid_price[0]) / i  
    
    return i, mu_a, sigma_a_sq, mu_b, sigma_b_sq


days = ['03','04','05','06','07','10','11','12','13','14','17',\
                  '18', '19','20','21','24','25','26','27','28','31']
stats = pd.DataFrame(columns=['day','No_intervals' ,'mu_a', 'sigma_a_sq', 'mu_b', 'sigma_b_sq'])

for i in range(len(days)):
    intervals, mu_a, sigma_a_sq, mu_b, sigma_b_sq = mean_variance(days[i],60e9)
    df = pd.DataFrame([[days[i], intervals, mu_a, sigma_a_sq, mu_b, sigma_b_sq]], columns =['day','No_intervals' ,'mu_a', 'sigma_a_sq', 'mu_b', 'sigma_b_sq'] )
    stats = stats.append(df)
    


df = create_df('data/1707'+ '07' +'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
                   
time  = np.array(df["time_from_open"])
ask_price = np.array(df["a1"])
bid_price = np.array(df["b1"])  
plt.plot(time,ask_price)