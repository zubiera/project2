import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


#setting the file path.
from sys import platform
path   = "/Users/amitarfan1/Documents/Phd/2yrReport/"
if platform == "linux" or platform == "linux2":
    path   = "/home/zubier/Desktop/Phd/2yrReport/"
elif platform == "darwin":
    path   = "/Users/amitarfan1/Documents/Phd/2yrReport/"
    

data = pd.HDFStore(path + 'data/170707#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
df = data['data']
df = df[:100]


#create a column with the unix time.
df['unixtime'] = df.index
df = df.drop(('seq_num',0), axis=1)
#reset the unix time column to start as seconds from open of market.



#A hack to remove the double row colum headings, and then rename it all as single row column headings
df  = df.as_matrix()   #makes df into np array
df  = pd.DataFrame(df) #To make it back into a DF

df  =df.rename(columns={ 0: "b1",  1: "bv1",  2: "a1",  3: "av1",\
                                                   4: "b2",  5: "bv2",  6: "a2",  7: "av2",\
                                                   8: "b3",  9: "bv3", 10: "a3", 11: "av3",\
                                                  12: "b4", 13: "bv4", 14: "a4", 15: "av4",\
                                                  16: "b5", 17: "bv5", 18: "a5", 19: "av5",\
                                                  20:"unix_time"})
        
    
    
    
market_open_time = df['unix_time'][0]
df["time_from_open"] = df["unix_time"] - market_open_time







#this is my working progress for restarting waiting time

#The below array is a mock time series to test any corner cases for my plotting.
time = np.array([0,1,2,3,4,5,6,7,11,20])
price = np.array([25,26,27,29,27,26,27,26,27,31])

plt.plot(time,price)
delta = 2

#initial quote at market open
quote = price[0]+ delta

#waiting time initilisations
past_hit_time = 0
waiting_time_sum =0
#count intililise
count=0
i=0
while i < len(time):
    if price[i] >= quote:
        waiting_time_sum += time[i] - past_hit_time
        past_hit_time = time[i]
        count +=1
        quote = price[i] + delta
    i+=1
        
waiting_time_sum
count




ts = [1,2,3,4,5,7,10,12,15,17,18,20]
S  = [1,1,1,1,0,0,0,1,0,0,0,0] 
N=4
dt = ts[-1]/N
counter = 0

for i in range(0,N):
    #S_a1 referes to the fixed 1 tick distance from ask_1 at time i*time_step

    for j in range(np.searchsorted(ts,dt*i), np.searchsorted(ts,dt*(i+1))):
        if S[j]==1:
            counter += 1                  #easier to debug and check intervals if needed
            break
    



time = np.array([0,1,2,3,4,5,6,7,11,20,21.5,24])
price = np.array([27,27,27,27,27,27,27,27,27,27,27,27])
N=4
delta=2
#size of time interval. ***** Might need to remove this when doing diff days
#this is because dt will vary on each day.
dt = time[len(time)-1]/N

#initilise parameters for the algorithm
#quote            = price[0] - delta    #The first quote posted
waiting_time_sum = 0                   #Sum the waiting times
count            = 0                   #Initilise the times a quote is hit

#traverse through the day 
for i in range(0,N):
    quote = price[np.searchsorted(time,dt*i)] + delta
    
    order_hit = False
    #traverse through each time interval and break when threshold is hit
    for j in range(np.searchsorted(time,dt*i)+1, np.searchsorted(time,dt*(i+1))+1):
        if price[j] >= quote:
            count += 1
            waiting_time_sum += time[j] - time[np.searchsorted(time,dt*i)]
            order_hit = True
            break
    if order_hit == False:
        waiting_time_sum += dt
        
        
    

count
waiting_time_sum


def Lambda_ST(time,price,delta,N):
    
    
    
    
    dt = time[len(time)-1]/N               # Size of Time interval
    
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    #Traverse through the day 
    for i in range(0,N):
        
        # np.searchsorted finds the closes time index to the time point dt*i
        quote = price[np.searchsorted(time,dt*i)] + delta
        
        order_hit = False                  # This is to dertimin the waiting time.
        
        #traverse through each time interval and break when threshold is hit
        for j in range(np.searchsorted(time,dt*i)+1, np.searchsorted(time,dt*(i+1))+1):
            if price[j] >= quote:
                count += 1
                waiting_time_sum += time[j] - time[np.searchsorted(time,dt*i)]
                order_hit = True
                break
        # if there is no order hit we choose waiting time to be the full dt.
        if order_hit == False:
            waiting_time_sum += dt
    
    return waiting_time_sum
    

time = np.array([0,1,2,3,4,5,6,7,11,20,21.5,24])
price = np.array([27,27,27,27,27,29,29,27,27,27,27,27])
N=4
delta=2
Lambda_ST(time,price,delta,N)











