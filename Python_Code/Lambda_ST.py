import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


from data_load import create_df


df = create_df('data/170707#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')


    
def Lambda_ST_ask(time,price,delta,N):
    
    """
    Parameters.
    df:        A data frame of the HFT trading data. Must contain a column for ref price, and time
    ref_price: Enter this as string, i.e for a1 enter "a1"
    delta:     This is how far from the ref price we set a quote.
    N:         Number of time Intervals 
    
    """
     
    dt = time[len(time)-1]/N               # Size of Time interval
    
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    #Traverse through the day 
    for i in range(0,N):
        
        # np.searchsorted finds the closes time index to the time point dt*i
        quote = price[np.searchsorted(time,dt*i)] + delta
        
        order_hit = False                  # This is to dertimin the waiting time.
        
        #traverse through each time interval and break when threshold is hit
        for j in range(np.searchsorted(time,dt*i), np.searchsorted(time,dt*(i+1))):
            if price[j] >= quote:
                count += 1
                waiting_time_sum += time[j] - time[np.searchsorted(time,dt*i)]
                order_hit = True
                break
        # if there is no order hit we choose waiting time to be the full dt.
        if order_hit == False:
            waiting_time_sum += dt
    
    return count/waiting_time_sum

def Lambda_ST_bid(time,price,delta,N):
    
    """
    Parameters.
    df:        A data frame of the HFT trading data. Must contain a column for ref price, and time
    ref_price: Enter this as string, i.e for a1 enter "a1"
    delta:     This is how far from the ref price we set a quote.
    N:         Number of time Intervals 
    
    """
     
    dt = time[len(time)-1]/N               # Size of Time interval
    
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    #Traverse through the day 
    for i in range(0,N):
        
        # np.searchsorted finds the closes time index to the time point dt*i
        quote = price[np.searchsorted(time,dt*i)] - delta
        
        order_hit = False                  # This is to dertimin the waiting time.
        
        #traverse through each time interval and break when threshold is hit
        for j in range(np.searchsorted(time,dt*i), np.searchsorted(time,dt*(i+1))):
            if price[j] <= quote:
                count += 1
                waiting_time_sum += time[j] - time[np.searchsorted(time,dt*i)]
                order_hit = True
                break
        # if there is no order hit we choose waiting time to be the full dt.
        if order_hit == False:
            waiting_time_sum += dt
    
    return count/waiting_time_sum



#time = np.array([0,1,2,3,4,5,6,7,11,20,21.5,24])
#price = np.array([27,27,27,27,27,29,29,27,27,27,27,27])



time  = np.array(df["time_from_open"])
ask_price = np.array(df["a1"])
bid_price = np.array(df["b1"])

L_vec_ask = [0]*10
L_vec_bid = [0]*10
for i in range (0,10):
    L_vec_ask[i] = Lambda_ST_ask(time,ask_price,25*(i+1),1380)
    L_vec_bid[i] = Lambda_ST_bid(time,bid_price,25*(i+1),1380)

#make into Numpy array and multiply by 1e9 to make nan0sec into sec. then multiply by 60 to make into minute
L_vec_ask = np.array(L_vec_ask)*1e9*60
L_vec_bid = np.array(L_vec_bid)*1e9*60


delta_vec = np.arange(25,275,25)
plt.plot(delta_vec[0:],L_vec_ask[0:],  label="$\hat\Lambda(\delta^a)$",color = 'blue')
plt.plot(delta_vec[0:],L_vec_bid[0:],  label="$\hat\Lambda(\delta^b)$", color = 'red')
plt.title("$\hat\Lambda(\delta)$ for $\Delta t = 60$sec")
plt.xlabel("$\delta$")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.legend()

#Log plot  
delta_vec = np.arange(25,275,25)
plt.plot(delta_vec,np.log(L_vec_ask),label="$\hat\Lambda(\delta^a)$")
plt.plot(delta_vec,np.log(L_vec_bid),label="$\hat\Lambda(\delta^b)$")
plt.title("$\hat\Lambda(\delta)$ for bid and ask sid of the LOB")
plt.xlabel("$\delta$")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.legend()
    

#linear regression
x = delta_vec[0:]
y_ask = np.log(L_vec_ask[0:])
y_bid = np.log(L_vec_bid)
#We can rewrite the line equation as y = Ap, where A = [[x 1]] and p = [[m], [c]]. Now use lstsq to solve for p:
A_ask = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A_ask, y_ask)[0]
#plt.plot(x, y_ask, 'o', label='Original data', markersize=10)
#plt.plot(x, m*x + c, 'r', label='Fitted line')
#plt.legend()



#Plotting the  regression fit 
A =  np.exp(c)
k = -m
regression_ask = A*np.exp(-k*x)
plt.plot(x,regression_ask,'--',color='red',label = "Regression")
plt.plot(x,L_vec_ask[0:],color= 'blue',label = "$\hat\Lambda$")
plt.legend()



k = (np.log(L_vec_ask[1]/L_vec_ask[3])) / (50) 
    



#*************************Plotting different time interval sizes****************
Ld1 = []
Ld2 = []
Ld3 = []
Ld4 = []


#create dt to increment from 1 second to 3600 seconds (1hour)
dt = np.arange(1,3601,1)*1e9
N = (np.round(time[-1]/dt))



#The above has duplicates, so the belwo removes thema and prserves order
#Note: pandas may be a lot faster at doing this!
_, idx = np.unique(N, return_index=True)
N = N[np.sort(idx)]  
dt = time[-1]/N    #re creating dt to amtch length of N.
   
for i in N:
    #Ld1.append(Lambda_ST_ask(time,ask_price,25,int(i)))
    Ld2.append(Lambda_ST_ask(time,ask_price,50,int(i)))
    Ld3.append(Lambda_ST_ask(time,ask_price,75,int(i)))
    Ld4.append(Lambda_ST_ask(time,ask_price,100,int(i)))


#N = np.arange(100000,2100000,100000)
#ld1_microseconds = Ld1

plt.plot(dt/1e9,Ld1,color='blue')
plt.title("$\hat\Lambda(\delta)$ for ask sid of the LOB for $\delta=25$")
plt.xlabel("$\Delta t$ (seconds)")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.grid()
plt.rc('grid', linestyle="--", color='grey')

plt.plot(dt/1e9,Ld2,color='blue')
plt.title("$\hat\Lambda(\delta)$ for ask sid of the LOB for $\delta=50$")
plt.xlabel("$\Delta t$ (seconds)")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.grid()
plt.rc('grid', linestyle="--", color='grey')

plt.plot(dt/1e9,Ld3,color='blue')
plt.title("$\hat\Lambda(\delta)$ for ask sid of the LOB for $\delta=75$")
plt.xlabel("$\Delta t$ (seconds)")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.grid()
plt.rc('grid', linestyle="--", color='grey')

plt.plot(dt/1e9,Ld4,color='blue')
plt.title("$\hat\Lambda(\delta)$ for ask sid of the LOB for $\delta=100$")
plt.xlabel("$\Delta t$ (seconds)")
plt.ylabel("$\hat\Lambda(\delta)$")
plt.grid()
plt.rc('grid', linestyle="--", color='grey')

Lambda_ST_ask(time,ask_price,25,2000000)

