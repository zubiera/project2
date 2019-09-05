import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

from data_load import create_df
               
#************************ Optimistic Estimation Functions********************              

def Ask_Opti(time_array,price_array, vol_array, vol_array2, delta, dt):
    """
    Parameters.
    time_array:   1d array of time in nanoseconds from market open time.
    price_array:  1d array of best ask price (level 1 ask)
    vol_array:    1d array of volume corresposning to level 1 ask
    vol_array2:   1d array of volume corresposning to level 2 ask
    delta:        This is how far from the ref price we set a quote.
    dt:           Time interval size in nanoseconds. i.e 60 seconds is 60e9 nanosecs
    
    """


    #This is used for debugging only.
    tracker  = np.array([['idx','time','quote',"delta","type","count","time"]])
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    
    
    i=0
    while np.searchsorted(time_array,dt*(i+1))<len(time_array):
        
        quote = price_array[np.searchsorted(time_array,dt*i)] + delta
        order_hit = False                  # This is to dertimin the waiting time
        
        #traverse through elements in
        for j in range(np.searchsorted(time_array,dt*i), np.searchsorted(time_array,dt*(i+1))):
            
            #Corner case for delta =0 (setting our quote at current best quote, 
            #this means we wait for some volume depletion befcore counting the execution)
            if (delta == 0):
                if ( (quote==price_array[j]==price_array[j+1]) \
                and vol_array[j] > vol_array[j+1] ) or price_array[j+1]>quote:
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote, delta, "0",count, waiting_time_sum]])))               
                    order_hit = True
                    break 
             #THe case if the ref price becomes equal to the quote, again wait for volume 
             #depletion before counting execution
            elif quote == price_array[j+1]:
                #The case an MO walks 1 level trhough the book, eating some 
                #volume where our quote is
                if (price_array[j+1]==quote) and (price_array[j] < quote) and (vol_array[j+1]<vol_array2[j]):
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"MO 1step",count,waiting_time_sum]])))
                    order_hit = True
                    break 
                #The case where our quote becomes equal to best ask, but no volume
                #has been been consumed just yet.
                elif price_array[j] == price_array[j+1] and vol_array[j] > vol_array[j+1]:
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"small MO",count,waiting_time_sum]])))
                    order_hit = True
                    break 
            #the case where a large MO results in the ref price becoming larger
            #than the quote (i.e it consumes our quote while walking the book)
            elif quote <  price_array[j+1]:
                count += 1
                waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"big MO",count,waiting_time_sum]])))
                order_hit = True
                break
        
        if order_hit == False:
            waiting_time_sum += dt
        i+=1
    
    return count/waiting_time_sum *  60e9 #this makes the intensity per minute. 

def Bid_Opti(time_array,price_array, vol_array, vol_array2, delta, dt):
    """
    Parameters.
    time_array:   1d array of time in nanoseconds from market open time.
    price_array:  1d array of best ask price (level 1 ask)
    vol_array:    1d array of volume corresposning to level 1 ask
    vol_array2:   1d array of volume corresposning to level 2 ask
    delta:        This is how far from the ref price we set a quote.
    dt:           Time interval size in nanoseconds. i.e 60 seconds is 60e9 nanosecs
    
    """


    #This is used for debugging only.
    tracker  = np.array([['idx','time','quote',"delta","type","count","time"]])
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    
    
    i=0
    while np.searchsorted(time_array,dt*(i+1))<len(time_array):
        
        quote = price_array[np.searchsorted(time_array,dt*i)] - delta
        order_hit = False                  # This is to dertimin the waiting time
        
        #traverse through elements in
        for j in range(np.searchsorted(time_array,dt*i), np.searchsorted(time_array,dt*(i+1))):
            
            #Corner case for delta =0 (setting our quote at current best quote, 
            #this means we wait for some volume depletion befcore counting the execution)
            if (delta == 0):
                if ( (quote==price_array[j]==price_array[j+1]) \
                and vol_array[j] > vol_array[j+1] ) or price_array[j+1]<quote:
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote, delta, "0",count, waiting_time_sum]])))               
                    order_hit = True
                    break 
             #THe case if the ref price becomes equal to the quote, again wait for volume 
             #depletion before counting execution
            elif quote == price_array[j+1]:
                #The case an MO walks a level through the book, eating some 
                #volume where our quote is
                if (price_array[j+1]==quote) and (price_array[j] > quote) and (vol_array[j+1]<vol_array2[j]):
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"MO 1step",count,waiting_time_sum]])))
                    order_hit = True
                    break 
                #The case where our quote becomes equal to best ask, but no volume
                #has been been consumed just yet.
                elif price_array[j] == price_array[j+1] and vol_array[j] > vol_array[j+1]:
                    count += 1
                    waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                    tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"small MO",count,waiting_time_sum]])))
                    order_hit = True
                    break 
            #the case where a large MO results in the ref price becoming larger
            #than the quote (i.e it consumes our quote while walking the book)
            elif quote >  price_array[j+1]:
                count += 1
                waiting_time_sum += time_array[j+1] - time_array[np.searchsorted(time_array,dt*i)]
                tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j+1]),time_array[j+1] ,quote,delta,"big MO",count,waiting_time_sum]])))
                order_hit = True
                break
        
        if order_hit == False:
            waiting_time_sum += dt
        i+=1
    
    return count/waiting_time_sum *  60e9 #this makes the intensity per minute. 



#**************************** Pessimistic estimation Functions *******************


def Ask_Pessi(time_array,price_array, vol_array, vol_array2, delta, dt):
    """
    Parameters.
    time_array:   1d array of time in nanoseconds from market open time.
    price_array:  1d array of best ask price (level 1 ask)
    vol_array:    1d array of volume corresposning to level 1 ask
    vol_array2:   1d array of volume corresposning to level 2 ask
    delta:        This is how far from the ref price we set a quote.
    dt:           Time interval size in nanoseconds. i.e 60 seconds is 60e9 nanosecs
    
    """


    #This is used for debugging only.
    tracker  = np.array([['idx','time','quote',"delta","type","count","time"]])
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    
    
    i=0
    while np.searchsorted(time_array,dt*(i+1))<len(time_array):
        
        quote = price_array[np.searchsorted(time_array,dt*i)] + delta
        order_hit = False                  # This is to dertimin the waiting time
        
        #traverse through elements in
        for j in range(np.searchsorted(time_array,dt*i), np.searchsorted(time_array,dt*(i+1))):
            if price_array[j] > quote:
                count += 1
                waiting_time_sum += time_array[j] - time_array[np.searchsorted(time_array,dt*i)]
                tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j]),time_array[j] ,quote, delta, "0",count, waiting_time_sum]])))               
                order_hit = True
                break 
             #THe case if the ref price becomes equal to the quote, again wait for volume 

        if order_hit == False:
            waiting_time_sum += dt
        i+=1
    
    return count/waiting_time_sum *  60e9 #this makes the intensity per minute. 

def Bid_Pessi(time_array,price_array, vol_array, vol_array2, delta, dt):
    """
    Parameters.
    time_array:   1d array of time in nanoseconds from market open time.
    price_array:  1d array of best ask price (level 1 ask)
    vol_array:    1d array of volume corresposning to level 1 ask
    vol_array2:   1d array of volume corresposning to level 2 ask
    delta:        This is how far from the ref price we set a quote.
    dt:           Time interval size in nanoseconds. i.e 60 seconds is 60e9 nanosecs
    
    """


    #This is used for debugging only.
    tracker  = np.array([['idx','time','quote',"delta","type","count","time"]])
    
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit
    
    
    
    i=0
    while np.searchsorted(time_array,dt*(i+1))<len(time_array):
        
        quote = price_array[np.searchsorted(time_array,dt*i)] - delta
        order_hit = False                  # This is to dertimin the waiting time
        
        #traverse through elements in
        for j in range(np.searchsorted(time_array,dt*i), np.searchsorted(time_array,dt*(i+1))):
            if price_array[j] < quote:
                count += 1
                waiting_time_sum += time_array[j] - time_array[np.searchsorted(time_array,dt*i)]
                tracker = np.concatenate((tracker, np.array([[np.searchsorted(time_array,time_array[j]),time_array[j] ,quote, delta, "0",count, waiting_time_sum]])))               
                order_hit = True
                break 
             #THe case if the ref price becomes equal to the quote, again wait for volume 

        if order_hit == False:
            waiting_time_sum += dt
        i+=1
    
    return count/waiting_time_sum * 60e9 #this makes the intensity per minute.  




#******************************Run the Functions here*********************

#
df = create_df('data/170707#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
time      = np.array(df["time_from_open"])  
ask = np.array(df["a1"])
avol   = np.array(df["av1"])
avol2  = np.array(df["av2"])
bid = np.array(df["b1"])
bvol   = np.array(df["bv1"])
bvol2  = np.array(df["bv2"])
dt = 60e9


A=Ask_Opti(time,ask,avol,avol2,50,dt)
B=Bid_Opti(time,bid,bvol,bvol2,50,dt)

Optimistic_ask_intensity  = [0]*5
Pessimistic_ask_intensity = [0]*5

for i in range(0,5):
    Optimistic_ask_intensity[i]  = Ask_Opti(time,ask,avol,avol2,i*25,dt)
    Pessimistic_ask_intensity[i] = Ask_Pessi(time,ask,avol,avol2,i*25,dt)
    
  