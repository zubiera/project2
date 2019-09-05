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
df = df[:1000]


#create a column with the unix time.
df['unixtime'] = df.index
df = df.drop(('seq_num',0), axis=1)
#reset the unix time column to start as seconds from open of market.



#A hack to remove the double row colum headings, and then rename it all as single row column headings
df  = df.as_matrix()   #makes df into np array
df  = pd.DataFrame(df) #To make it back into a DF

df  = df.rename(columns={ 0: "b1",  1: "bv1",  2: "a1",  3: "av1",\
                                                   4: "b2",  5: "bv2",  6: "a2",  7: "av2",\
                                                   8: "b3",  9: "bv3", 10: "a3", 11: "av3",\
                                                  12: "b4", 13: "bv4", 14: "a4", 15: "av4",\
                                                  16: "b5", 17: "bv5", 18: "a5", 19: "av5",\
                                                  20:"unix_time"})
  
market_open_time = df['unix_time'][0]
df["time_from_open"] = df["unix_time"] - market_open_time





#Parameters used to estimate Lambda 
time = df["time_from_open"]
price = df["a1"]
plt.plot(time,price)


delta = 25  #this is the tick size


#initilise parameters for the algo
quote            = price[0] + delta    #The first quote posted
past_hit_time    = 0                   #Track the past order hit
waiting_time_sum = 0                   #Sum the waiting times
count            = 0                   #Initilise the times a quote is hit





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

Lambda_hat = count/waiting_time_sum



def Lambda_hat_ST_ask(df, ref_price, delta):
    
    """ 
    This function requires
    
    df: the data frame must have a time column and a price column
    ref_price:   This is the price we observe to see when it reaches our quoted price
    delta:       This is how far our quote is from the reference price. 
                 since the tick size for out data is 25, delta shoulde be multiples of 25. 
                 i.e if if our quote is 1 level away from ref_price, then delta is 1*25
    
    """
    
    #Extracting price and time vectors from df
    time  = df["time_from_open"]
    price = df[ref_price]
    
    
    #initilise parameters for the algo
    quote            = price[0] + delta    #The first quote posted
    past_hit_time    = 0                   #Track the past order hit
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit

    i=0 #initlise i
    #traverses throught he dataframe unti final time point
    while i < len(time):
        
        #only counts if the 'if' condition is met
        if price[i] >= quote:
            waiting_time_sum += time[i] - past_hit_time
            past_hit_time = time[i]
            count +=1
            quote = price[i] + delta
        #iterate i by one time point
        i+=1 
        
    #returns Lambda_hat = count/waiting_time_sum
    return count/waiting_time_sum
a = Lambda_hat_ST_ask(df, "a1", 25)
b = Lambda_hat_ST_ask(df, "a1", 50)
c = Lambda_hat_ST_ask(df, "a1", 75)
d = Lambda_hat_ST_ask(df, "a1", 100)
e = Lambda_hat_ST_ask(df, "a1", 125)
f = Lambda_hat_ST_ask(df, "a1", 150)
g = Lambda_hat_ST_ask(df, "a1", 175)
h = Lambda_hat_ST_ask(df, "a1", 200)


#plotting Lambda for different deltas
L_hat_ask = [a,b,c,d,e,f,g,h]
vary_delt = np.array([25,50,75,100,125,150,175,200])
plt.plot(vary_delt,L_hat_ask,"-")
plt.title("$\hat\Lambda(\delta^a)$ for ask sid of the LOB")
plt.xlabel("$\delta$")
plt.ylabel("$\hat\Lambda(\delta)$")






#•••••••••••••••••••• The same as the above but for bid •••••••••••••••••••••••••



def Lambda_hat_ST_bid(df, ref_price, delta):
    
    """ 
    This function requires
    
    df: the data frame must have a time column and a price column
    ref_price:   This is the price we observe to see when it reaches our quoted price
    delta:       This is how far our quote is from the reference price. 
                 since the tick size for out data is 25, delta shoulde be multiples of 25. 
                 i.e if if our quote is 1 level away from ref_price, then delta is 1*25
    
    """
    
    #Extracting price and time vectors from df
    time  = df["time_from_open"]
    price = df[ref_price]
    
    
    #initilise parameters for the algo
    quote            = price[0] - delta    #The first quote posted
    past_hit_time    = 0                   #Track the past order hit
    waiting_time_sum = 0                   #Sum the waiting times
    count            = 0                   #Initilise the times a quote is hit

    i=0 #initlise i
    #traverses throught he dataframe unti final time point
    while i < len(time):
        
        #only counts if the 'if' condition is met
        if price[i] <= quote:
            waiting_time_sum += time[i] - past_hit_time
            past_hit_time = time[i]
            count +=1
            quote = price[i] - delta
        #iterate i by one time point
        i+=1 
        
    #returns Lambda_hat = count/waiting_time_sum
    return count/waiting_time_sum

b1 = Lambda_hat_ST_bid(df, "b1", 25)
b2 = Lambda_hat_ST_bid(df, "b1", 50)
b3 = Lambda_hat_ST_bid(df, "b1", 75)
b4 = Lambda_hat_ST_bid(df, "b1", 100)
b5 = Lambda_hat_ST_bid(df, "b1", 125)
b6 = Lambda_hat_ST_bid(df, "b1", 150)
b7 = Lambda_hat_ST_bid(df, "b1", 175)
b8 = Lambda_hat_ST_bid(df, "b1", 200)


# Plotting the different  Lambdas
L_hat_bid = [b1,b2,b3,b4,b5,b6,b7,b8]
vary_delt = np.array([25,50,75,100,125,150,175,200])
plt.plot(vary_delt,L_hat_bid,"-")
plt.title("$\hat\Lambda(\delta^b)$ for bid sid of the LOB")
plt.xlabel("$\delta$")
plt.ylabel("$\hat\Lambda(\delta)$")



plt.plot([b1,b2,b3,b4,b5,b6])
plt.plot([a,b,c,d,e,f])
plt.title("$\hat\Lambda(\delta)$ for bid and ask sid of the LOB")
plt.xlabel("$\delta$")
plt.ylabel("$\hat\Lambda(\delta)$")


#********************** Comments to self about the code above ***************

#   I think the code is assuming a pessimistic case. I.e for level one quote
#   we assume we are at the back of the queue since we wait for price to jump to level
#   two before we count it as being hit.

