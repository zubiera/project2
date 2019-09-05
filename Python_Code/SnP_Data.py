import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt



from sys import platform
path   = "/Users/amitarfan1/Documents/Phd/2yrReport/"
if platform == "linux" or platform == "linux2":
    path   = "/home/zubier/Desktop/Phd/2yrReport/"
elif platform == "darwin":
    path   = "/Users/amitarfan1/Documents/Phd/2yrReport/"
    
    

#*********************************Grabbng the Data*************************************


data = pd.HDFStore(path + 'data/170707#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
df = data['data']
#df = df[:1000]
pd.to_datetime(df.index[0], utc=True)
pd.to_datetime(df.index[-1], utc=True)
df = df.drop(('seq_num',0), axis=1) 
#Add in a time conversion 



#•••••••••••••••••••••••••••••••••••••••  Clearning up the column names ••••••••••••••••••••••••••
#A hack to remove the double row colum headings, and then rename it all as single row column headings
df  = df.as_matrix()   #makes df into np array
df  = pd.DataFrame(df) #To make it back into a DF


#This is to rename the colums becasue the two lines above
df  =df.rename(columns={ 0: "b1",  1: "bv1",  2: "a1",  3: "av1",\
                                               4: "b2",  5: "bv2",  6: "a2",  7: "av2",\
                                               8: "b3",  9: "bv3", 10: "a3", 11: "av3",\
                                              12: "b4", 13: "bv4", 14: "a4", 15: "av4",\
                                              16: "b5", 17: "bv5", 18: "a5", 19: "av5"})
    
    
    
#create a mide price.   
df['midprice'] = ( df['b1'] + df['a1'] ) / 2   


#Define a column that shows total volume on bid and ask side
df['total_bid_vol'] = (df['bv1'] + df['bv2'] + df['bv3'] + df['bv4'] + df['bv5'])
df['total_ask_vol'] = (df['av1'] + df['av2'] + df['av3'] + df['av4'] + df['av5'])


#•••••••••••••••••Some calculations to help estimate order arrival intensity••••••••••••••
#Labels 1 if first level of ask and bid volume is depleted in a SINGLE TRADE
df['a1count'] = 0
df['a1count'][ ((df['a1']-df['a1'].shift(1))/ 25 == 1) ] = 1
df['b1count'] = 0
df['b1count'][ ((df['b1']-df['b1'].shift(1))/ 25 == -1) ] = 1

#Labels 1 if first two levels of ask and bid volume is depleted
df['a2count'] = 0
df['a2count'][ ((df['a1']-df['a1'].shift(1))/ 25 == 2) ] = 1
df['b2count'] = 0
df['b2count'][ ((df['b1']-df['b1'].shift(1))/ 25 == -2) ] = 1

#Labels 1 if first two levels of ask and bid volume is depleted
df['a3count'] = 0
df['a3count'][ ((df['a3']-df['a3'].shift(1))/ 25 == 3) ] = 1
df['b3count'] = 0
df['b3count'][ ((df['b3']-df['b3'].shift(1))/ 25 == -3) ] = 1

#Labels 1 if first two levels of ask and bid volume is depleted
df['a4count'] = 0
df['a4count'][ ((df['a4']-df['a4'].shift(1))/ 25 == 4) ] = 1
df['b4count'] = 0
df['b4count'][ ((df['b4']-df['b4'].shift(1))/ 25 == -4) ] = 1




#Intensity Estimates
lambda_b1 = sum(df['b1count']) / len(df['b1count'])
lambda_b2 = sum(df['b2count']) / len(df['b2count'])
lambda_b3 = sum(df['b3count']) / len(df['b3count'])
lambda_b4 = sum(df['b4count']) / len(df['b4count'])

lambda_b    = [0]*4
lambda_b[0] = lambda_b1
lambda_b[1] = lambda_b2
lambda_b[2] = lambda_b3
lambda_b[3] = lambda_b4


df.index[df['b2count']==1].tolist()


