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
    
#••••••••••••••••••••••••••••••Below is  just data preprocessing •••••••••••••••••••••
    #Find a way to no have to keep duplicating this.
    #find a way to pull this using a header file
    

    
def create_df(file): 
    
    #grabbing the Data from location
    data = pd.HDFStore(path + file)
    df = data['data']
    #df = df[:100]  #uncomment if want to open a smaller DF

 
    #create a column with the unix time.
    df['unixtime'] = df.index
    df = df.drop(('seq_num',0), axis=1)
    
    
    #•••••••••••••••••••••••••••••••••••••••  Clearning up the column names ••••••••••••••••••••••••••
    #A hack to remove the double row colum headings, and then rename it all as single row column headings
    df  = np.array(df)#.as_matrix()   #makes df into np array
    df  = pd.DataFrame(df) #To make it back into a DF
    

    
    
    #This is to rename the colums becasue the two lines above
    df  =df.rename(columns={ 0: "b1",  1: "bv1",  2: "a1",  3: "av1",\
                                                   4: "b2",  5: "bv2",  6: "a2",  7: "av2",\
                                                   8: "b3",  9: "bv3", 10: "a3", 11: "av3",\
                                                  12: "b4", 13: "bv4", 14: "a4", 15: "av4",\
                                                  16: "b5", 17: "bv5", 18: "a5", 19: "av5",\
                                                  20:"unix_time"})
    market_open_time = df['unix_time'][0]
    df["time_from_open"] = df["unix_time"] - market_open_time    
        
        
    #create a mide price.   
    df['midprice'] = ( df['b1'] + df['a1'] ) / 2   
    
    
    #Define a column that shows total volume on bid and ask side
    #df['total_bid_vol'] = (df['bv1'] + df['bv2'] + df['bv3'] + df['bv4'] + df['bv5'])
    #df['total_ask_vol'] = (df['av1'] + df['av2'] + df['av3'] + df['av4'] + df['av5'])
    
    return df



def create_3d_df(delta_str, T, N, Q):
    
    """
    Creats a 3d dataframe using a long string separated by commas, and 3 parameters
    T:
    N:
    Q:
    """
    delta = delta_str.split(",")            #Split the string by comma's
    delta = list(map(float,delta[:-1]))    #turns the strings into floats and removes the last space
    delta = np.array(delta)  
    
    df  = np.zeros((T,N,Q)) #3d array set
    
    counter  = 0
    for t in range(T):
        for j in range(N):
            for q in range(Q):
                df[t][j][q] =  delta[counter]
                counter+=1 
    return  df

