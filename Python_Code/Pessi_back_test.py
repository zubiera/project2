import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


from data_load import create_df

#enter the days of the files to get results for.
data_file_list = ['17','18', '19','20','21','24','25','26','27','28','31']

Max_inv=0    #P0aramteer to check max_inv 

PnL  = np.array([['No_Trades','Cash',"maxinv"]])

for  file in  data_file_list:
    #**************************Process of the raw data**************************
    df = create_df('data/1707'+str(file)+'#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')
    df_small = df[:1000]
    
                
    time      = np.array(df["time_from_open"])  
    ask_price = np.array(df["a1"])
    ask_vol   = np.array(df["av1"])
    ask_vol2   = np.array(df["av2"])
    
    bid_price = np.array(df["b1"])
    bid_vol   = np.array(df["bv1"])
    bid_vol2   = np.array(df["bv2"])
    
    
    #***************************** PDE Results *************************************
    #This does not really need to be inside the loop!!!
    #Run the PDE (check file to make sure the right code is in place)
    #import os
    #a  = os.system("/Users/amitarfan1/Documents/Phd/3yr/NAG\ Code/buildXCode/Debug/nagtest2.out 10000000000000 0.0 0.3")
    
    
    #load PDE results
    folder       = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/vary_vol/"



    ask_deltas = "sigma_4seg_10d_ave_Pessi_delta_a.csv"
    bid_deltas = "sigma_4seg_10d_ave_Pessi_delta_b.csv"
    #turn the data sheets into dataframes
    df_a    =   pd.read_csv(folder + ask_deltas ) #,header=None
    df_b    =   pd.read_csv(folder + bid_deltas ) #,header=None
    
    #Overwrite the Strategy Dataframe swith chosen numbers
    for i in range(-10,11):
        df_a["q="+str(i)] = 25
        df_b["q="+str(i)] = 25  
        
    # *******************function to round the quotes to nearest tick***************
    def myround(x, base=25):
        return int(base * round(float(x)/base))
   
    for i in range(-10,11):
        df_a["q="+str(i)] = df_a["q="+str(i)].apply(lambda x: myround(x, base=25))
        df_b["q="+str(i)] = df_b["q="+str(i)].apply(lambda x: myround(x, base=25))
    df_a["q=-10"] = 500    #setting the Q_max so large that it doesnt go over inv bounds
    df_b["q=10"]  = 500
    #***************************** Run Strategy **************************************
    
    #The algo posts quotes every dt.  (60billion nano-seconds = 1 minute)
    dt = 60e9                     
    
    
    
    
    #initlisise paramters
    q   = 0      #Initial inventory 
    x   = 0      #Initial cashflow
    S_a = [0]    #vector tracking quoted ask price at each dt
    S_b = [0]
    idx = [0]    #vector tracking index locations
    event_tracker  = np.array([['minute','time_idx','event','inventory','cash','ex_price','delta']])
    quote_tracker  = np.array([['time_idx','ask_quote','bid_quote']])
    
    
    i=0
    while np.searchsorted(time,dt*(i+1))<len(time):
        
        
        # can manually overwride this to determin fixed strategy vs DataFrame strategy
        delta_a = df_a.at[i,'q=' + str(q)]
        delta_b = df_b.at[i,'q=' + str(q)]
        
        #setting the Limit order quotes
        LO_a = ask_price[np.searchsorted(time,dt*i)] + delta_a
        LO_b = bid_price[np.searchsorted(time,dt*i)] - delta_b
        
        #keeping track of what is quoted at what time.    
        quote_tracker = np.concatenate((quote_tracker, np.array([[np.searchsorted(time,dt*i), LO_a, LO_b]])))
        
        is_b_hit = False
        is_a_hit = False
        
        
        
        #iterate within each timestep
        for j in range( np.searchsorted(time,dt*i), np.searchsorted(time,dt*(i+1)) ):
            
        #************************below checks for ask side executions****************************
            #check if delta implies a market order  (if so, insntant execution)
            if (delta_a < 0) and (is_a_hit == False):
                x = x + bid_price[j]  #execute as market order at best bid
                q = q - 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'MO_sell', int(q), float(x), bid_price[j],delta_a ]])))
                is_a_hit = True
    

            if LO_a <  ask_price[j+1] and (is_a_hit == False):
                x = x + LO_a
                q = q - 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'a at MO walking over our price', int(q), float(x), LO_a, delta_a ]])))
                is_a_hit = True
                    
                    

                
        #*************Below checks for  bid sides executions****************************************
            #check if delta implies a market order
            if (delta_b < 0) and (is_b_hit == False):
                x = x - ask_price[j]  #execute as market order at best bid
                q = q + 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'MO_buy', int(q), float(x), ask_price[j], delta_b ]])))
                is_b_hit = True
                
            if LO_b >  bid_price[j+1] and (is_b_hit == False):
                x = x - LO_b
                q = q + 1
                event_tracker = np.concatenate((event_tracker, np.array([[i,j, 'b at MO walking over our price', int(q), float(x), LO_b, delta_b ]])))
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
        
    

    #convnert event_tracker to a dataframe
    event_tracker = pd.DataFrame(event_tracker)
    event_tracker.columns = event_tracker.iloc[0]  #make row zero column headings
    event_tracker = event_tracker.reindex(event_tracker.index.drop(0))  #drop row zero
    event_tracker['inventory'] = event_tracker['inventory'].astype('float64') 
    Max_inv = max(Max_inv, max(abs(min(event_tracker["inventory"])), abs(max(event_tracker["inventory"])) ) )
    
    PnL = np.concatenate((PnL, np.array([[len(event_tracker),cash, Max_inv]])))

#print(bid_vol[-1]," ", ask_vol[-1])

PnL[-1][-1]







