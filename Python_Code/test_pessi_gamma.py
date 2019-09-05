import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


from data_load import create_df

gamma =  ["0.001",\
          "0.002",\
          "0.003", \
          "0.004", \
          "0.005",\
          "0.006",\
          "0.007",\
          "0.008",\
          "0.009",\
          "0.01",\
          "0.011",\
          "0.012", \
          "0.013", \
          "0.014",\
          "0.015",\
          "0.016",\
          "0.017",\
          "0.018",\
          "0.019",\
          "0.02"
          ]

output = np.array([['gamma','mean_trades', "variance_trades",'mean_Cash',"variance_cash"]])
for g in range (len(gamma)):
    #enter the days of the files to get results for.
    data_file_list = ['03','04','05','06','07','10','11','12','13','14','17',\
                      '18', '19','20','21','24','25','26','27','28','31']
    
    
    Max_inv=0    #P0aramteer to check max_inv 
    
    PnL  = np.array([['No_Trades','Cash']])
    
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
        folder       = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/average_folder/3_to17_jul/"
        ask_deltas = "delta_a_0.26468_10.99_"+gamma[g]+"_1.19156_0.04535_ave_opti_ask.csv"
        bid_deltas = "delta_b_0.26468_10.99_"+gamma[g]+"_1.19156_0.04535_ave_opti_ask.csv"
        #turn the data sheets into dataframes
        df_a    =   pd.read_csv(folder + ask_deltas ) #,header=None
        df_b    =   pd.read_csv(folder + bid_deltas ) #,header=None
        
        #Overwrite the Strategy Dataframe swith chosen numbers
    #    for i in range(-10,11):
    #        df_a["q="+str(i)] = 0
    #        df_b["q="+str(i)] = 0  
            
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
            
        PnL = np.concatenate((PnL, np.array([[len(event_tracker),cash]])))
    
        #convnert event_tracker to a dataframe
        event_tracker = pd.DataFrame(event_tracker)
        event_tracker.columns = event_tracker.iloc[0]  #make row zero column headings
        event_tracker = event_tracker.reindex(event_tracker.index.drop(0))  #drop row zero
        event_tracker['inventory'] = event_tracker['inventory'].astype('float64') 
        Max_inv = max(Max_inv, max(abs(min(event_tracker["inventory"])), abs(max(event_tracker["inventory"])) ) )
    
    PnL  = pd.DataFrame(PnL)
    PnL.columns = PnL.iloc[0] 
    PnL = PnL.reindex(PnL.index.drop(0))
    PnL['Cash'] = PnL['Cash'].astype('float64')
    PnL['No_Trades'] = PnL['No_Trades'].astype('float64')
        
    Mean = np.mean(PnL["Cash"])
    SD = np.var(PnL["Cash"])
    Mean_trades = np.mean(PnL["No_Trades"])
    SD_trades = np.var(PnL["No_Trades"])
        
    output = np.concatenate((output, np.array([[gamma[g],Mean_trades,SD_trades,Mean,SD]])))


    
output  = pd.DataFrame(output)
output.columns = output.iloc[0] 
output = output.reindex(output.index.drop(0))
output['gamma'] = output['gamma'].astype('float64')
output['mean_Cash'] = output['mean_Cash'].astype('float64')
output['variance_cash'] = output['variance_cash'].astype('float64')
output['mean_trades'] = output['mean_trades'].astype('float64')
output['variance_trades'] = output['variance_trades'].astype('float64')
    
#Opti_output  = pd.DataFrame(Opti_output)
#Opti_output.columns = Opti_output.iloc[0] 
#Opti_output = Opti_output.reindex(Opti_output.index.drop(0))
#Opti_output['gamma'] = Opti_output['gamma'].astype('float64')
#Opti_output['mean_Cash'] = Opti_output['mean_Cash'].astype('float64')
#Opti_output['variance_cash'] = Opti_output['variance_cash'].astype('float64')
#Opti_output['mean_trades'] = Opti_output['mean_trades'].astype('float64')
#Opti_output['variance_trades'] = Opti_output['variance_trades'].astype('float64')



#plotting some graphs of vairance. and    

plt.subplot(2,2,3)
plt.plot(output['gamma'],output['mean_Cash'],linewidth=0.5, color = 'blue')
plt.xlabel("Risk Averasion $\gamma$")
plt.ylabel("Average Profit (cents)")
plt.title("Pessimistic")


plt.subplot(2,2,4)
plt.plot(output['gamma'],output['variance_cash'],linewidth=0.5, color = 'blue')
plt.xlabel("Risk Averasion $\gamma$")
plt.ylabel("Variance of Profit")
plt.title("Pessimistic")


plt.subplot(2,2,1)
plt.plot(Opti_output['gamma'],Opti_output['mean_Cash'],linewidth=0.5, color = 'blue')
plt.xlabel("Risk Averasion $\gamma$")
plt.ylabel("Average Profit (cents)")
plt.title("Optimistic")

plt.subplot(2,2,2)
plt.plot(Opti_output['gamma'],Opti_output['variance_cash'],linewidth=0.5, color = 'blue')
plt.xlabel("Risk Averasion $\gamma$")
plt.ylabel("Variance of Profit")
plt.title("Optimisitc")


plt.tight_layout()




plt.plot(output['gamma'],output['mean_trades'])
plt.plot(output['gamma'],output['variance_trades'])



plt.plot(Opti_output['gamma'],Opti_output['mean_Cash'])
plt.plot(Opti_output['gamma'],Opti_output['variance_cash'])
plt.plot(Opti_output['gamma'],Opti_output['mean_trades'])
plt.plot(Opti_output['gamma'],Opti_output['variance_trades'])




def maxDiff(arr, arr_size): 
    max_diff = arr[0] - arr[1] 
    max_element = arr[0]
      
    for i in range( 1, arr_size ): 
        if (arr[i] > max_element): 
            max_element = arr[i] 
            
        if (max_element - arr[i] > max_diff): 
            max_diff = max_element - arr[i] 
      

    return max_diff 

arr = [1, 10, 6, 101, 100,99,92,91.5] 
size = len(arr) 
print ("Maximum difference is",  
        maxDiff(arr, size)) 

        