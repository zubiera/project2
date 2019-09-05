# -*- coding: utf-8 -*-
"""
This file runs the PDE for different paramters
"""
import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

#Parameters (argv)
spread = "150000"
mu     = "0"
sigma  = "13.06"
gamma  = "0.01 "
A      = "28.69"
k      = "0.0272"
case   = "opti"

#These lists run the PDE using the calibration associated with each date for a given dt.
spread = "15000"  

#day    = "Comapre_w_stochastic"
#mu     =  "0"
#sigma  = "4"
#gamma  =  "0.01"  
#A      = "140"     
#k      = "1.5"
#case   = "base"
#sigma_1= "4"            
#sigma_2= "4"


day    = "2_seg"
mu     =  "0"
sigma  = "4"
gamma  =  "0.001"  
A      = "1.212"     
k      = "0.0448"
case   = "opti"
sigma_1= "5.3"            
sigma_2= "13.8"

        
for i in range(len(A)):
   parameters = spread + " " + mu + " " + sigma + " " + gamma + " " + A + " " + k + " " +  case + " "+ day + " " + sigma_1 + " " + sigma_2  # therese are the argv[]
   file       = "/Users/amitarfan1/Documents/Phd/3yr/NAG\ Code/buildXCode/Debug/nagtest2.out "
   os.system(file + parameters) 
   print(parameters)






#**************************  Opening PDe results for comparisons *********************

#comaparing the pessimisssmitc v optimistc with drift.
main  = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/minmax_sigma/"
main2 = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/vary_vol/gamm001/"


case1      = "vol_sigma_10d_ave_Opti_delta_a.csv"
df_case1   = pd.read_csv(main2 + case1)
case2      = "const_sigma_10d_ave_Opti_delta_a.csv"
df_case2   = pd.read_csv(main2 + case2)     
case3      = "fixed_sigma_3_delta_a.csv"
df_case3   = pd.read_csv(main2 + case3)
case4     = "fixed_sigma_4_delta_a.csv"
df_case4   = pd.read_csv(main2 + case4)  
case5     = "fixed_sigma_5_delta_a.csv"
df_case5   = pd.read_csv(main2 + case5)  




for i in range(-5,6):
    
    plt.plot(df_case1["Time_Step"],df_case1["q="+str(i)],linewidth = 0.6, color = 'red')
    plt.plot(df_case2["Time_Step"],df_case2["q="+str(i)],linewidth = 0.6, color = 'blue')
    plt.plot(df_case3["Time_Step"],df_case3["q="+str(i)],linewidth = 0.6, color = 'green')
    plt.plot(df_case4["Time_Step"],df_case4["q="+str(i)],linewidth = 0.6, color = 'black')
    plt.plot(df_case5["Time_Step"],df_case5["q="+str(i)],linewidth = 0.6, color = 'pink')

  

plt.title("$\sigma_{min}$ vs $\sigma_{max}$ ")
plt.xlabel("Time (minutes)")
plt.ylabel("$\delta$")





def myround(x, base=25):
    return int(base * round(float(x)/base))

#round the solutins to nearest 25   
for i in range(-10,11):
    df_opti["q="+str(i)] = df_opti["q="+str(i)].apply(lambda x: myround(x, base=25))
    df_pessi["q="+str(i)] = df_pessi["q="+str(i)].apply(lambda x: myround(x, base=25))


for i in range(-9,10):
    plt.plot(df_opti["Time_Step"],df_opti["q="+str(i)],linewidth = 0.6, color = 'blue')
    plt.plot(df_pessi["Time_Step"],df_pessi["q="+str(i)],linewidth = 0.6, color='grey')

plt.title("Optimistic vs Pessimsitic strategy")
plt.xlabel("Time (minutes)")
plt.ylabel("$\delta$")






#just to plot an illustration of the optimal deltas
tester = pd.read_csv("/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/Tester_b.csv")


tester = pd.read_csv("/Users/amitarfan1/Documents/Phd/3yr/Code/Output/back_test/excl_delta0/vary_vol/07_opti_delta_a19augtest.csv")
fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(-6,7,2):
    plt.plot(tester["Time_Step"],tester["q="+str(i)],linewidth = 0.6)

import matplotlib.patches as mpatches
blue = mpatches.Patch(color='blue', label='q=-6')
orange = mpatches.Patch(color='orange', label='q=-4')
green = mpatches.Patch(color='green', label='q=-2')
red = mpatches.Patch(color='red', label='q=0')
purp = mpatches.Patch(color='purple', label='q=2')
brown = mpatches.Patch(color='brown', label='q=4')
pink = mpatches.Patch(color='pink', label='q=6')
ax.legend(handles=[blue,orange, green, red, purp, brown, pink])
ax.set_xlabel('time')
ax.set_ylabel('$\delta^b$')
