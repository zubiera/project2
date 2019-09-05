import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

#Import my own functions
import ST_Calibration_Fns as cb
from data_load import create_df






def Empirical_Inentsity(tick_size, depth, dt, day_string):
    
    """
    tick_size: this si the tick size in the data
    depth:     This is asking for how many levels we want an intensity for
    dt:        Time interval size
    day_string:the day in form 'dd' we want from out data set of SP500 in july2017
     
    """
    
    
    
    #******************Import the data needed for the below.
    df = create_df('data/1707' + day_string + '#ES.F.GLOB.1709#I#310#SBE##AB#t0.chi#top5.h5')    
    
    time  = np.array(df["time_from_open"])   #Time from market open in nano seconds
     
    ask    = np.array(df["a1"])     #best ask price
    avol   = np.array(df["av1"])    #volume at best ask
    avol2  = np.array(df["av2"])    #volume at level 2
    
    bid    = np.array(df["b1"])
    bvol   = np.array(df["bv1"])
    bvol2  = np.array(df["bv2"])
    
    
    
    #Creat vectors to stor the intensity at each distance from reference price.
    Optimistic_ask_intensity  = [0]*depth
    Pessimistic_ask_intensity = [0]*depth
    Optimistic_bid_intensity  = [0]*depth
    Pessimistic_bid_intensity = [0]*depth
    
    for i in range(0,depth):
        Optimistic_ask_intensity[i]  = cb.Ask_Opti(time, ask,avol,avol2,i*tick_size,dt)
        Pessimistic_ask_intensity[i] = cb.Ask_Pessi(time,ask,avol,avol2,i*tick_size,dt)
        Optimistic_bid_intensity[i]  = cb.Bid_Opti(time, bid,bvol,bvol2,i*tick_size,dt)
        Pessimistic_bid_intensity[i] = cb.Bid_Pessi(time,bid,bvol,bvol2,i*tick_size,dt)
    
    #Convert the Above lists into numoy array for math calculations
    Opti_a  = np.array(Optimistic_ask_intensity)
    Pessi_a = np.array(Pessimistic_ask_intensity)
    
    Opti_b  = np.array(Optimistic_bid_intensity)
    Pessi_b = np.array(Pessimistic_bid_intensity)
        
        
    return     (Opti_a, Pessi_a, Opti_b, Pessi_b)






def lstsq_regression(intensity_vector,tick_size,delta0_on_off=0):
    
    
    """
    This function does least square regression on a vector of values.
    intensity_vector : here input the vecor that is returned by def Empirical_Inentsity()
    tick_size:          tick size
    delta0_on_off :     chooses where delta starts from, i.e whether we include delta=0 
                        or not
    """
    
    y = np.log(intensity_vector)
    x = np.arange(0,(len(y))*tick_size,tick_size)
    
    #This cuts off the first intensity at delta=0 comment out if not needed.
    y = y[delta0_on_off:]
    x = x[delta0_on_off:]
    
    #We can rewrite the line equation as y = Wp, where W = [[x 1]] and p = [[m], [c]]. Now use lstsq to solve for p:
    W = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(W, y)[0]
    
    #Paramters from exponential intensity function
    A =  np.exp(c)
    k = -m
    
    return A,k



#•••••••••••••••••••••••••••• Running a loop to do regression on multiple days intensity ••••••••••••••••••••••••••••••••••

tick  = 25
depth = 6
dt    = 60e9
#day_list   = ["07"]
day_list = ['03','04','05','06','07','10','11','12','13','14','17',\
                  '18', '19','20','21','24','25','26','27','28','31']
vector_of_As = [0]*len(day_list)
vector_of_ks = [0]*len(day_list)



Parameters   = np.array([['Day','A_oa','k_oa','A_pa','k_pa','A_ob','k_ob','A_pb','k_pb']])

for i in range(len(day_list)):
    Opti_a, Pessi_a, Opti_b, Pessi_b = Empirical_Inentsity(tick,depth,dt,day_list[i])
    
    #Run a regression fit on the above outputs
    A_oa, k_oa = lstsq_regression(Opti_a,tick,1)
    A_pa, k_pa = lstsq_regression(Pessi_a,tick,1)
    
    A_ob, k_ob = lstsq_regression(Opti_b,tick,1)
    A_pb, k_pb = lstsq_regression(Pessi_b,tick,1)
    
    
    Parameters  = np.concatenate((Parameters, np.array([[day_list[i], float(A_oa), k_oa, A_pa, k_pa,\
                                                        A_ob, k_ob, A_pb, k_pb]])))

    
    
import csv 
path = '/Users/amitarfan1/Documents/Phd/3yr/Excel_Files/'
with  open(path +'intensity_param_30sec.csv', 'w',newline='') as f:
    thewriter =  csv.writer(f)
    thewriter.writerow(['Day','A_oa','k_oa','A_pa','k_pa','A_ob','k_ob','A_pb','k_pb'])
    
    rows = Parameters.shape[0]
    for i in range (1,rows):
        thewriter.writerow([Parameters[i][0],Parameters[i][1],Parameters[i][2],Parameters[i][3],\
                            Parameters[i][4],Parameters[i][5],Parameters[i][6],Parameters[i][7],\
                            Parameters[i][8]])        




    
    
    
#************************************** Creating Intensity Graphs   **********************************
day_list = ['03','04','05','06','07','10','11','12','13','14']
 
Opti_a_average = [0]*6 
Pessi_a_average = [0]*6
Opti_b_average = [0]*6
Pessi_b_average = [0]*6   
delta = np.arange(0,6*25,25)
for i in range(len(day_list)):   
    Opti_a, Pessi_a, Opti_b, Pessi_b = Empirical_Inentsity(tick,6,60e9,day_list[i])
    Opti_a_average  += Opti_a 
    Pessi_a_average += Pessi_a
    Opti_b_average  += Opti_b
    Pessi_b_average += Pessi_b
    
    plt.subplot(2,2,1)
    plt.plot(delta[1:],Opti_a[1:],linewidth=0.5)
    
    plt.subplot(2,2,2)
    plt.plot(delta[1:],Opti_b[1:],linewidth=0.5)
    
    plt.subplot(2,2,3)
    plt.plot(delta,Pessi_b,linewidth=0.5)
    
    plt.subplot(2,2,4)
    plt.plot(delta,Pessi_b,linewidth=0.5)

plt.subplot(2,2,1)
plt.plot(delta[1:],Opti_a_average[1:]/len(day_list), color = 'k',linewidth=2)
plt.title('Optimistic Ask Intensity w/out $\delta=0$')
plt.ylabel('$\Lambda(\delta^a)$ per min')
plt.xlabel('$\delta^a$')

plt.subplot(2,2,2)
plt.plot(delta[1:],Opti_b_average[1:]/len(day_list), color = 'k',linewidth=2)
plt.title('Optimistic Bid Intensity w/out $\delta=0$')
plt.ylabel('$\Lambda(\delta^b)$ per min')
plt.xlabel('$\delta^b$')

plt.subplot(2,2,3)
plt.plot(delta,Pessi_a_average/len(day_list), color = 'k',linewidth=2)
plt.title('Pessimistic Ask Intensity')
plt.ylabel('$\Lambda(\delta^a)$ per min')
plt.xlabel('$\delta^a$')

plt.subplot(2,2,4)
plt.plot(delta,Pessi_b_average/len(day_list), color = 'k',linewidth=2)
plt.title('Pessimistic Bid Intensity')
plt.ylabel('$\Lambda(\delta^b)$ per min')
plt.xlabel('$\delta^b$')

plt.tight_layout()



#Plotting only optimistic case with vs without delta=0
day_list = ['03','04','05','06','07','10','11','12','13','14']
 
Opti_a_average = [0]*6 
Pessi_a_average = [0]*6
Opti_b_average = [0]*6
Pessi_b_average = [0]*6   
delta = np.arange(0,6*25,25)
for i in range(len(day_list)):   
    Opti_a, Pessi_a, Opti_b, Pessi_b = Empirical_Inentsity(tick,6,60e9,day_list[i])
    
    Opti_a_average  += Opti_a 
    Opti_b_average  += Opti_b
    
    #plotting with delta=0
    plt.subplot(1,2,1)
    plt.plot(delta[0:],Opti_a[0:],linewidth=0.5)
    
    #plotting without delta=0
    plt.subplot(1,2,2)
    plt.plot(delta[1:],Opti_a[1:],linewidth=0.5)

#Adding on the mean plot with labels and cosmetic shit
plt.subplot(1,2,1)
plt.plot(delta[0:],Opti_a_average[0:]/len(day_list), color = 'k',linewidth=2)
plt.title('Optimistic Ask Intensity')
plt.ylabel('$\Lambda(\delta^a)$ per min')
plt.xlabel('$\delta^a$')

plt.subplot(1,2,2)
plt.plot(delta[1:],Opti_a_average[1:]/len(day_list), color = 'k',linewidth=2)
plt.title('Optimistic Ask Intensity')
plt.ylabel('$\Lambda(\delta^a)$ per min')
plt.xlabel('$\delta^a$')

plt.tight_layout()




#*********************** creat a plot of the fit and choosing what to plot it with************
A,k = lstsq_regression(Opti_a_average/10,tick,1)

x =  np.arange(25,126)
fit = A*np.exp(-k*x)
plt.subplot(1,2,2)
plt.plot(delta[1:],Opti_a[1:])
plt.plot(x,fit,color='red')

#plot including delta = 0
A,k = lstsq_regression(Opti_a_average/10,tick)

x =  np.arange(0,126)
fit = A*np.exp(-k*x)
plt.subplot(1,2,1)
#plt.plot(delta[0:],Opti_a[0:])
plt.plot(x,fit,color='red')
