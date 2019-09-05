import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt




simulation_df = pd.DataFrame()



mu    = 0.814627081824764
sigma = 12.3
 
No_sims = 1000

T  = 1380
N  = 1380
dt = T/N 



for j in range(No_sims):
    S = np.zeros(1380)
    S[0] =244650
    
    for i in range(1,T):
        S[i] = S[i-1] + mu*dt + sigma*np.random.normal(0, 1)
    simulation_df[j] = S   
    plt.plot(S)
 
    
average = simulation_df.mean(axis=1) 
variance_check = np.mean(np.sqrt(np.var(simulation_df)/1380))  
np.var(average)
np.sqrt(np.var(average)/1380)
plt.plot(ask["14"],color="black")