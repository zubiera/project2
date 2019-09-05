from lmfit import minimize, Parameters

import numpy as np
import pandas as pd


S0    =  241025.0   # 100      
V0    =  1281.0734463283723 #4        
dt    =     1 #0.005
V = np.zeros(581)
S = np.zeros(581)
V[0] = V0
S[0] = S0

Vo = np.array(RVask60)

plotter = pd.DataFrame()
plotter["real"] = Vo
store = 0

storageParams = []
storageResid = []

def simulator(ps,simulation):
    theta = ps["theta"].value
    alpha = ps["alpha"].value
    eta   = ps["eta"].value
    #rho   = ps["rho"].value
    np.random.seed(simulation)
    rho = -0.05
    for i in range(1,len(V)):
        Zv   =  np.random.normal(0, 1)
        Zs   =  rho*Zv   + np.sqrt(1-rho*rho)*np.random.normal(0, 1)
        #V_plus   = max(0,V[i-1])
        V[i] = V[i-1] + theta*(alpha - max(0,V[i-1]))*dt + eta*np.sqrt(max(0,V[i-1]))*Zv*np.sqrt(dt)
        S[i] = S[i-1] + np.sqrt(max(0,V[i-1]))*Zs*np.sqrt(dt)
    if(store == 1):
        plotter[str(simulation)] = V
    return V
        


def minimizeSimu(ps):
    simulations = 10
    residuals = 0
    for i in range(simulations):
        linked = simulator(ps,i);
        diff = linked-Vo
        residuals += (diff*diff).sum()
    storageParams.append(ps)
    storageResid.append(residuals)
    print(residuals)
    return residuals

params = Parameters()
params.add('theta', value=0.0001, min=0.001, max=0.1)
params.add('alpha', value=5000, min=20, max=100000)
params.add('eta', value=5, min=2, max=50)
#params.add('rho', value=0.5, min=0.01, max=1.)

result = minimize(minimizeSimu, params, method='nelder')

store = 1
fd = minimizeSimu(result.params)
plotter.plot()
