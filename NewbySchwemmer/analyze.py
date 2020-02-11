import sys
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from NewbySchwemmer import polar2XY

import scipy.integrate as integrate

def Tbarnum(P, D, rhom, rhop):
    
    P["D"] = D
    
    gbar = lambda r: -(P["gamma"]/P["D"])*r*(r**2 - 1) + 1/r
    f = lambda r: P["omega"] + P["omega"]*P["gamma"]*P["c"]*(1-r)**2
    
    integrand1 = lambda x: np.exp(-integrate.quad(gbar, x, rhop)[0])
    integrand2 = lambda x: f(x)*integrand1(x)
    
    Nominator = 2*np.pi*integrate.quad(integrand1, rhom, rhop)[0]
    Denominator = integrate.quad(integrand2, rhom, rhop)[0]
    
    return Nominator / Denominator


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python" + sys.argv[0] +  "<Filepath>")
        exit()

    FILEPATH = Path(sys.argv[1])
    FILENAMES = glob.glob(str(FILEPATH / "*.h5"))

    Diff = np.zeros(len(FILENAMES))
    Tbar = np.zeros(len(FILENAMES))
    P = 0

    for i in range(len(FILENAMES)):

        store = pd.HDFStore(FILENAMES[i])
        P = store.get("parameters")
        DFkeys = store.keys()[1::]
        Emeans = pd.Series()

        # Ensemble mean
        for k in DFkeys:
            D = store.get(k)
            #print(D)
            Phis = D[D.columns[pd.Series(D.columns).str.startswith('Phi')]]
            #print(Phis.iloc[-1].mean())
            Emeans = Emeans.append(pd.Series([Phis.iloc[-1].mean()]),\
                    ignore_index=True)
        
        Diff[i] = P["D"]
        Tbar[i] = 2*np.pi*P["T"] / Emeans.mean()
        print("Mean Value for D = %f: %f" %(P["D"], Tbar[i]))
        store.close()
    
    Ds = np.linspace(0.001, 1.5, num=20)
    Tbars = np.zeros(Ds.shape)
    for j in range(len(Tbars)):
        Tbars[j] = Tbarnum(P, Ds[i], 0.01, 1.5)

    plt.figure()
    #plt.plot(np.log(Diff), 2*np.pi / Tbar, '*')
    plt.plot(np.log(Ds), 2*np.pi / Tbars, '*')
    plt.show()
