import sys
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from NewbySchwemmer import polar2XY
from Tbar import Phidotbarnum

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python" + sys.argv[0] +  "<Filepath>")
        exit()

    FILEPATH = Path(sys.argv[1])
    FILENAMES = glob.glob(str(FILEPATH / "*.h5"))

    Diff = np.zeros(len(FILENAMES))
    Tbar = np.zeros(len(FILENAMES))

    for i in range(len(FILENAMES)):

        store = pd.HDFStore(FILENAMES[i])
        P = store.get("parameters")
        DFkeys = store.keys()[1::]
        Emeans = pd.Series([])

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
    
    P = {
        "D": 0.198,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
        }    

    Ds = np.linspace(0.01, 1.0, num=100)
    PhiDots = np.zeros(Ds.shape)
    for m in range(len(Ds)):
        PhiDots[m] = Phidotbarnum(P, Ds[m], 0.01, 1.5)

    fig, ax = plt.subplots()
    ax.set_title(r"$N = 20$, $\rho_{-} = 0.01$, $\rho_{+} = 1.5$")
    
    ax.plot(np.log(Diff[1::]), 2*np.pi / Tbar[1::], '*')
    ax.plot(np.log(Ds), PhiDots, 'r')
    
    ax.set_xlabel(r"$log(D)$")
    ax.set_ylabel(r"$2\,\pi / \overline{T}$")
    
    ax.spines["right"].set_position("zero")
    ax.spines["left"].set_color("none")
    ax.yaxis.tick_right()
    ax.spines["top"].set_position("zero")
    ax.spines["bottom"].set_color("none")
    ax.xaxis.tick_top()

    ax.legend(["in silico", "numeric"], frameon = False, loc='lower left')

    fig.savefig(str(FILEPATH / "phidotbar.pdf"))
