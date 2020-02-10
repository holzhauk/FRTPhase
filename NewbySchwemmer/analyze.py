import sys
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from NewbySchwemmer import polar2XY

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

    plt.figure()
    plt.plot(Diff, Tbar, '*')
    plt.show()
