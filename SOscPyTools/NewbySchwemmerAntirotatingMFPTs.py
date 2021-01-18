#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 18:10:55 2020

@author: konstantin
"""

import matplotlib

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

MODEL_NAME = "NewbySchwemmer"

from matplotlib import pyplot as plt
from matplotlib import cm

from SOscFilePy import *
from pathlib import Path
from SimConfigPy import *

if __name__ == "__main__":

    if (len(sys.argv)) != 2:
        print("Usage: " + sys.argv[0] + " <SimConfigFile.json>")
        os._exit(os.EX_IOERR)

    simConfig = SimConfig()
    simConfig.load_from_file(Path(sys.argv[1]))
    print("Input: " + str(simConfig.Paths["In"].resolve()))
    print("Output: " + str(simConfig.Paths["Out"].resolve()))

    isochroneset = IsochroneSet(simConfig.Paths["In"].resolve(), \
            simConfig.ModelName)
    mfptset = MFPTSet(simConfig.Paths["Out"].resolve(), \
            simConfig.ModelName)
    
    latex_textwidth_in = 6.50127
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(w=latex_textwidth_in, h=0.5*latex_textwidth_in)
    
    evenly_spaced_interval = np.linspace(0.0, 1.0, num=7)
    colors = [cm.Set1(x) for x in evenly_spaced_interval]
    
    for i, I in enumerate(mfptset):
        if i != 7 :
            (pset, rhos, phis) = isochroneset.__get_Isochrone__(I.key)
            ax1.plot(I.Rho_init, I.MFPT, '.', color = colors[i],\
                     label="D = " + str(pset["D"][0])) 
            ax2.plot(rhos*np.cos(phis), rhos*np.sin(phis), color = colors[i])
        
    phis = np.linspace(0.0, 2*np.pi, num=100)
    rm = np.amin(rhos)*np.ones(phis.shape)
    rp = np.amax(rhos)*np.ones(phis.shape)
    ax2.plot(rm*np.cos(phis), rm*np.sin(phis), ":k")
    ax2.plot(rp*np.cos(phis), rp*np.sin(phis), "--k")
 

    ax2.axis("equal")
    ax2.set_title("Isochronen")
    ax1.legend(loc=3)
    ax1.set_xlabel(r"$\rho$")
    ax1.set_ylabel(r"MFPT")
    fig.suptitle(r"N = $" + str(simConfig.Simulation["Ensemble Size"]) + \
            "$, dt = $" + str(simConfig.Simulation["dt"]) + "$")

    plt.savefig("NewbySchwemmerMFPTsIsochronen.pgf")
    plt.savefig("NewbySchwemmerMFPTsIsochronen.pdf")
    
    del isochroneset
    del mfptset 
 
