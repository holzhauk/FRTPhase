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

MFPTpath = Path("../../SimData/MFPTs/NewbySchwemmerAntirotating.h5")
Isochronepath = Path("../../Isochrones/" + MODEL_NAME + \
                     "/AntirotatingNewbySchwemmerIsochroneSet.h5")
isochroneset = IsochroneSet(Isochronepath, MODEL_NAME)
mfptset = MFPTSet(MFPTpath, MODEL_NAME)

latex_textwidth_in = 6.50127
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(w=latex_textwidth_in, h=0.5*latex_textwidth_in)

evenly_spaced_interval = np.linspace(0.0, 1.0, num=7)
colors = [cm.Set1(x) for x in evenly_spaced_interval]

Esize = 0


for i, I in enumerate(mfptset):
    (pset, rhos, phis) = isochroneset.__get_Isochrone__(I.key)
    ax1.plot(I.Rho_init, I.MFPT, '.', color = colors[i],\
             label="D = " + str(pset["D"][0]))
    ax2.plot(rhos*np.cos(phis), rhos*np.sin(phis), color = colors[i])
    Esize = I.EnsembleSize
    
ax2.axis("equal")
ax2.set_title("Isochronen")
ax1.legend(loc=3)
ax1.set_xlabel(r"$\rho$")
ax1.set_ylabel(r"MFPT")
fig.suptitle(r"N = $" + str(Esize[0]) + "$")

del isochroneset
del mfptset

plt.savefig("../../../masterthesis/Plots/NewbySchwemmerMFPTsIsochronen.pgf")
    