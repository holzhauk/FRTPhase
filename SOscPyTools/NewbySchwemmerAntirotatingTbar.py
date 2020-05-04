#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:19:24 2020

@author: konstantin
"""
import matplotlib

matplotlib.use("pdf")

MODEL_NAME = "NewbySchwemmer"
from matplotlib import pyplot as plt
from SOscFilePy import *
from pathlib import Path
from scipy.integrate import quad


# parameters
P = {
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
        }

# radial distribution function
def RDF(r, D):
    return np.exp( (P["gamma"] / (2*D)) * r**2 * (1 - r**2 / 2) )

def OmegaBar(D, rm, rp):
    
    i1 = lambda r: r*(1-r)**2*RDF(r, D)
    i2 = lambda r: r*RDF(r, D)
    I1, errI1 = quad(i1, rm, rp)
    I2, errI2 = quad(i2, rm, rp)
    
    omegabar = P["omega"] + P["omega"]*P["gamma"]*P["c"]*(I1 / I2)
    # gaussian error propagation formula
    err = ((P["omega"]*P["gamma"]*P["c"]) / I2) * errI1 + \
            P["omega"]*P["gamma"]*P["c"]*(I1 / I2**2) * errI2
    
    return (omegabar, err)

MFPTpath = Path("../../SimData/MFPTs/NewbySchwemmerAntirotating.h5")
Isochronepath = Path("../../Isochrones/" + MODEL_NAME + \
                     "/AntirotatingNewbySchwemmerIsochroneSet.h5")
ISet = IsochroneSet(Isochronepath, MODEL_NAME)
FPTSet = MFPTSet(MFPTpath, MODEL_NAME)

noSets = sum(1 for _ in FPTSet)
Ds = np.zeros((noSets, ))
Tbars = np.zeros((noSets, ))


rho_min = 0.0
rho_max = 0.0
for i, I in enumerate(FPTSet, start=0):
    (pSet, IsoRho, IsoPhi) = ISet.__get_Isochrone__(I.key)
    rho_min = np.amin(IsoRho)
    rho_max = np.amax(IsoRho)
    Ds[i] = pSet["D"]
    Tbars[i] = np.mean(I.Tbar)
    
    
D = np.linspace(Ds[1], 100, num=1000)
obars_D0 = np.ones(D.shape)
obars_D_inf = (P["omega"] + P["omega"]*P["gamma"]*P["c"]*\
         (rho_max**2*(rho_max**2/2 - 4*rho_max/3 + 1) - \
          rho_min**2*(rho_min**2/2 - 4*rho_min/3 + 1)) / \
             (rho_max**2 - rho_min**2))*np.ones(D.shape)
obars = np.zeros(D.shape)
for i in range(len(D)):
    obars[i] = OmegaBar(D[i], rho_min, rho_max)[0]
    
latex_textwidth_in = 6.50127
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=latex_textwidth_in, h=0.5*latex_textwidth_in)
    
ax.plot(np.log(Ds[1:]), 2*np.pi / Tbars[1:], "x")
ax.plot(np.log(D), obars, "k")
ax.plot(np.log(D), obars_D_inf, "k--")
ax.plot(np.log(D), obars_D0, "k:")
ax.set_title(r"$T=100.0, N=100, dt=0.0001$")
ax.set_xlabel(r"$\log{D}$")
ax.set_ylabel(r"$\overline{\omega}$")
ax.legend(["Simulation", "Analytisch", r"$D \rightarrow \infty$",\
            r"$D \rightarrow 0$"])
    
del FPTSet
del ISet

plt.savefig("OmegaBar.pdf")