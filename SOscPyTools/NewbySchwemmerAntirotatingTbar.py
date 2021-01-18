#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:19:24 2020

@author: konstantin
"""
import matplotlib

matplotlib.use("pdf")

import os, sys
from pathlib import Path

MODEL_NAME = "NewbySchwemmer"
from matplotlib import pyplot as plt
from SOscFilePy import *
from pathlib import Path
from scipy.integrate import quad

from SimConfigPy import *


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


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <SimConfigFile>.json")
        os._exit(os.EX_IOERR)

    simConfig = SimConfig()
    simConfig.load_from_file(Path(sys.argv[1]))
    print("Input: " + str(simConfig.Paths["In"].resolve()))
    print("Output: " + str(simConfig.Paths["Out"].resolve()))

    ISet = IsochroneSet(simConfig.Paths["In"].resolve(), simConfig.ModelName)
    FPTSet = MFPTSet(simConfig.Paths["Out"].resolve(), simConfig.ModelName)

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

    del ISet
    del FPTSet
    
    
    D = np.linspace(Ds[1], 100, num=1000)
    obars_D0 = np.ones(D.shape)
    obars_D_inf = (P["omega"] + P["omega"]*P["gamma"]*P["c"]*\
            (rho_max**2*(rho_max**2/2 - 4*rho_max/3 + 1) - \
            rho_min**2*(rho_min**2/2 - 4*rho_min/3 + 1)) / \
            (rho_max**2 - rho_min**2))*np.ones(D.shape)
    obars = np.zeros(D.shape)
    for i in range(len(D)):
        obars[i] = OmegaBar(D[i], rho_min, rho_max)[0]

    HighNoiseSimConfig = SimConfig()
    HighNoiseSimConfig.load_from_file(Path("../configs/config_cluster_high_noise.json"));
    ISet = IsochroneSet(HighNoiseSimConfig.Paths["In"].resolve(), HighNoiseSimConfig.ModelName)
    FPTSet = MFPTSet(HighNoiseSimConfig.Paths["Out"].resolve(), HighNoiseSimConfig.ModelName)
    noSets = sum(1 for _ in FPTSet)
    HighNoiseTbars = np.zeros((noSets, ))
    HighNoiseVarTs = np.zeros((noSets, ))
    for i, I in enumerate(FPTSet, start=0):
        HighNoiseTbars[i] = np.mean(I.Tbar)

    VarTs = np.array([0.000190951, 0.000168459, 0.000190211, 0.00016833, 0.000178149, 0.000192491,\
                        0.000164856, 0.000174417, 0.000172473, 0.00017444, 0.000177766, 0.000176331, \
                        0.000172935, 0.000168798, 0.000178058, 0.00015719, 0.000162182, 0.000177715, \
                        0.000179366, 0.000160911], dtype=np.double)
    HighNoiseVarTs[0] = np.mean(VarTs)
    

    del ISet
    del FPTSet

    HighNoise2Tbars = np.array(np.mean(np.array([ -0.323264, -0.322959, -0.322945, -0.322403, -0.32411, -0.323167, \
                        -0.322175, -0.322915, -0.323712, -0.321833, -0.322604, -0.323077, -0.324169, \
                        -0.32533, -0.323273, -0.320134, -0.323755, -0.319503, -0.325758, -0.322986 ], dtype=np.double)))
    HighNoise2VarTs = np.array(np.mean(np.array([ 7.84507e-05, 4.77232e-05, 8.04538e-05, 7.61505e-05, 8.40962e-05, \
                        3.95456e-05, 7.22033e-05, 8.21001e-05, 5.89792e-05, 7.39548e-05, 7.54628e-05, 7.00665e-05, \
                        6.80643e-05, 0.000117195, 7.17906e-05, 4.88781e-05, 5.80263e-05, 6.82364e-05, 7.10541e-05, 8.87985e-05 ], dtype=np.double)))
        
    latex_textwidth_in = 6.50127
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=latex_textwidth_in, h=0.5*latex_textwidth_in) 

    ax.plot(np.log10(D), obars, "k")
    ax.plot(np.log10(D), obars_D_inf, "k--")
    ax.plot(np.log10(D), obars_D0, "k:")

    ax.plot(np.log10(Ds[1:6]), 2*np.pi / Tbars[1:6], "x")
    ax.plot(np.array([np.log10(100.0) - 0.05]), 2*np.pi / HighNoiseTbars, marker="x", color="orange")
    ax.errorbar(np.array([np.log10(100.0) - 0.05]), 2*np.pi / HighNoiseTbars,\
            yerr=(2*np.pi*np.sqrt(HighNoiseVarTs)/(HighNoiseTbars**2)), mec="orange")
    ax.plot(np.array([np.log10(100.0) + 0.05]), 2*np.pi / HighNoise2Tbars, marker="x", color="green")
    ax.errorbar(np.array([np.log10(100.0) + 0.05]), 2*np.pi / HighNoise2Tbars,\
            yerr=(2*np.pi*np.sqrt(HighNoise2VarTs)/(HighNoise2Tbars**2)), mec="red")

    ax.set_title(r"$T=" + str(simConfig.Simulation["T"]) + ", N=" + \
            str(simConfig.Simulation["Ensemble Size"]) + ", dt=" + \
            str(simConfig.Simulation["dt"]) + "$")
    ax.set_xlabel(r"$\log{D}$")
    ax.set_ylabel(r"$\overline{\omega}$")
    ax.legend([ "Analytisch", r"$D \rightarrow \infty$", r"$D \rightarrow 0$", \
            "Simulation",r"$dt=0.00001, T=500, N=1000$", r"$dt=0.000001, T=500, N=100$"])
        
    plt.savefig("OmegaBar.pgf")
    plt.savefig("OmegaBar.pdf")
