import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapz
from scipy.special import erf

from SOscFilePy import *
from SimConfigPy import *

class Model():

    def __init__(self):
        self.Parameters = {}

    def f(self, r):
        pass

    def g(self, r):
        pass

    def gbar(self, r):
        pass

    def get_omegabar(self):
        pass

    def get_PhiIso(self):
        pass

    def get_DetPhiIso(self):
        pass


class NewbySchwemmer(Model):

    def __init__(self, parameters={
                                    "D": 0.198,
                                    "omega": 1.0,
                                    "gamma": 15.0,
                                    "c": -15.0
                                    },\
                        rm=0.25,\
                        rp=1.35,\
                        num=200):
        self.Parameters = parameters
        self.rm = rm
        self.rp = rp
        self.num = num
        self.Rhos = np.linspace(self.rm, \
                                self.rp, \
                                num=self.num)
        
    def f(self, r):
        P = self.Parameters
        return P["omega"]*(1 + P["gamma"]*P["c"]*(1 - r)**2)

    def g(self, r):
        P = self.Parameters
        return -P["gamma"]*r*(r**2 - 1)

    def gbar(self, r, D):
        P = self.Parameters
        return -P["gamma"]*r*(r**2 - 1) / D + 1 / r

    def set_Parameters(self, parameters):
        self.Parameters = parameters

    def get_Rhos(self):
        return self.Rhos

    def get_omegabar(self):
        P = self.Parameters
        R = self.Rhos

        f_r = lambda x: x*(1-x)**2
        a = P["gamma"] / (4*P["D"])
        I1 = trapz(f_r(R)*np.exp(-a*(R**2 - 1)**2), R)
        sqrt_a = np.sqrt(a)

        return P["omega"]*(1 + (P["gamma"]*P["c"]*\
                sqrt_a*4/np.sqrt(np.pi))*\
                   (I1 / (erf(sqrt_a*(self.rp**2 - 1)) \
                         - erf(sqrt_a*(self.rm**2 - 1)))))

    def get_PhiIso(self):
        P = self.Parameters
        R = self.Rhos
        Phi = np.zeros(R.shape)
        
        obar = self.get_omegabar()
        a = P["gamma"]/(4*P["D"])
        
        for i in range(len(R)):
            r1 = R[:(i+1)]
            i1 = np.zeros(r1.shape)
    
            for j in range(len(r1)):
                r2 = r1[:(j+1)]
                i2 = np.zeros(r2.shape)
                    
                i1[j] = trapz((self.f(r2) - obar)*r2*\
                            np.exp(-a*(r2**2 - 1)**2), r2)
                
            Phi[i] = trapz(i1*np.exp(a*(r1**2 - 1)**2) / r1, r1)
            
        return Phi / P["D"]
    
    def get_DetPhiIso(self):
        P = self.Parameters
        R = self.Rhos
        Phi = P["omega"]*P["c"]*np.log(((self.rm + 1)**2 / self.rm) \
                * (R / (R + 1)**2))
        return Phi


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <SimConfigFile>.json")
        os._exit(os.EX_IOERR)
        
    simConfig = SimConfig()
    simConfig.load_from_file(Path(sys.argv[1]))


    ##################################################
    # define Isochrone parameters
    NewbySchwemmerAntirotating_parameters = {
                "D": 0.198,
                "omega": 1.0,
                "gamma": 15.0,
                "c": -15.0
            }

    rm = 0.25
    rp = 1.35
    num = 200

    dtype = np.double
    #################################################

    NewbySchwemmerAntirotating = \
        NewbySchwemmer(NewbySchwemmerAntirotating_parameters, \
                           rm, \
                           rp, \
                           num)

    IsoFile = IsochroneSet(simConfig.Paths["In"], 
                simConfig.ModelName)

    NewbySchwemmerAntirotating_parameters["D"] = 0.0
    NewbySchwemmerAntirotating.set_Parameters(\
            NewbySchwemmerAntirotating_parameters)
    IsoFile.add_Isochrone(NewbySchwemmerAntirotating_parameters, \
            NewbySchwemmerAntirotating.get_Rhos(), \
            NewbySchwemmerAntirotating.get_DetPhiIso())

    Ds = np.array([0.1, 0.2, 0.5, 1.0, 10.0, 100.0])
    for D in Ds:
        print("Write D:", D)
        NewbySchwemmerAntirotating_parameters["D"] = D
        NewbySchwemmerAntirotating.set_Parameters(\
                NewbySchwemmerAntirotating_parameters)
        IsoFile.add_Isochrone(NewbySchwemmerAntirotating_parameters, \
                NewbySchwemmerAntirotating.get_Rhos(), \
                NewbySchwemmerAntirotating.get_PhiIso())
            
    NewbySchwemmerAntirotating_parameters["D"] = 0.5
    NewbySchwemmerAntirotating.set_Parameters(\
        NewbySchwemmerAntirotating_parameters)
    IsoFile.add_Isochrone(NewbySchwemmerAntirotating_parameters, \
        NewbySchwemmerAntirotating.get_Rhos(), \
        np.zeros(NewbySchwemmerAntirotating.get_Rhos().shape, \
                 dtype=np.double))
    

    del IsoFile
    print("Everything written")
    del NewbySchwemmerAntirotating
