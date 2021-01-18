import unittest
import numpy as np
from SPhaseFile import *

class TestSPhaseFileFunctionality(unittest.TestCase):

    def test_WriteIsoSurfaceFile(self):
        isoSurfaceFile = IsoSurfaceFile("theTestModel")
        c1 = isoSurfaceFile.createCurve("Curve0")
        c1.parameterSet["D"] = 0.1256
        c1.parameterSet["gamma"] = 15.0
        c1.rho = np.linspace(0.1, 1.6, num=100)
        c1.phi = np.exp(-c1.parameterSet["gamma"]*c1.rho)

        c2 = isoSurfaceFile.createCurve("Curve1")
        c2.parameterSet["alpha"] = 2.71
        c2.parameterSet["beta"] = 3.141592
        c2.rho = np.linspace(0.3, 3.67, num=50)
        c2.phi = c2.parameterSet["beta"]*c2.rho
        
        isoSurfaceFile.write("../../../pythonTestFile.h5")

    def test_ReadIsoSurfaceFile(self):
        isoSurfaceFile = IsoSurfaceFile("theTestModel")




if __name__=="__main__":
    unittest.main()