import unittest
from Toolbox.SPhaseFile import *

class TestSPhaseFileFunctionality(unittest.TestCase):

    def test_WriteReadIsoSurfaceFile(self):
        isoSurfaceFile = IsoSurfaceFile("theTestModel")
        c1 = isoSurfaceFile.createCurve("Curve0")
        c1.parameterSet["D"] = 0.1256
        c1.parameterSet["gamma"] = 15.0
        c1.rho = np.linspace(0.1, 1.6, num=100)
        c1.phi = np.exp(-c1.parameterSet["gamma"]*c1.rho)
        c1.omegaBar = 0.1324

        c2 = isoSurfaceFile.createCurve("Curve1")
        c2.parameterSet["alpha"] = 2.71
        c2.parameterSet["beta"] = 3.141592
        c2.rho = np.linspace(0.3, 3.67, num=50)
        c2.phi = c2.parameterSet["beta"]*c2.rho
        c2.omegaBar = -5.424
        
        isoSurfaceFile.write("pythonTestFile.h5")

        readFile = IsoSurfaceFile("theTestModel")
        readFile.read("pythonTestFile.h5")
        self.assertEqual(isoSurfaceFile, readFile)
        os.remove("pythonTestFile.h5")

    def test_WriteReadCorrelationFile(self):
        noLags = 10
        N = 10000
        offset = 100
        corrWriteFile = SerialCorrFile("theTestModel")
        corrWriteFile.isoSurfaceFilePath = "theIsoSurfaceFile.h5"
        corrWriteFile.configFilePath = "theConfigFile.h5"
        c1 = corrWriteFile.createIsoSurfaceCorr("Curve0")
        c1.N = np.array([N], dtype=np.uint32)
        c1.offset = np.array([offset], dtype=np.uint32)
        c1.rho_k = np.random.random(size=(noLags,))

        c2 = corrWriteFile.createIsoSurfaceCorr("Curve1")
        c2.N = np.array([N], dtype=np.uint32)
        c2.offset = np.array([offset], dtype=np.uint32)
        c2.rho_k = np.random.random(size=(noLags,))
        corrWriteFile.write("pythonTestFile.h5")

        corrReadFile = SerialCorrFile("theTestModel")
        corrReadFile.read("pythonTestFile.h5")
        self.assertEqual(corrWriteFile, corrReadFile)
        os.remove("pythonTestFile.h5")




if __name__=="__main__":
    unittest.main()