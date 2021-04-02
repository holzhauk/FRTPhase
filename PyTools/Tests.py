import unittest
from Toolbox.SPhaseFile import *
from Toolbox.Domain import *
from Toolbox.SDEIntegrators import *


class TestSPhaseFileFunctionality(unittest.TestCase):

    def test_WriteReadIsoSurfaceFile(self):
        isoSurfaceFile = IsoSurfaceFile("theTestModel")
        c1 = isoSurfaceFile.createCurve("Curve0")
        c1.parameterSet["D"] = 0.1256
        c1.parameterSet["gamma"] = 15.0
        c1.rho = np.linspace(0.1, 1.6, num=100)
        c1.phi = np.exp(-c1.parameterSet["gamma"] * c1.rho)
        c1.omegaBar = 0.1324

        c2 = isoSurfaceFile.createCurve("Curve1")
        c2.parameterSet["alpha"] = 2.71
        c2.parameterSet["beta"] = 3.141592
        c2.rho = np.linspace(0.3, 3.67, num=50)
        c2.phi = c2.parameterSet["beta"] * c2.rho
        c2.omegaBar = -5.424

        isoSurfaceFile.write("pythonTestFile.h5")

        readFile = IsoSurfaceFile("theTestModel")
        readFile.read("pythonTestFile.h5")
        self.assertEqual(isoSurfaceFile, readFile)
        os.remove("pythonTestFile.h5")

    def test_WriteReadFRTDataFile(self):
        frtDataFile_write = FRTDataFile("theTestModel")
        dSet1 = frtDataFile_write.createFRTData("Isochron1")
        dSet1.rho0 = np.random.rand(100)
        dSet1.phi0 = np.random.rand(100)
        dSet1.mFRT = np.random.rand(100)
        dSet1.varFRT = np.random.rand(100)
        dSet2 = frtDataFile_write.createFRTData("Isochron100")
        dSet2.rho0 = np.random.rand(5)
        dSet2.phi0 = np.random.rand(5)
        dSet2.mFRT = np.random.rand(5)
        dSet2.varFRT = np.random.rand(5)
        frtDataFile_write.write("FRTTestFile.h5")

        frtDataFile_read = FRTDataFile("theTestModel")
        frtDataFile_read.read("FRTTestFile.h5")
        self.assertEqual(frtDataFile_write, frtDataFile_read)
        os.remove("FRTTestFile.h5")

    def test_WriteReadTbarDataFile(self):
        tbarDataFile_write = TbarDataFile("theTestModel")
        dSet1 = tbarDataFile_write.createTbarData("Isochron0")
        dSet1.Tbar = 2.71
        dSet2 = tbarDataFile_write.createTbarData("Isovariant0")
        dSet2.Tbar = 0.001
        tbarDataFile_write.write("TbarFileTest.h5")

        tbarDataFile_read = TbarDataFile("theTestModel")
        tbarDataFile_read.read("TbarFileTest.h5")
        self.assertEqual(tbarDataFile_write, tbarDataFile_read)
        os.remove("TbarFileTest.h5")

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

    def test_ReadSimConfigFile(self):
        # read test
        simConfigFile = SimConfigFile()
        simConfigFile.read("../ExampleFiles/configs/config_verification.json")
        # test Model content
        self.assertEqual(simConfigFile.modelName, "NewbySchwemmer")
        self.assertEqual(len(simConfigFile.pSets), 1)
        self.assertEqual(simConfigFile.pSets[0]["D"], 0.5)
        self.assertEqual(simConfigFile.pSets[0]["omega"], 1.0)
        self.assertEqual(simConfigFile.pSets[0]["gamma"], 15.0)
        self.assertEqual(simConfigFile.pSets[0]["c"], -15.0)
        # test Domain content
        self.assertEqual(simConfigFile.domainName, "ReflectiveAnnulus")
        self.assertEqual(len(simConfigFile.domainDims), 1)
        self.assertEqual(simConfigFile.domainDims[0]["rho_min"], 0.5)
        self.assertEqual(simConfigFile.domainDims[0]["rho_max"], 1.5)
        # test Simulation content
        simP = simConfigFile.simulation
        self.assertEqual(simP["dt"], 0.001)
        self.assertEqual(simP["t0"], 0.0)
        self.assertEqual(simP["T"], 1000.0)
        self.assertEqual(simP["Ensemble Size"], 1000)
        self.assertEqual(simP["Sample Size"], 5)
        # test Path variable content
        self.assertEqual(simConfigFile.paths["In"], \
            Path("../../../FRTPhase/ExampleFiles/IsoSurfaceFiles" + \
                 "/NewbySchwemmer/AntirotatingD05_IsoSurfacePair.h5"))
        self.assertEqual(simConfigFile.paths["Out"], \
                         Path("../../../FRTPhase/ExampleFiles/SimData" + \
                              "/NewbySchwemmer/AntirotatingD05_FRTFile.h5"))


class TestStochasticIntegration(unittest.TestCase):

    def test_reflectiveAnnulus(self):
        rho_min = 0.5
        rho_max = 1.5
        domain = ReflectiveAnnulus(rho_min, rho_max)
        x0 = np.array([0.6, 0.0])
        x = np.array([0.8, 0.02])
        x_p = domain.apply_boundary_conditions(x0, x)
        self.assertTrue(np.array_equal(x_p, x))
        x = np.array([0.4, 0.02])
        x_p = domain.apply_boundary_conditions(x0, x)
        x_c = np.array([0.6, 0.02])
        self.assertTrue(np.array_equal(x_p, x_c))
        x = np.array([1.7, 0.6])
        x_p = domain.apply_boundary_conditions(x0, x)
        x_c = np.array([1.3, 0.6])
        self.assertTrue(np.array_equal(x_p, x_c))
        x = np.array([3.8, 0.6])
        x_p = domain.apply_boundary_conditions(x0, x)
        self.assertTrue(x_p[0] - 1.2 < 0.00001)
        self.assertEqual(x_p[1], 0.6)
        x = np.array([-3.8, 0.6])
        x_p = domain.apply_boundary_conditions(x0, x)
        self.assertTrue(x_p[0] - 0.8 < 0.00001)
        self.assertEqual(x_p[1], 0.6)


if __name__ == "__main__":
    unittest.main()
