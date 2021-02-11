import sys
from Toolbox.SPhaseFile import *
from Toolbox.ModelZoo import *
from Toolbox.ItoIsoSurfaces import *


import numpy as np

rho_min = 0.5
rho_max = 1.2
doF = 400

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <SimConfigFile>.json")
        os._exit(os.EX_IOERR)

    config_file_path = Path(sys.argv[1])
    config = SimConfigFile()
    config.read(config_file_path)

    theModelFactory = ModelFactory()
    isoSurfaceFile = IsoSurfaceFile(config.modelName)
    index = 0

    for pSet in config.pSets:

        theModel = theModelFactory.create(config.modelName, parameters=pSet)
        rho = np.linspace(rho_min, rho_max, num=doF)
        NewbySIsochron = ItoIsochron(theModel, rho)
        NewbySIsovariant = ItoIsovariant(theModel, rho)

        isochronCurve = isoSurfaceFile.createCurve("Isochron" + str(index))
        isochronCurve.parameterSet = theModel.pSet
        isochronCurve.omegaBar = NewbySIsochron.get_OmegaBar()
        (rho, phi) = NewbySIsochron.get_curve()
        isochronCurve.rho = rho.copy()
        isochronCurve.phi = phi.copy()

        isovariantCurve = isoSurfaceFile.createCurve("Isovariant" + str(index))
        isovariantCurve.parameterSet = theModel.pSet
        isovariantCurve.omegaBar = NewbySIsovariant.get_OmegaBar()
        (rho, phi) = NewbySIsovariant.get_curve()
        isovariantCurve.rho = rho.copy()
        isovariantCurve.phi = phi.copy()

        index += 1

    isoSurfaceFile.write(config.paths["In"])
