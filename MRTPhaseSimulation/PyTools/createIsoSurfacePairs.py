import sys, os
from pathlib import Path
import numpy as np
from Toolbox import *

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <SimConfigFile>.json")
        os._exit(os.EX_IOERR)

    config_file_path = Path(sys.argv[1])
    config = SimConfigFile()
    config.read(config_file_path)

    antirotatingPSet = {
        "D": 1.0,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
    }
    NewbySModel = NewbySchwemmer(parameters=antirotatingPSet)
    rho = np.linspace(0.5, 1.2, num=200)
    NewbySIsochron = ItoIsochron(NewbySModel, rho)
    NewbySIsovariant = ItoIsovariant(NewbySModel, rho)

    isoSurfaceFile = IsoSurfaceFile(NewbySModel.get_name())
    isochronCurve = isoSurfaceFile.createCurve("Isochron")
    isochronCurve.parameterSet = NewbySModel.pSet
    (rho, phi) = NewbySIsochron.get_curve()
    isochronCurve.rho = rho.copy()
    isochronCurve.phi = phi.copy()

    isovariantCurve = isoSurfaceFile.createCurve("Isovariant")
    isovariantCurve.parameterSet = NewbySModel.pSet
    (rho, phi) = NewbySIsovariant.get_curve()
    isovariantCurve.rho = rho.copy()
    isovariantCurve.phi = phi.copy()

    isoSurfaceFile.write(config.paths["In"])
