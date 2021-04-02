import numpy as np
from .ItoIsoSurfaces import *
from .Exceptions import *


class FirstPassageDetector():

    def __init__(self, isoSurface: IsoSurface):
        self.isoSurface = isoSurface
        if (self.isoSurface.phi == None):
            self.isoSurface.run_calculations()

    def is_first_return_event(self, x: np.array, pos_sense_of_rotation: bool) -> bool:

        rho = self.isoSurface.rho
        phi = self.isoSurface.phi
        dof = rho.size
        rho_min = np.amin(rho)
        rho_max = np.amax(rho)

        if ((x[0] > rho_max) or (x[0] < rho_min)):
            raise OutOfDomain("FirstPassageDetector: the realization lies out of the domain")

        k = int(np.floor((dof - 1.0) * (x[0] - rho_min) / (rho_max - rho_min)))

        if (pos_sense_of_rotation):
            return (x[1] >= phi[k] + 2 * np.pi)
        else:
            return (x[1] <= phi[k] - 2 * np.pi)