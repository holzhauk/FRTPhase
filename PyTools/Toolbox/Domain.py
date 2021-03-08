import numpy as np


class Domain():

    def apply_boundary_conditions(self, x_init: np.array, x_final: np.array) -> np.array:
        pass


class ReflectiveAnnulus(Domain):

    def __init__(self, rho_min: float, rho_max: float):
        if (rho_max < rho_min):
            self.rho_min = rho_max
            self.rho_max = rho_min
        else:
            self.rho_min = rho_min
            self.rho_max = rho_max

        self.width = np.abs(rho_max - rho_min)


    def apply_boundary_conditions(self, x_init: np.array, x_final: np.array) -> np.array:
        if (x_final[0] < self.rho_min):
            Drho = np.abs(x_final[0] - self.rho_min)
            m = np.floor(Drho / self.width)
            if ((m % 2) == 0):
                return np.array([self.rho_min + (Drho - m*self.width), x_final[1]])
            else:
                return np.array([self.rho_max - (Drho - m*self.width), x_final[1]])
        else:
            if (x_final[0] > self.rho_max):
                Drho = np.abs(x_final[0] - self.rho_max)
                m = np.floor(Drho / self.width)
                if ((m % 2) == 0):
                    return np.array([self.rho_max - (Drho - m*self.width), x_final[1]])
                else:
                    return np.array([self.rho_min + (Drho - m*self.width), x_final[1]])
            else:
                return x_final