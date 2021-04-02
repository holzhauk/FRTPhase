import numpy as np
from scipy.integrate import trapz, cumtrapz
from .ModelZoo import IsotropicPlanarSSDE


class IsoSurface():
    def __init__(self, model: IsotropicPlanarSSDE, rho: np.array):
        self.model = model
        self.rho = rho
        self.phi = None

    def run_calculations(self):
        pass 

    def get_curve(self) -> (np.array, np.array):
        pass


class ItoIsochron(IsoSurface):
    def __init__(self, model: IsotropicPlanarSSDE, rho: np.array):
        super(ItoIsochron, self).__init__(model, rho)
        self.OmegaBar = None

    def calc_OmegaBar(self):

        rho = self.rho
        i1 = cumtrapz(self.model.g(rho) / (self.model.q_rho(rho) ** 2), rho, initial=0.0)
        i1 = i1[-1] - i1
        E = np.exp(-2 * i1) / (self.model.q_rho(rho) ** 2)
        return trapz(self.model.f(rho) * E, rho) / trapz(E, rho)

    def calc_phi(self):

        rho = self.rho
        if self.OmegaBar == None:
            self.OmegaBar = self.calc_OmegaBar()
        Obar = self.OmegaBar
        phi = np.zeros(rho.shape)

        for i in range(len(rho)):
            r1 = rho[:(i + 1)]
            i1 = np.zeros(r1.shape)

            for j in range(len(r1)):
                r2 = r1[:(j + 1)]
                i2 = cumtrapz(self.model.g(r2) / (self.model.q_rho(r2) ** 2), r2, \
                              initial=0.0)
                i2 = i2[-1] - i2

                i1[j] = trapz((self.model.f(r2) - Obar) * np.exp(-2 * i2) / \
                              (self.model.q_rho(r2) ** 2), r2)

            phi[i] = 2 * trapz(i1, r1)

        return phi

    def run_calculations(self):
        self.OmegaBar = self.calc_OmegaBar()
        self.phi = self.calc_phi()

    def get_OmegaBar(self):
        if self.OmegaBar == None:
            self.OmegaBar = self.calc_OmegaBar()
        return self.OmegaBar

    def get_curve(self) -> (np.array, np.array):
        if self.phi == None:
            self.phi = self.calc_phi()
        return (self.rho, self.phi)


class ItoIsovariant(IsoSurface):

    def __init__(self, model: IsotropicPlanarSSDE, rho: np.array):
        super(ItoIsovariant, self).__init__(model, rho)
        self.isochron = ItoIsochron(model, rho)
        self.DeltaBar = None
        self.OmegaBar = None
        self.phi = None

    def calc_DeltaBar(self):

        rho = self.rho
        if self.OmegaBar == None:
            self.OmegaBar = self.isochron.calc_OmegaBar()
        OB = self.OmegaBar

        iE = cumtrapz(self.model.g(rho) / self.model.q_rho(rho) ** 2, rho, \
                      initial=0.0)
        iE = iE[-1] - iE

        i1 = np.zeros(rho.shape)
        for i in range(len(rho)):
            r1 = rho[:i + 1]
            ie = cumtrapz(self.model.g(r1) / self.model.q_rho(r1) ** 2, r1, \
                          initial=0.0)
            ie = ie[-1] - ie
            i1[i] = 2 * trapz(((self.model.f(r1) / OB) - 1) * np.exp(-2 * ie) / \
                              self.model.q_rho(r1) ** 2, r1)

        i_numerator = (i1 ** 2 + \
                       (self.model.q_phi(rho) / (self.model.q_rho(rho) * OB)) ** 2) * \
                      np.exp(-2 * iE)
        I_numerator = trapz(i_numerator, rho)

        I_denominator = \
            trapz((self.model.f(rho) / self.model.q_rho(rho) ** 2) * np.exp(-2 * iE), rho)

        return 2 * np.pi * I_numerator / I_denominator

    def calc_phi(self):

        rho = self.rho
        if self.OmegaBar == None:
            self.OmegaBar = self.isochron.calc_OmegaBar()
        OB = self.OmegaBar
        if self.DeltaBar == None:
            self.DeltaBar = self.calc_DeltaBar()
        DB = self.DeltaBar

        self.phi = np.zeros(rho.shape)
        for i in range(len(self.phi)):
            r1 = rho[:i + 1]
            i1 = np.zeros(r1.shape)
            for j in range(len(i1)):
                r2 = r1[:j + 1]
                iE = cumtrapz(self.model.g(r2) / (self.model.q_rho(r2) ** 2), r2, \
                              initial=0.0)
                iE = iE[-1] - iE

                i3 = np.zeros(r2.shape)
                for k in range(len(i3)):
                    r3 = r2[:k + 1]
                    ie = cumtrapz(self.model.g(r3) / (self.model.q_rho(r3) ** 2), r3, \
                                  initial=0.0)
                    ie = ie[-1] - ie
                    i3[k] = \
                        2 * trapz((self.model.f(r3) / OB - 1) * np.exp(-2 * ie) \
                                  / self.model.q_rho(r3) ** 2, r3)

                i2 = (self.model.f(r2) - \
                      (2 * np.pi / DB) * ((self.model.q_rho(r2) * i3) ** 2) + \
                      (self.model.q_phi(r2) / OB) ** 2) * np.exp(-2 * iE) / \
                     self.model.q_rho(r2) ** 2
                i1[j] = trapz(i2, r2)
            self.phi[i] = 2 * trapz(i1, r1)

        return self.phi

    def run_calculations(self):
        self.OmegaBar = self.isochron.calc_OmegaBar()
        self.DeltaBar = self.calc_DeltaBar()
        self.phi = self.calc_phi()

    def get_OmegaBar(self):
        if self.OmegaBar == None:
            self.OmegaBar = self.isochron.calc_OmegaBar()
        return self.OmegaBar

    def get_DeltaBar(self):
        if self.DeltaBar == None:
            self.DeltaBar = self.calc_DeltaBar()
        return self.DeltaBar

    def get_curve(self) -> (np.array, np.array):
        if self.phi == None:
            self.phi = self.calc_phi()
        return (self.rho, self.phi)
