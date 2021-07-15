from .ModelZoo import IsotropicPlanarSSDE
from scipy.integrate import trapz, cumtrapz
import numpy as np

class StatPDF():
    def __init__(self, model: IsotropicPlanarSSDE, rho: np.array):
        self.model = model
        self.rho = rho
        self.p = None

    def calc_p(self):
        qFactor = 1 / self.model.q_rho(self.rho)**2
        i_exp = qFactor*self.model.g(self.rho)
        I_exp = 2.0*cumtrapz(i_exp, self.rho, initial=0.0)
        A = 2*np.pi*trapz(qFactor*np.exp(I_exp), self.rho)
        self.p = qFactor*np.exp(I_exp) / A

    def get_curve(self):
        if type(self.p) == type(None):
            self.calc_p()
        return (self.rho, self.p)