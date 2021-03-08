import numpy as np
from .ModelZoo import IsotropicPlanarSSDE
from .Domain import Domain

class SDEIntegrator():

    class Config():
        x0 = np.array([0.0, 0.0])
        t0 = 0.0
        T = 100.0
        dt = 0.001

    def __init__(self, config: Config):
        self.config = config
        self.x = config.x0
        self.t = config.t0
        self.T = config.T
        self.dt = config.dt

    def evolve(self) -> tuple:
        pass

    def configure(self, config: Config):
        self.config = config
        self.x = config.x0
        self.t = config.t0
        self.T = config.T
        self.dt = config.dt

    def reset(self, x0: np.array, t0: float):
        self.x = x0
        self.t = t0

    def is_in_time(self) -> bool:
        return self.t < self.T

    def integrate(self) -> tuple:
        while (self.t < self.T):
            (self.x, self.t) = self.evolve()
        return (self.x.copy(), self.t)


class ItoEulerIntegrator(SDEIntegrator):

    def __init__(self, domain: Domain, model: IsotropicPlanarSSDE):
        super(ItoEulerIntegrator, self).__init__()
        self.domain = domain
        self.model = model

    def evolve(self) -> tuple:
        x0 = self.x.copy()
        dW = np.random.normal(0.0, 1.0, size=2)*np.sqrt(self.dt)
        self.x[0] = x0[0] + self.model.g(x0[0])*self.dt + self.model.q_rho(x0[0])*dW[0]
        self.x[1] = x0[1] + self.model.f(x0[0])*self.dt + self.model.q_phi(x0[0])*dW[1]

        self.x = self.domain.apply_boundary_conditions(x0, self.x)
        self.t += self.dt
        return (self.x.copy(), self.t)


