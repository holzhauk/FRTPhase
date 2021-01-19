import numpy as np

class IsotropicPlanarSSDE():
    modelName = ""
    pSet = {}
    def __init__(self, parameters: dict = {}):
        self.pSet = parameters

    def g(self, rho: np.array) -> np.array:
        pass

    def f(self, rho: np.array) -> np.array:
        pass

    def q_rho(self, rho: np.array) -> np.array:
        pass

    def q_phi(self, rho: np.array) -> np.array:
        pass

    def get_name(self):
        return self.modelName

class NewbySchwemmer(IsotropicPlanarSSDE):
    modelName = "NewbySchwemmer"
    pSet = {
        "D": 0.2,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
    }

    def g(self, rho: np.array) -> np.array:
        P = self.pSet
        return -P["gamma"]*rho*(rho**2 - 1.0) + P["D"] / rho

    def f(self, rho: np.array) -> np.array:
        P = self.pSet
        return P["omega"]*(1 + P["gamma"]*P["c"]*(1 - rho)**2)

    def q_rho(self, rho: np.array) -> np.array:
        P = self.pSet
        return np.ones(rho.shape)*np.sqrt(2*P["D"])

    def q_phi(self, rho: np.array) -> np.array:
        P = self.pSet
        return np.sqrt(2*P["D"]) / rho

