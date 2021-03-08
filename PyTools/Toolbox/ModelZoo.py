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
        return -P["gamma"] * rho * (rho ** 2 - 1.0) + P["D"] / rho

    def f(self, rho: np.array) -> np.array:
        P = self.pSet
        return P["omega"] * (1.0 + P["gamma"] * P["c"] * (1.0 - rho) ** 2)

    def q_rho(self, rho: np.array) -> np.array:
        P = self.pSet
        return np.ones(rho.shape) * np.sqrt(2 * P["D"])

    def q_phi(self, rho: np.array) -> np.array:
        P = self.pSet
        return np.sqrt(2 * P["D"]) / rho


class SchwabedalPikovsky(IsotropicPlanarSSDE):
    modelName = "SchwabedalPikovsky"
    pSet = {
        "sigma": 0.05,
        "omega": 0.025,
        "delta": 0.01,
        "c": 2.1
    }

    def g(self, rho: np.array) -> np.array:
        P = self.pSet
        return rho * (1.0 - rho) * (3.0 - rho) * (P["c"] - rho) + rho * P["sigma"] ** 2 / 2

    def f(self, rho: np.array) -> np.array:
        P = self.pSet
        return P["omega"] + P["delta"] * (rho - 2.0) - (1.0 - rho) * (3.0 - rho)

    def q_rho(self, rho: np.array) -> np.array:
        P = self.pSet
        return P["sigma"] * rho

    def q_phi(self, rho: np.array) -> np.array:
        return np.zeros(rho.shape)


class ModelFactory():

    def __init__(self):
        self.NewbyS = NewbySchwemmer()
        self.SchwabP = SchwabedalPikovsky()

    def create(self, modelName: str, parameters: dict = {}):
        if (modelName == self.NewbyS.modelName):
            return NewbySchwemmer(parameters)
        if (modelName == self.SchwabP.modelName):
            return SchwabedalPikovsky(parameters)
