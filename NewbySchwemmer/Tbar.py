import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate

P = {
        "D": 0.198,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
        }

def Tbarnum(P, D, rhom, rhop):

    P["D"] = D

    gbar = lambda r: -(P["gamma"] / P["D"])*(r*(r**2 - 1)) + 1 / r 
    f = lambda r: P["omega"]*(1 + P["gamma"]*P["c"]*(1 - r)**2)

    integrand1 = lambda x: np.exp(-integrate.quad(gbar, x, rhop)[0])
    integrand2 = lambda x: f(x)*integrand1(x)

    Nominator = 2*np.pi*integrate.quad(integrand1, rhom, rhop)[0]
    Denominator = integrate.quad(integrand2, rhom, rhop)[0]

    return Nominator / Denominator

def Phidotbarnum(P, D, rhom, rhop):

    P["D"] = D

    gbar = lambda r: -(P["gamma"] / P["D"])*(r*(r**2 - 1)) + 1 / r 
    f = lambda r: P["omega"]*(1 + P["gamma"]*P["c"]*(1 - r)**2)

    integrand1 = lambda x: np.exp(-integrate.quad(gbar, x, rhop)[0])
    integrand2 = lambda x: f(x)*integrand1(x)

    Denominator = integrate.quad(integrand1, rhom, rhop)[0]
    Nominator = integrate.quad(integrand2, rhom, rhop)[0]

    return Nominator / Denominator

if __name__ == "__main__":

    Ds = np.linspace(0.01, 1.0, num=100)
    PhiDots = np.zeros(Ds.shape)
    for i in range(len(Ds)):
        PhiDots[i] = Phidotbarnum(P, Ds[i], 0.01, 1.5)

    plt.figure()
    plt.plot(np.log(Ds), PhiDots)
    plt.show()
