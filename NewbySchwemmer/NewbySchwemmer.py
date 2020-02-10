from datetime import datetime

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import multiprocessing as mp

t0 = 0.0
T = 10.0
dt = 0.001

N = 1000 # ensemble size

R0 = 1.5
Phi0 = np.pi / 1.2
x0 = np.array([R0, Phi0])

P = { # model parameters
        "D": 0.198,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
        }

def polar2XY(Rho, Phi):

    X = Rho*np.cos(Phi)
    Y = Rho*np.sin(Phi)

    return (X, Y)

def NewbySchwemmer(x, P):

    R = -P["gamma"]*x[0]*(x[0]**2 - 1.0)
    Phi = P["omega"]*(1 + P["gamma"]*P["c"]*(1.0 - x[0])**2)

    return np.array([R, Phi])


def StuartLandau(x, P):

    R = (1.0 - x[0]**2)*x[0]
    Phi = P["omega"]

    return np.array([R, Phi])

def PlanarPolarHeunIteration(x, F, P, dt, dW):
    # x - state vector
    # F - handler of a function with argument x and P
    #     returning vector
    # dt - time step
    # dw - Wiener Process Differential

    # stochastic increment
    dS = np.zeros(2)
    dS[0] = (P["D"] / x[0])*dt + np.sqrt(2*P["D"])*dW[0]
    dS[1] = np.sqrt(2*P["D"])*dW[1] / x[0]
    
    # predictor step
    yp = np.zeros(2)
    dF = dt*F(x, P)
    yp = x + dF + dS

    # corrector step
    y = np.zeros(2)
    dFp = dt*F(yp, P)
    y = x + (dFp + dF) / 2 + dS

    return y


D = pd.DataFrame()

def SolveSDE(i):
   
    print("solve", i)

    Ndof = int((T - t0) / dt)
    x = np.zeros((Ndof , 2))
    x[0] = x0
    dW = np.random.normal(loc=0.0, scale=1.0, size=(Ndof - 1, 2))*np.sqrt(dt)
    
    for t in range(1, Ndof):
        x[t] = PlanarPolarHeunIteration(x[t-1], NewbySchwemmer, P, dt, dW[t-1])

    return (i, x)

def Clbk(result):

    global D

    i = result[0]
    x = result[1]

    D["R" + str(i)] = x[:, 0]
    D["Phi" + str(i)] = x[:, 1]


def Errclbk(R):
    print("error callback")


if __name__ == "__main__":         
    
    pool = mp.Pool(mp.cpu_count())
    for i in range(N):
        pool.apply_async(SolveSDE, args=(i,), callback=Clbk)#, error_callback=Errclbk)
    pool.close()
    pool.join()

    print(D)
