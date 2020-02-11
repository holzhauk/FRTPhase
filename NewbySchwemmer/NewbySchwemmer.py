from datetime import datetime

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import multiprocessing as mp

def polar2XY(Rho, Phi):

    X = Rho*np.cos(Phi)
    Y = Rho*np.sin(Phi)

    return (X, Y)

def genTstamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

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


def SolveSDE(i, P, x0, dt, t0, T):
   
    #print("solve", i)

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

    if len(sys.argv) != 3:
        print("Usage: python " + sys.argv[0] + " <Ensemble size> <Filename>")
        exit()

    N = int(sys.argv[1])
    PATH = Path(sys.argv[2])

    P = { # model parameters
        "D": 0.198,
        "omega": 1.0,
        "gamma": 15.0,
        "c": -15.0
        }

    t0 = 0.0
    T = 50.0
    dt = 0.01
    
    Ds = np.linspace(0.0, 1.0, num=50)

    for d in Ds:

        FILENAME = str(PATH) + "Esize" + str(N) + "_" + genTstamp()\
                + ".h5"
        store = pd.HDFStore(FILENAME)
        P["D"] = d
        print("D = ", d)
        Pseries = pd.Series(P)
        Pseries = Pseries.append(pd.Series({"t0": t0, "T": T, "dt": dt}))
        store.put("parameters", Pseries)
        D = pd.DataFrame()
        
        Rs = np.linspace(0.01, 1.5, num=20)
        x0 = np.zeros((len(Rs), 2))
        x0[:, 0] = Rs
        x0[:, 1] = np.pi / 2
         
        for j in range(len(Rs)):
            pool = mp.Pool(mp.cpu_count())
            D = pd.DataFrame()
            for i in range(N):
                pool.apply_async(SolveSDE, args=(i, P, x0[j], dt, t0, T),\
                        callback=Clbk)#, error_callback=Errclbk)
            pool.close()
            pool.join()
            store.put("ts" + str(j), D)
            del D


        store.close()
