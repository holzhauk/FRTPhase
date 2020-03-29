import matplotlib as mpl
mpl.use('pdf')

import glob

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

from StuartLandauSOsc import Ensemble
from StuartLandauSOsc import XY2polar
from StuartLandauSOsc import polar2XY

def ind(x, y):
    # calculates the instantaneous winding number of a trajectory
    ind = np.ones(len(x))
    ind[0] = 0.0

    for i in range(0, len(x)-1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        r = x[i]**2 + y[i]**2
        ind[i+1] = ind[i] + (x[i]*dy - y[i]*dx) / (2*np.pi*r)

    return ind


def calcCartesianMFPT(path, filename):

    Trajectories = Ensemble()
    Trajectories.readFromFile(path, filename)
    dfCartesian = Trajectories.getTimeSeries()

    (Rho_isochrone, Phi_isochrone) = XY2polar(dfCartesian["x1"].iloc[0],\
            dfCartesian["y1"].iloc[0])

    Ncol = (int)(dfCartesian.shape[1] // 2)
    FP_times = pd.Series([])
    for i in range(1, Ncol + 1):
        ind_i = ind(dfCartesian['x' + str(i)], \
                dfCartesian['y' + str(i)])
        mfpt_index = np.argmax(ind_i > 1.0)

        if mfpt_index != 0:
            FP_times = FP_times.append(pd.Series(\
                    dfCartesian['t'].iloc[mfpt_index]\
                    ), ignore_index=True)

    return (Rho_isochrone, Phi_isochrone, FP_times.mean())
    


def calcPolarMFPT(path, filename):

    Trajectories = Ensemble()
    Trajectories.readFromFile(path, filename)
    dfPolar = Trajectories.getTimeSeries()
    
    Phi_isochrone = dfPolar['Phi1'].iloc[0]
    Rho_isochrone = dfPolar['Rho1'].iloc[0]
    Phi_FP = Phi_isochrone + 2*np.pi
    Ncol = (int)(dfPolar.shape[1] // 2)
    
    FP_times = pd.Series([])
    for i in range(1, Ncol + 1):
        
        dselect = dfPolar.loc[dfPolar['Phi' + str(i)] > Phi_FP]
        if dselect.shape[0] != 0:
            FP_times = FP_times.append(pd.Series([dselect['t'].iloc[0]]),\
                    ignore_index=True)
    
    return (Rho_isochrone, Phi_isochrone, FP_times.mean())


if __name__ == "__main__":

    fig, axs = plt.subplots(2, 2)

    fnames = glob.glob('simulationPolar*')
    Pmfpts = np.zeros(len(fnames))
    Prhos = np.zeros(len(fnames))

    i = 0
    for name in fnames:

        (Prhos[i], phi, Pmfpts[i]) = calcPolarMFPT('', name)
        i += 1

    fnames = glob.glob('simulationCartesian*')
    Cmfpts = np.zeros(len(fnames))
    Crhos = np.zeros(len(fnames))
    
    j = 0
    for name in fnames:
        (Crhos[j], phi, Cmfpts[j]) = calcCartesianMFPT('', name)
        j += 1

    fig.suptitle(r'Mean first passage times; $N = 100$')

    # Cartesian MFPTS subplot
    mnC = np.ones(5)*Cmfpts.mean()
    axs[1,0].plot(np.linspace(0.0, np.amax(Prhos)+0.3, num=5), mnC, 'r')
    axs[1,0].plot(Crhos, Cmfpts, '*')
    axs[1,0].axis([0.0, np.amax(Prhos)+0.3, np.amax(Cmfpts)-3.0, np.amax(Cmfpts)+3.0])
    axs[1,0].set_xlabel(r'$\rho$')
    axs[1,0].set_ylabel(r'$\overline{\tau_{FP}}$')

    # Polar MFPTs subplot
    mnP = np.ones(5)*Pmfpts.mean()
    axs[1,1].plot(np.linspace(0.0, np.amax(Prhos)+0.3, num=5), mnP, 'r')
    axs[1,1].plot(Prhos, Pmfpts, '*')
    axs[1,1].axis([0.0, np.amax(Prhos)+0.3, np.amax(Pmfpts)-3.0, np.amax(Pmfpts)+3.0])
    axs[1,1].set_xlabel(r'$\rho$')
    axs[1,1].set_ylabel(r'$\overline{\tau_{FP}}$')

    Trajectories = Ensemble()
    Trajectories.readFromFile('', 'simulationCartesian_2.h5')
    dfC = Trajectories.getTimeSeries()
    Tx = dfC['x1']
    Ty = dfC['y1']
    ind_T = ind(Tx, Ty)
    I = np.argmax(ind_T > 1.0)
    axs[0,0].plot(Tx[:I+30], Ty[:I+30], 'b')
    axs[0,0].plot(Tx[0], Ty[0], '*r')
    axs[0,0].plot(Tx[I], Ty[I], '*r')
    
    R_iso = np.linspace(0.0, 2.0, num=10)
    Phi_iso = np.ones(len(R_iso))*np.arctan2(Ty[0], Tx[0])
    (X_iso, Y_iso) = polar2XY(R_iso, Phi_iso)
    axs[0,0].plot(X_iso, Y_iso, 'r')

    axs[0,0].axis('equal')
    axs[0,0].set_title(r'Cartesian')
    axs[0,0].set_xlabel(r'x')
    axs[0,0].set_ylabel(r'y')


    Trajectories = Ensemble()
    Trajectories.readFromFile('', 'simulationPolar_2.h5')
    dfP = Trajectories.getTimeSeries()
    TRho = dfP['Rho3']
    TPhi = dfP['Phi3']
    (Tx, Ty) = polar2XY(TRho, TPhi)

    Phi_iso = TPhi.iloc[0] + 2*np.pi
    I = np.argmax(TPhi > Phi_iso)
    axs[0,1].plot(Tx[:I], Ty[:I], 'b')
    axs[0,1].plot(Tx[0], Ty[0], '*r')
    axs[0,1].plot(Tx[I], Ty[I], '*r')
    axs[0,1].plot(X_iso, Y_iso, 'r')
    axs[0,1].axis('equal')
    axs[0,1].set_title(r'Polar')
    axs[0,1].set_xlabel(r'x')
    axs[0,1].set_ylabel(r'y')



    fig.savefig('MFPT_StuartLandauSOsc.pdf')
