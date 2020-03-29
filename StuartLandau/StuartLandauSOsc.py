from datetime import datetime

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

class Ensemble:

    def __init__(self):

        self.config = pd.Series() 
        self.timeSeries = pd.DataFrame()

    def setConfig(self, config):

        keys = [k for k in config]
        values = np.array([v for v in config.values()])

        self.config = pd.Series(values, index=keys)

    def getConfig(self):
        return self.config.to_dict()

    def add(self, data, indices):

        if data.ndim == 1:
            self.timeSeries[indices[0]] = data
        else:
            (rows, cols) = data.shape
            for c in range(cols):
                self.timeSeries[indices[c]] = data[:, c]

    def getTimeSeries(self):
        return self.timeSeries

    def saveToFile(self, path, tag):

        store = pd.HDFStore(path + tag + '.h5')
        store['config'] = self.config
        store['timeSeries'] = self.timeSeries
        store.close()

    def readFromFile(self, path, name):

        store = pd.HDFStore(path + name)
        self.config = store['config']
        self.timeSeries = store['timeSeries']
        store.close()


class Model:

    def __init__(self, parameters):

        self.parameters = parameters
        self.Trajectories = Ensemble()

    def writeResults(self, path, tag):
        self.Trajectories.setConfig(self.parameters)
        self.Trajectories.saveToFile(path, tag)


def FCartesian(parameters, x, y):

    p = parameters

    fx = (p['a']- (x**2 + y**2))*x\
            - p['omega']*y
    fy = (p['a'] - (x**2 + y**2))*y\
            + p['omega']*x

    return (fx, fy)

class CartesianModel(Model):
    
    def solve(self, t0, dt, T, x0):

        Ndof = (int)(( T - t0 ) // dt) # integer division

        t = np.linspace(t0, T, Ndof)

        self.Trajectories.add(t, ['t'])

        if x0.ndim == 1:
            NInst = 1 # number of instances
            dim = x0.shape[0]

        else:
            (NInst, dim) = x0.shape

        for j in range(NInst):
         
            x = np.zeros((Ndof, dim))
            dWx = np.random.normal(loc=0.0, scale=1.0, size=Ndof)*np.sqrt(dt)
            dWy = np.random.normal(loc=0.0, scale=1.0, size=Ndof)*np.sqrt(dt)
            
            x[0, :] = x0[j, :]
            g = np.sqrt(2*self.parameters['D'])
            
            for i in range(1, Ndof):
            
                # euler predictor substep
                (fx, fy) = FCartesian(self.parameters, x[i-1, 0], x[i-1, 1])
                x_pred = x[i-1, 0] + dt*fx + g*dWx[i-1]
                y_pred = x[i-1, 1] + dt*fy + g*dWy[i-1]
                
                # corrector substep
                (f_predx, f_predy) = FCartesian(self.parameters, x_pred, y_pred)
                x[i, 0] = x[i-1, 0] + dt*(f_predx + fx)/2\
                         + g*dWx[i-1]
                x[i, 1] = x[i-1, 1] + dt*(f_predy + fy)/2\
                         + g*dWy[i-1]
                
            self.Trajectories.add(x, ['x'+str(j + 1), 'y'+str(j + 1)])


def FPolar(parameters, rho, phi):

    p = parameters

    frho = (p['a'] - rho**2)*rho + p['D']/rho
    fphi = p['omega']

    return (frho, fphi)


class PolarModel(Model):

    def solve(self, t0, dt, T, x0):

        Ndof = (int)(( T - t0 ) // dt) # integer division

        t = np.linspace(t0, T, Ndof)

        self.Trajectories.add(t, ['t'])

        if x0.ndim == 1:
            NInst = 1 # number of instances
            dim = x0.shape[0]

        else:
            (NInst, dim) = x0.shape

        for j in range(NInst):
         
            x = np.zeros((Ndof, dim))
            dWx = np.random.normal(loc=0.0, scale=1.0, size=Ndof)*np.sqrt(dt)
            dWy = np.random.normal(loc=0.0, scale=1.0, size=Ndof)*np.sqrt(dt)
            
            x[0, :] = x0[j, :]
            g = np.sqrt(2*self.parameters['D'])
            
            for i in range(1, Ndof):
            
                # euler predictor substep
                (fRho, fPhi) = FPolar(self.parameters, x[i-1, 0], x[i-1, 1])
                rho_pred = x[i-1, 0] + dt*fRho + g*dWx[i-1]
                phi_pred = x[i-1, 1] + dt*fPhi + g*dWy[i-1] / x[i-1, 1]
                
                # corrector substep
                (f_predRho, f_predPhi) = FPolar(self.parameters, rho_pred, phi_pred)
                x[i, 0] = x[i-1, 0] + dt*(f_predRho + fRho)/2\
                         + g*dWx[i-1]
                x[i, 1] = x[i-1, 1] + dt*(f_predPhi + fPhi)/2\
                         + g*dWy[i-1] / x[i-1, 1]
                
            self.Trajectories.add(x, ['Rho'+str(j + 1), 'Phi'+str(j + 1)])


def genTstamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def polar2XY(Rho, Phi):

    X = Rho*np.cos(Phi)
    Y = Rho*np.sin(Phi)

    return (X, Y)

def XY2polar(X, Y):

    Phi = np.arctan2(Y, X)
    Rho = np.sqrt(X**2 + Y**2)

    return (Rho, Phi)



if __name__ == "__main__":
    
    parameters = {
            'a': 1.0,
            'omega': 0.5,
            'D': 0.01 
            }

    N = 100 # Ensemble size

    t0 = 0.0
    T = 30.0
    dt = 0.01

    Phi_isochrone = -parameters['omega']
    Rhos_isochrone = np.linspace(0.1, 3.0, num=5)

    for i in range(len(Rhos_isochrone)):

        xp0 = np.zeros((N, 2))
        xp0[:, 0] = Rhos_isochrone[i]
        xp0[:, 1] = Phi_isochrone
    
        osp = PolarModel(parameters)
        osp.solve(t0, dt, T, xp0)
        osp.writeResults('', 'simulationPolar_' + str(i))

        (X_isochrone, Y_isochrone) = polar2XY(Rhos_isochrone[i], Phi_isochrone)
        x0 = np.zeros((N, 2))
        x0[:, 0] = X_isochrone
        x0[:, 1] = Y_isochrone

        osc = CartesianModel(parameters)
        osc.solve(t0, dt, T, x0)
        osc.writeResults('', 'simulationCartesian_' + str(i))

   # x0 = np.array([[1.0, 3.0],
   #                 [0.01, 0.03]])
   # 
   # (Rho0, Phi0) = XY2polar(x0[:, 0], x0[:, 1])
   # xp0 = np.array([Rho0, Phi0]).T

   # tstamp = genTstamp()

   # osc = CartesianModel(parameters)
   # print("solve cartesian SDE...")
   # osc.solve(t0, dt, T, x0)
   # osc.writeResults('', 'simulationCartesian_' + tstamp)

   # osp = PolarModel(parameters)
   # print("solve polar SDE...")
   # osp.solve(t0, dt, T, xp0)
   # osp.writeResults('', 'simulationPolar_' + tstamp)

   # Trajectories = Ensemble()
   # Trajectories.readFromFile('', 'simulationCartesian_' + tstamp + '.h5')
   # dfCartesian = Trajectories.getTimeSeries()
   # Trajectories = Ensemble()
   # Trajectories.readFromFile('', 'simulationPolar_' + tstamp + '.h5')
   # dfPolar = Trajectories.getTimeSeries()
   # 
   # print("print trajectories")
   # fig, axs = plt.subplots(1, 2)
   # 
   # axs[0].plot(dfCartesian['x1'], dfCartesian['y1'])
   # axs[0].plot(dfCartesian['x2'], dfCartesian['y2'])
   # axs[0].axis('equal')
   # 

   # (X1, Y1) = polar2XY(dfPolar['Rho1'], dfPolar['Phi1'])
   # axs[1].plot(X1, Y1)
   # (X2, Y2) = polar2XY(dfPolar['Rho2'], dfPolar['Phi2'])
   # axs[1].plot(X2, Y2)
   # axs[1].axis('equal')
   # #plt.show()

   # # test for FP
   # phi_fp = xp0[1,1] + 2*np.pi
   # dselect = dfPolar.loc[dfPolar['Phi2'] > phi_fp]
   # print(dselect['t'].iloc[0])
