from StuartLandauSOsc import *
import numpy as np
from matplotlib import pyplot as plt

Trajectories = Ensemble()
Trajectories.readFromFile('', 'simulationCartesian_0.h5')
df = Trajectories.getTimeSeries()

x = df["x1"]
y = df["y1"]

ind = np.ones(len(x))
ind[0] = 0.0

for i in range(0, len(x)-1):
    dx = x[i + 1]-x[i]
    dy = y[i + 1]-y[i]
    r = x[i]**2 + y[i]**2
    ind[i + 1] = ind[i] + (x[i]*dy - y[i]*dx) / (2*np.pi*r)

plt.figure()
plt.plot(df["t"], ind)
plt.show()
