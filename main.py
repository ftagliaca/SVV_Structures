from aileronProperties import Aileron
from internalLoadsStresses import solveInternal, v, w
from aero_loads import AerodynamicLoad
import numpy as np
from matplotlib import pyplot as plt
import time

t0 = time.time()
A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
# _ = A320.crossArea()
# print(_)
# _ = A320.stringersPosition()
# print(_)
# _ = A320.zCentroid()
# print(_)
# _ = A320.momInertia()
# print(_)

X = np.linspace(0, A320.l_a, 10)
V_p = np.zeros(10)
W_p = np.zeros(10)
for i,x in enumerate(X):
    V_p[i] = v(x)
    W_p[i] = w(x)

V = V_p*np.cos(A320.theta) + W_p*np.sin(A320.theta)
W = W_p*np.cos(A320.theta) + V_p*np.sin(A320.theta)

plt.subplot(121)
plt.plot(X,V)
plt.subplot(122)
plt.plot(X,W)

plt.show()




t1 = time.time()
dt = t1-t0
print("Time taken to execute program ", dt/60, "min")
