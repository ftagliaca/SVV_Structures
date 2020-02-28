from aileronProperties import Aileron
from internalLoadsStresses_validation import solveInternal, v, w, phi
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
n = 100
X = np.linspace(0, A320.l_a, n)
V_p = np.zeros(n)
W_p = np.zeros(n)
P_p = np.zeros(n)
for i,x in enumerate(X):
    V_p[i] = v(x)
    W_p[i] = w(x)
    P_p[i] = phi(x)

V = V_p*np.cos(A320.theta) + W_p*np.sin(A320.theta)
W = W_p*np.cos(A320.theta) + V_p*np.sin(A320.theta)
W = W*-1


plt.subplot(221)
plt.plot(X,V)
plt.ylabel("V(x) [m]")
plt.xlabel("x [m]")
plt.title("Deflection in the y direction")
plt.subplot(222)
plt.plot(X,W)
plt.ylabel("W(x) [m]")
plt.xlabel("x [m]")
plt.title("Deflection in the z direction")
plt.subplot(223)
plt.plot(X,P_p)
plt.ylabel("phi(x) [rad]")
plt.xlabel("x [m]")
plt.title("Rotation around the hinge line")

plt.show()

data = np.vstack((X, V, W, P_p))
np.savetxt("data_validation.dat", data, delimiter=",")




t1 = time.time()
dt = t1-t0
print("Time taken to execute program ", dt/60, "min")
