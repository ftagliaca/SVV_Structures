from aileronProperties import Aileron
from internalLoadsStresses import *
from aero_loads import AerodynamicLoad
import numpy as np
from matplotlib import pyplot as plt
import time
from math import cos, sin

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
n = 50
#_,X,_ = AerodynamicLoad(A320, "data/aerodynamicloada320.dat").interpolate_predefined_grid()
X = np.linspace(0, A320.l_a, n)


W_I = np.load("verification_data/defy.npy")
V_I = np.load("verification_data/defx.npy")
P_I = np.load("verification_data/defz.npy")
X_I = np.linspace(0, A320.l_a, len(W_I))
My_I = np.load("verification_data/Myarray.npy")
Mz_I = np.load("verification_data/Mzarray.npy")
Sy_I = np.load("verification_data/Syarray.npy")
Sz_I = np.load("verification_data/Szarray.npy")
T_I  = np.load("verification_data/Tarray.npy")
X_II = np.linspace(0, A320.l_a, len(My_I))


print(np.shape(My_I), np.shape(Mz_I), np.shape(Sy_I), np.shape(Sz_I))


offset = phi(X[0]) - P_I[0]
print(offset/(np.max(P_I)-np.min(P_I))*100)

V = np.zeros(len(X))
W = np.zeros(len(X))
P = np.zeros(len(X))
My = np.zeros(len(X))
Mz = np.zeros(len(X))
Sy = np.zeros(len(X))
Sz = np.zeros(len(X))
Tx  = np.zeros(len(X))

dt_f = 0
for i,x in enumerate(X):

    pd = i/len(X)
    print("Plotting, {0}% done, ETA {1} seconds".format(int(pd*100), round(dt_f*(1-pd)*len(X)),0), end="\r")
    t0_f = time.time()
    P[i]  = phi(x)# - offset
    V[i]  = v(x) + P[i]*(-0.215+A320.r)
    W[i]  = w(x)
    My[i] = M_y(x)
    Mz[i] = M_z(x)
    Sy[i] = S_y(x)
    Sz[i] = S_z(x)
    Tx[i]  = T(x)
    t1_f = time.time()
    dt_f = t1_f-t0_f

P = (P-offset).reshape(len(X),1)

t1 = time.time()
dt = t1-t0
print("Time taken to execute program ", dt/60, "min")

plt.figure("Deflections", figsize =(22,12), dpi = 100)
plt.subplot(221)
plt.plot(X,V)
plt.plot(X_I,V_I)
plt.ylabel("V(x) [m]")
plt.xlabel("x [m]")
plt.title("Deflection in the y direction")
plt.subplot(222)
plt.plot(X,W)
plt.plot(X_I,W_I)
plt.ylabel("W(x) [m]")
plt.xlabel("x [m]")
plt.title("Deflection in the z direction")
plt.subplot(223)
plt.plot(X,P)
plt.plot(X_I,P_I)
plt.ylabel("phi(x) [rad]")
plt.xlabel("x [m]")
plt.title("Rotation around the hinge line")

plt.savefig('figures/Deflections.png')
plt.show()

plt.figure("Moments_Shear", figsize =(22,12), dpi = 100)
plt.subplot(231)
plt.plot(X,-1*My)
plt.plot(X_II,My_I)
plt.ylabel("My(x) [Nm]")
plt.xlabel("x [m]")
plt.title("Moment around the y direction")
plt.subplot(232)
plt.plot(X,Mz)
plt.plot(X_II,Mz_I)
plt.ylabel("Mz(x) [Nm]")
plt.xlabel("x [m]")
plt.title("Moment around the z direction")
plt.subplot(233)
plt.plot(X,Tx)
plt.plot(X_II,T_I)
plt.ylabel("T(x) [Nm]")
plt.xlabel("x [m]")
plt.title("Torque")
plt.subplot(234)
plt.plot(X,Sy)
plt.plot(X_II,Sy_I)
plt.ylabel("Sy(x) [N]")
plt.xlabel("x [m]")
plt.title("Shear force in y direction")
plt.subplot(235)
plt.plot(X,-1*Sz)
plt.plot(X_II,Sz_I)
plt.ylabel("Sz(x) [N]")
plt.xlabel("x [m]")
plt.title("Shear force in z direction")

plt.savefig('figures/Moments_Shear.png')
plt.show()
