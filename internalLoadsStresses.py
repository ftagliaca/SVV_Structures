from math import sqrt, cos, sin, tan
import numpy as np
from tools import macaulay, integrate2D, solveInternal
from aileronProperties import Aileron
from aero_loads import AerodynamicLoad
from test import FiveIntegral

A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)

aeroLoad = AerodynamicLoad(A320, "data/aerodynamicloada320.dat")
q = aeroLoad.get_values_grid

try:
    cF = np.genfromtxt("reactionForces.dat", delimiter=",", comments = "#")
except OSError as e:
    cF = solveInternal(A320, q)

cF = solveInternal(A320, q)

def normalStress(y, z, Aileron, M_z, M_y):
    '''
    Input:

    y = y-coordinate
    z = z-coordinate
    Aileron = aileron class containing geometrical and material properties
    M_z = moment around z-axis
    M_x = moment around x-axis

    Output:

    sigma_x = normal stress along x axis, in Pa
    '''
    sigma_x = (M_z*Aileron.Iyy*y + M_y*Aileron.Izz*(z-Aileron.zCentroid))/(Aileron.Izz*Aileron.Iyy)
    return sigma_x

def vonMises(sigma, tau):
    '''
    Input:

    sigma = list with len = 3 containing [sigma_xx, sigma_yy, sigma_zz]
    tau = list with len = 3 containing [tau_xy, tau_xz, tau_yz]

    Output:

    sigma_vm = Von Mises stress
    '''
    sigma_vm = sqrt(((sigma[0]-sigma[1])**2 + (sigma[1]-sigma[2])**2 + (sigma[0]-sigma[2])**2)/2+3*(tau[0]**2+tau[1]**2+tau[2]**2))
    return sigma_vm

def v(x, aileron = A320):
    v  = cF[5]/6*macaulay(x,aileron.x_1)**3
    v += cF[11]/6*macaulay(x,aileron.x_I)**3*sin(aileron.theta)
    v += cF[7]/6*macaulay(x,aileron.x_2)**3
    v += -aileron.P/6*macaulay(x,aileron.x_II)**3*sin(aileron.theta)
    v += cF[9]/6*macaulay(x,aileron.x_3)**3
    v += -FiveIntegral(x)
    v *= -1/(aileron.E*aileron.Izz)
    v += cF[0]*x+cF[1]

    return v

def w(x, aileron = A320):
    W  = -cF[6]/6*macaulay(x,aileron.x_1)**3
    W += -cF[11]/6*macaulay(x,aileron.x_I)**3*cos(aileron.theta)
    W += -cF[8]/6*macaulay(x,aileron.x_2)**3
    W += aileron.P/6*macaulay(x,aileron.x_II)**3*cos(aileron.theta)
    W += -cF[10]/6*macaulay(x,aileron.x_3)**3
    W *= -1/(aileron.E*aileron.Iyy)
    W += cF[2]*x+cF[3]

    return W

def S_y(x, aileron = A320):
    S_y_tot  = -cF[6]*macaulay(x,aileron.x_1)**0 if macaulay(x,aileron.x_1)>0 else 0
    S_y_tot += -cF[11]*cos(aileron.theta)*macaulay(x,aileron.x_I)**0 if macaulay(x,aileron.x_I)>0 else 0
    S_y_tot += -cF[8]*macaulay(x,aileron.x_2)**0 if macaulay(x,aileron.x_2)>0 else 0
    S_y_tot += -cF[10]*macaulay(x,aileron.x_3)**0 if macaulay(x,aileron.x_3)>0 else 0
    S_y_tot += aileron.P*cos(aileron.theta)*macaulay(x,aileron.x_II)**0 if macaulay(x,aileron.x_II)>0 else 0

    return S_y_tot

def S_z(x, aileron = A320):
    S_z_tot  = cF[5]*macaulay(x,aileron.x_1)**0 if macaulay(x,aileron.x_1)>0 else 0
    S_z_tot += cF[11]*sin(aileron.theta)*macaulay(x,aileron.x_I)**0 if macaulay(x,aileron.x_I)>0 else 0
    S_z_tot += cF[7]*macaulay(x, aileron.x_2)**0 if macaulay(x,aileron.x_2)>0 else 0
    S_z_tot += cF[9]*macaulay(x, aileron.x_3)**0 if macaulay(x,aileron.x_3)>0 else 0
    S_z_tot += -aileron.P*sin(aileron.theta)*macaulay(x, aileron.x_II)**0 if macaulay(x,aileron.x_II)>0 else 0
    S_z_tot += -integrate2D(q, -aileron.C_a, 0, 0, x, 10, 10, p=1)

    return S_z_tot

def phi(x, aileron = A320):
    z_hat = -0.215
    T = cos(aileron.theta)*aileron.r-sin(aileron.theta)*z_hat

    phi_tot  = cF[11]*macaulay(x, aileron.x_I)*T
    phi_tot += -aileron.P*macaulay(x, x_II)*T
    phi_tot += -integrate2D(q, aileron.C_a, 0, 0, x, 10, 10, p=2)
    phi_tot *= (1/aileron.G*aileron.J)
    phi_tot += cF[4]

    return phi_tot

def M_y(x, aileron = A320):
    M_y_tot  = -cF[6]*macaulay(x, aileron.x_1)
    M_y_tot += -cF[11]*cos(aileron.theta)*macaulay(x, aileron.x_I)
    M_y_tot += -cF[8]*macaulay(x, aileron.x_2)
    M_y_tot += -cF[10]*macaulay(x, aileron.x_3)
    M_y_tot += aileron.P*cos(aileron.theta)*macaulay(x, aileron.x_II)

    return M_y_tot

def M_z(x, aileron = A320):
    M_z_tot  = cF[5]*macaulay(x, aileron.x_1)
    M_z_tot += cF[11]*sin(aileron.theta)*macaulay(x, aileron.x_I)
    M_z_tot += cF[7]*macaulay(x, aileron.x_2)
    M_z_tot += cF[9]*macaulay(x, aileron.x_3)
    M_z_tot += -aileron.P*sin(aileron.theta)*macaulay(x, aileron.x_II)
    M_z_tot += -integrate2D(q, aileron.C_a, 0, 0, x, 0, 0, p=2)

    return M_y_tot

def T(x, aileron = A320):
    def dtau(z, x):
        return q(z,x)*(z-z_hat)

    z_hat = -0.215
    T_c = cos(aileron.theta)*aileron.r-sin(aileron.theta)*z_hat
    T_tot  = cF[11]*T_c*macaulay(x, aileron.x_I)**0
    T_tot += -aileron.P*T_c*macaulay(x, aileron.x_II)**0
    T_tot += -integrate2D(dtau, aileron.C_a, 0, 0, x, 10, 10, p=1)

    return T_tot
