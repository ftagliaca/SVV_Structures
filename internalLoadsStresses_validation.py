from math import sqrt, cos, sin, tan
import numpy as np
from tools_validation import macaulay, solveInternal
from aileronProperties import Aileron
from aero_loads import AerodynamicLoad
from integrals import FiveIntegral, TripleIntegralZSC, DoubleIntegral, DoubleIntegralZSC, ThreeIntegral

A320 = Aileron(0.605, 2.661, 0.171,1.211,2.591,35,20.5,1.1,2.8,1.2,1.6,1.9,15,1.154,1.840,28,97.4)
Q_coord = np.array([ 0.,0.26610001,0.53220001,0.79829999,  1.06440002,1.3305,1.59659998,1.86269995,2.12880005,2.3948999,2.661])
Q_z = -0.04875
Q = np.array([-0.737, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -0.737])

cF = solveInternal(A320)

def v(x, aileron = A320):
    v  = cF[5]/6*macaulay(x,aileron.x_1)**3
    v += cF[11]/6*macaulay(x,aileron.x_I)**3*sin(aileron.theta)
    v += cF[7]/6*macaulay(x,aileron.x_2)**3
    v += -aileron.P/6*macaulay(x,aileron.x_II)**3*sin(aileron.theta)
    v += cF[9]/6*macaulay(x,aileron.x_3)**3
    for i,q in enumerate(Q):
        v += q*macaulay(x,Q_coord[i])**3
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
    S_z_tot += -DoubleIntegral(x)

    return S_z_tot

def phi(x, aileron = A320):
    z_hat = -0.24023
    J = 0.00022293131689593327
    T = cos(aileron.theta)*aileron.r-sin(aileron.theta)*z_hat

    phi_tot  = cF[11]*macaulay(x, aileron.x_I)*T
    phi_tot += -aileron.P*macaulay(x, aileron.x_II)*T
    phi_tot += -cF[5]*macaulay(x, aileron.x_1)*(z_hat+aileron.r)
    phi_tot += -cF[7]*macaulay(x, aileron.x_2)*(z_hat+aileron.r)
    phi_tot += -cF[9]*macaulay(x, aileron.x_3)*(z_hat+aileron.r)
    phi_tot += -TripleIntegralZSC(x, z_hat)
    phi_tot *= (1/aileron.G*J)
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
    M_z_tot += -ThreeIntegral(x)

    return M_y_tot

def T(x, aileron = A320):
    def dtau(z, x):
        return q(z,x)*(z-z_hat)

    z_hat = -0.215
    T_c = cos(aileron.theta)*aileron.r-sin(aileron.theta)*z_hat
    T_tot  = cF[11]*T_c*macaulay(x, aileron.x_I)**0
    T_tot += -aileron.P*T_c*macaulay(x, aileron.x_II)**0
    T_tot += -DoubleIntegralZSC(x, z_hat)

    return T_tot
