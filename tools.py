import numpy as np
from math import cos, sin, tan
from aileronProperties import Aileron
from math import cos, sin, sqrt, tan
from integrals import TripleIntegralZSC, DoubleIntegralZSC, Integral
import sympy

def macaulay(x,x1):
    return x-x1 if (x-x1)>0 else 0

def solveInternal(alr: Aileron):
    '''
    Input:

    alr = aileron class, containing all the geometrical properties
    q = distributed load q (function)

    Output:

    x = numpy matrix of size (11,1) containing all the unknown as follows
        C1, C2, C3, C4, C5, F_1y, F_1z, F_2y, F_2z, F_3y, F_3z, P_I
    '''
    theta = alr.theta
    x_1 = alr.x_1
    x_2 = alr.x_2
    x_3 = alr.x_3
    x_I = alr.x_I
    x_II = alr.x_II
    x_a = alr.x_a
    h_a = alr.h
    I_yy = alr.Iyy
    I_zz = alr.Izz
    J = 0.00024311681258111343
    #J = alr.J
    G = alr.G
    E = alr.E
    l_a = alr.l_a
    C_a = alr.C_a
    r = alr.r
    d_1 = alr.d_1
    d_3 = alr.d_3
    z_hat = -0.215
    #z_hat = alr.z_hat
    P = alr.P

    T = -sin(theta)*z_hat-cos(theta)*r

    '''
    A = np.matrix([[x_1,1,0,0,z_hat+r,0,0,0,0,0,0,0],#v_1+p_1*(z+r) = d_1cos(theta)
                   [0,0,x_1,1,-r,0,0,0,0,0,0,0], #w_1-p_1*r = -d_1sin(theta)
                   [x_2,1,0,0,z_hat+r,-(x_2-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_2-x_1)/(G*J),0,0,0,0,0,-sin(theta)*(x_2-x_I)**3/(6*E*I_zz)-T*(x_2-x_I)/(G*J)*(z_hat+r)],#v_2+p_2*(z+r) = 0
                   [0,0,x_2,1,-r,r*(z_hat+r)*(x_2-x_1)/(G*J),(x_2-x_1)**3/(6*E*I_yy),0,0,0,0,cos(theta)*(x_2-x_I)**3/(6*E*I_yy)+r*T*(x_2-x_I)/(G*J)], #w_2-p_2*r = 0
                   [x_3,1,0,0,z_hat+r,-(x_3-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_1)/(G*J),0,-(x_3-x_2)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_2)/(G*J),0,0,0,-sin(theta)*(x_3-x_I)**3/(6*E*I_zz)-T*(x_3-x_I)/(G*J)*(z_hat+r)], #v_3+p_3*(z+r) = d_3cos(theta)
                   [0,0,x_3,1,0,0,(x_3-x_1)**3/(6*E*I_yy),0,(x_3-x_2)**3/(6*E*I_yy),0,0,cos(theta)*(x_3-x_I)**3/(6*E*I_yy)], #w_3 = -d_3sin(theta)
                   #[0,0,x_3,1,-r,r*(z_hat+r)*(x_3-x_1)/(G*J),(x_3-x_1)**3/(6*E*I_yy),r*(z_hat+r)*(x_3-x_2)/(G*J),(x_3-x_2)**3/(6*E*I_yy),0,0,cos(theta)*(x_3-x_I)**3/(6*E*I_yy)+r*T*(x_3-x_I)/(G*J)], #w_3-p_3*r = -d_3sin(theta)
                   [x_I*sin(theta),sin(theta),x_I*cos(theta),cos(theta),(z_hat+r)*sin(theta)-r*cos(theta),-sin(theta)*((x_I-x_1)**3/(6*E*I_zz)+(x_I-x_1)*(z_hat+r)**2/(G*J))+cos(theta)*r*(z_hat+r)/(G*J),cos(theta)*(x_I-x_1)**3/(6*E*I_yy),0,0,0,0,0],#sin(theta)*(v_I+p_I*(z+r))+cos(theta)*(w_I-p_I*r) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #S_z(l_a) = 0
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #S_y(l_a) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,sin(theta)*(l_a-x_I)], #M_z(l_a) = 0
                   [0,0,0,0,0,0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(theta)*(l_a-x_I)], #M_x(l_a) = 0
                   [0,0,0,0,1,0,0,0,0,0,0,0]],dtype='float') #phi(0) = 0.0037654882662205846
                   #[0,0,0,0,0,-(z_hat+r),0,-(z_hat+r),0,-(z_hat+r),0,-T]],dtype='float') #T(l_a) = 0
    '''
    A = np.matrix([[x_1,1,0,0,z_hat+r,0,0,0,0,0,0,0],#v_1+p_1*(z+r) = d_1cos(theta)
                   [0,0,x_1,1,0,0,0,0,0,0,0,0], #w_1 = -d_1sin(theta)
                   [x_2,1,0,0,z_hat+r,-(x_2-x_1)**3/(6*E*I_zz)+(z_hat+r)**2*(x_2-x_1)/(G*J),0,0,0,0,0,-sin(theta)*(x_2-x_I)**3/(6*E*I_zz)+T*(x_2-x_I)/(G*J)*(z_hat+r)],#v_2+p_2*(z+r) = 0
                   [0,0,x_2,1,0,0,-(x_2-x_1)**3/(6*E*I_yy),0,0,0,0,-cos(theta)*(x_2-x_I)**3/(6*E*I_yy)], #w_2 = 0
                   [x_3,1,0,0,z_hat+r,-(x_3-x_1)**3/(6*E*I_zz)+(z_hat+r)**2*(x_3-x_1)/(G*J),0,-(x_3-x_2)**3/(6*E*I_zz)+(z_hat+r)**2*(x_3-x_2)/(G*J),0,0,0,-sin(theta)*(x_3-x_I)**3/(6*E*I_zz)+T*(x_3-x_I)/(G*J)*(z_hat+r)], #v_3+p_3*(z+r) = d_3cos(theta)
                   [0,0,x_3,1,0,0,-(x_3-x_1)**3/(6*E*I_yy),0,-(x_3-x_2)**3/(6*E*I_yy),0,0,-cos(theta)*(x_3-x_I)**3/(6*E*I_yy)], #w_3 = -d_3sin(theta)
                   [x_I*sin(theta)*(z_hat+r),sin(theta)*(z_hat+r),x_I*cos(theta),cos(theta),(z_hat+r),-sin(theta)*(x_I-x_1)**3/(6*E*I_zz)+(x_I-x_1)*(z_hat+r)**2/(G*J),-cos(theta)*(x_I-x_1)**3/(6*E*I_yy),0,0,0,0,0],#sin(theta)*(v_I+p_I*(z+r))+cos(theta)*(w_I) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #R_z(l_a) = 0
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #R_y(l_a) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,sin(theta)*(l_a-x_I)], #M_z(l_a) = 0
                   [0,0,0,0,0,0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(theta)*(l_a-x_I)], #M_y(l_a) = 0
                   [0,0,0,0,0,-(z_hat+r),0,-(z_hat+r),0,-(z_hat+r),0,-T]],dtype='float') #T(l_a) = 0

    #'''
    print(np.linalg.matrix_rank(A))
    _, inds = sympy.Matrix(A).T.rref()
    print(inds)
    '''
    b = np.matrix([[d_1*cos(theta)-FiveIntegral(x_1)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_1,z_hat+r)/(G*J)],
                   [-d_1*sin(theta)+r*TripleIntegralZSC(x_1,z_hat+r)/(G*J)],
                   [-FiveIntegral(x_2)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_2,z_hat+r)/(G*J)],
                   [r*TripleIntegralZSC(x_2, z_hat+r)/(G*J)],
                   [d_3*cos(theta)-FiveIntegral(x_3)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_3,z_hat+r)/(G*J)-P*(sin(theta)*(x_3-x_II)**3/(6*E*I_zz)+T*(z_hat+r)*(x_3-x_II)/(G*J))],
                   [-d_3*sin(theta)+P*cos(theta)*(x_3-x_II)**3/(6*E*I_yy)],
                   #[-d_3*sin(theta)+r*TripleIntegralZSC(x_1,z_hat+r)/(G*J)+P*(cos(theta)*(x_3-x_II)**3/(6*E*I_yy)+r*T*(x_3-x_II)/(G*J))],
                   [(-(z_hat+r)*sin(theta)+r*cos(theta))*TripleIntegralZSC(x_I,z_hat+r)/(G*J)-FiveIntegral(x_I)*sin(theta)/(6*E*I_zz)],
                   [P*sin(theta)+DoubleIntegral(l_a)],
                   [-P*cos(theta)],
                   [+P*sin(theta)*(l_a-x_II)+ThreeIntegral(l_a)],
                   [-P*cos(theta)*(l_a-x_II)],
                   [0.0037654882662205846]],dtype='float')
                   #[-P*T-DoubleIntegralZSC(l_a,z_hat+r)]],dtype='float')
    '''
    b = np.matrix([[d_1*cos(theta)-Integral(x_1,5)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_1,z_hat+r)/(G*J)],
                   [-d_1*sin(theta)],
                   [-Integral(x_2,5)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_2,z_hat+r)/(G*J)],
                   [0],
                   [d_3*cos(theta)-Integral(x_3,5)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_3,z_hat+r)/(G*J)-P*(sin(theta)*(x_3-x_II)**3/(6*E*I_zz)+T*(z_hat+r)*(x_3-x_II)/(G*J))],
                   [-d_3*sin(theta)-P*cos(theta)*(x_3-x_II)**3/(6*E*I_yy)],
                   [-(z_hat+r)/(G*J)*TripleIntegralZSC(x_I,z_hat+r)-Integral(x_I,5)*sin(theta)*(z_hat+r)/(6*E*I_zz)],
                   [P*sin(theta)+Integral(l_a,2)],
                   [-P*cos(theta)],
                   [P*sin(theta)*(l_a-x_II)+Integral(l_a,3)],
                   [P*cos(theta)*(l_a-x_II)],
                   [-P*T-DoubleIntegralZSC(l_a,z_hat+r)]],dtype='float')
    #'''
    x = np.linalg.solve(A,b)
    header = 'C1, C2, C3, C4, C5, F_1y, F_1z, F_2y, F_2z, F_3y, F_3z, P_I'
    np.savetxt("reactionForces.dat", x, delimiter=",", header = header)
    print(np.allclose(np.dot(A, x), b))
    return x
