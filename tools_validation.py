import numpy as np
from math import cos, sin, tan
from aileronProperties import Aileron
from math import cos, sin, sqrt, tan
from integrals import  FiveIntegral, TripleIntegralZSC, DoubleIntegral, DoubleIntegralZSC, ThreeIntegral
import sympy

B737 = Aileron(0.605, 2.661, 0.171,1.211,2.591,35,20.5,1.1,2.8,1.2,1.6,1.9,15,1.154,1.840,28,97.4)
Q_coord = np.array([[ 0.,0.,-0.04875],
 [ 0.26610001,0.,-0.04875   ],
 [ 0.53220001,0.,-0.04875   ],
 [ 0.79829999,  0.,-0.04875   ],
 [ 1.06440002,0.,-0.04875   ],
 [ 1.3305,0.,-0.04875   ],
 [ 1.59659998,0.,-0.04875   ],
 [ 1.86269995,0.,-0.04875   ],
 [ 2.12880005,0.,-0.04875   ],
 [ 2.3948999,0.,-0.04875   ],
 [ 2.661,0.,-0.04875   ]])
Q = [-0.737, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -1.474, -0.737]

def macaulay(x,x1):
    return x-x1 if (x-x1)>0 else 0

def solveInternal(alr: Aileron, q):
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
    G = alr.G
    E = alr.E
    l_a = alr.l_a
    C_a = alr.C_a
    r = alr.r
    d_1 = alr.d_1
    d_3 = alr.d_3
    z_hat = 0.215
    P = alr.P

    T = sin(theta)*z_hat-cos(theta)*r


    '''
    A = np.matrix([[x_1,1,0,0,z_hat+r,0,0,0,0,0,0,0],#v_1+p_1*(z+r) = d_1cos(theta)
                   [0,0,x_1,1,r,0,0,0,0,0,0,0], #w_1-p_1*r = -d_1sin(theta)
                   [x_2,1,0,0,z_hat+r,-(x_2-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_2-x_1)/(G*J),0,0,0,0,0,-sin(theta)*(x_2-x_I)**3/(6*E*I_zz)-T*(x_2-x_I)/(G*J)*(z_hat+r)],#v_2+p_2*(z+r) = 0
                   [0,0,x_2,1,r,-r*(z_hat+r)*(x_2-x_1)/(G*J),(x_2-x_1)**3/6*E*I_yy,0,0,0,0,cos(theta)*(x_2-x_I)**3/(6*E*I_yy)+-r*T*(x_2-x_I)/(G*J)], #w_2-p_2*r = 0
                   [x_3,1,0,0,z_hat+r,-(x_3-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_1)/(G*J),0,-(x_3-x_2)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_2)/(G*J),0,0,0,-sin(theta)*(x_3-x_I)**3/(6*E*I_zz)-T*(x_3-x_I)/(G*J)*(z_hat+r)], #v_3+p_3*(z+r) = d_3cos(theta)
                   [0,0,x_3,1,r,-r*(z_hat+r)*(x_3-x_1)/(G*J),(x_3-x_1)**3/6*E*I_yy,-r*(z_hat+r)*(x_3-x_2)/(G*J),(x_3-x_2)**3/6*E*I_yy,0,0,cos(theta)*(x_3-x_I)**3/(6*E*I_yy)+-r*T*(x_3-x_I)/(G*J)], #w_3-p_3*r = -d_3sin(theta)
                   [x_I*sin(theta),sin(theta),x_I*cos(theta),cos(theta),(z_hat+r)*sin(theta)-r*cos(theta),-sin(theta)*((x_I-x_1)**3/(6*E*I_zz)+(x_I-x_1)*(z_hat+r)**2/(G*J))+cos(theta)*r*(z_hat+r)/(G*J),cos(theta)*(x_I-x_1)**3/(6*E*I_yy),0,0,0,0,0],#sin(theta)*(v_I+p_I*(z+r))+cos(theta)*(w_I-p_I*r) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #S_z(l_a) = 0
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #S_y(l_a) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,sin(theta)*(l_a-x_I)], #M_z(l_a) = 0
                   [0,0,0,0,0,0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(theta)*(l_a-x_I)], #M_x(l_a) = 0
                   [0,0,0,0,1,0,0,0,0,0,0,0]],dtype='float') #phi(0) = 0
                   #[0,0,0,0,0,-(z_hat+r),0,-(z_hat+r),0,-(z_hat+r),0,-T]],dtype='float') #T(l_a) = 0

    b = np.matrix([[d_1*cos(theta)-FiveIntegral(x_1)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [-d_1*sin(theta)-r*TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [0],
                   [-r*TripleIntegralZSC(x_2, z_hat)/(G*J)],
                   [d_3*cos(theta)-FiveIntegral(x_3)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_3,z_hat)/(G*J)-P*(sin(theta)*(x_3-x_II)**3/(6*E*I_zz)+T*(z_hat+r)*(x_3-x_II)/(G*J))],
                   [-d_3*sin(theta)-r*TripleIntegralZSC(x_1,z_hat)/(G*J)+P*(cos(theta)*(x_3-x_II)**3/(6*E*I_yy)-r*T*(x_3-x_II)/(G*J))],
                   [(-(z_hat+r)*sin(theta)+r*cos(theta))*TripleIntegralZSC(x_I,z_hat)/(G*J)-FiveIntegral(x_I)*sin(theta)/(6*E*I_zz)],
                   [P*sin(theta)+DoubleIntegral(l_a)],
                   [-P*cos(theta)],
                   [+P*sin(theta)*(l_a-x_II)+ThreeIntegral(l_a)],
                   [-P*cos(theta)*(l_a-x_II)],
                   [0]],dtype='float')
    '''
    A = np.matrix([[x_1,1,0,0,z_hat+r,0,0,0,0,0,0,0],#v_1+p_1*(z+r) = d_1cos(theta)
                   [0,0,x_1,1,-r,0,0,0,0,0,0,0], #w_1-p_1*r = -d_1sin(theta)
                   [x_2,1,0,0,z_hat+r,-(x_2-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_2-x_1)/(G*J),0,0,0,0,0,-sin(theta)*(x_2-x_I)**3/(6*E*I_zz)-T*(x_2-x_I)/(G*J)*(z_hat+r)],#v_2+p_2*(z+r) = 0
                   [0,0,x_2,1,-r,r*(z_hat+r)*(x_2-x_1)/(G*J),(x_2-x_1)**3/6*E*I_yy,0,0,0,0,cos(theta)*(x_2-x_I)**3/(6*E*I_yy)+r*T*(x_2-x_I)/(G*J)], #w_2-p_2*r = 0
                   [x_3,1,0,0,z_hat+r,-(x_3-x_1)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_1)/(G*J),0,-(x_3-x_2)**3/(6*E*I_zz)-(z_hat+r)**2*(x_3-x_2)/(G*J),0,0,0,-sin(theta)*(x_3-x_I)**3/(6*E*I_zz)-T*(x_3-x_I)/(G*J)*(z_hat+r)], #v_3+p_3*(z+r) = d_3cos(theta)
                   [0,0,x_3,1,-r,r*(z_hat+r)*(x_3-x_1)/(G*J),(x_3-x_1)**3/6*E*I_yy,r*(z_hat+r)*(x_3-x_2)/(G*J),(x_3-x_2)**3/6*E*I_yy,0,0,cos(theta)*(x_3-x_I)**3/(6*E*I_yy)+r*T*(x_3-x_I)/(G*J)], #w_3-p_3*r = -d_3sin(theta)
                   [x_I*sin(theta),sin(theta),x_I*cos(theta),cos(theta),(z_hat+r)*sin(theta)-r*cos(theta),-sin(theta)*((x_I-x_1)**3/(6*E*I_zz)+(x_I-x_1)*(z_hat+r)**2/(G*J))+cos(theta)*r*(z_hat+r)/(G*J),cos(theta)*(x_I-x_1)**3/(6*E*I_yy),0,0,0,0,0],#sin(theta)*(v_I+p_I*(z+r))+cos(theta)*(w_I-p_I*r) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #S_z(l_a) = 0
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #S_y(l_a) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,sin(theta)*(l_a-x_I)], #M_z(l_a) = 0
                   [0,0,0,0,0,0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(theta)*(l_a-x_I)], #M_x(l_a) = 0
                   [0,0,0,0,1,0,0,0,0,0,0,0]],dtype='float') #phi(0) = 0
                   #[0,0,0,0,0,-(z_hat+r),0,-(z_hat+r),0,-(z_hat+r),0,-T]],dtype='float') #T(l_a) = 0


    print(np.linalg.matrix_rank(A))
    _, inds = sympy.Matrix(A).T.rref()
    print(inds)

    b = np.matrix([[d_1*cos(theta)-FiveIntegral(x_1)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [-d_1*sin(theta)+r*TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [-FiveIntegral(x_2)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_2,z_hat)/(G*J)],
                   [r*TripleIntegralZSC(x_2, z_hat)/(G*J)],
                   [d_3*cos(theta)-FiveIntegral(x_3)/(6*E*I_zz)-(z_hat+r)*TripleIntegralZSC(x_3,z_hat)/(G*J)-P*(sin(theta)*(x_3-x_II)**3/(6*E*I_zz)+T*(z_hat+r)*(x_3-x_II)/(G*J))],
                   [-d_3*sin(theta)+r*TripleIntegralZSC(x_1,z_hat)/(G*J)+P*(cos(theta)*(x_3-x_II)**3/(6*E*I_yy)+r*T*(x_3-x_II)/(G*J))],
                   [(-(z_hat+r)*sin(theta)+r*cos(theta))*TripleIntegralZSC(x_I,z_hat)/(G*J)-FiveIntegral(x_I)*sin(theta)/(6*E*I_zz)],
                   [P*sin(theta)+DoubleIntegral(l_a)],
                   [-P*cos(theta)],
                   [+P*sin(theta)*(l_a-x_II)+ThreeIntegral(l_a)],
                   [-P*cos(theta)*(l_a-x_II)],
                   [0]],dtype='float')
                   #[-P*T-DoubleIntegralZSC(l_a,z_hat)]],dtype='float')
    #'''
    '''
    A = np.matrix([[0,0,x_1,1,-r,0,0,0,0,0,0,0], #w(x_1) - p(x_1)*r = -d_1*sin(theta)
                   [0,0,x_2,1,-r,r*(z_hat-r)*(x_2-x_1)/(G*J),((x_2-x_1)**3)/(6*E*I_yy),0,0,0,0,cos(theta)*((x_2-x_1)**3)/(6*E*I_yy)-T*r*(x_2-x_I)/(G*J)], #w(x_2) - p(x_2)*r = 0
                   [0,0,x_3,1,-r,r*(z_hat-r)*(x_3-x_1)/(G*J),((x_3-x_1)**3)/(6*E*I_yy),r*(z_hat-r)*(x_3-x_2)/(G*J),((x_3-x_2)**3)/(6*E*I_yy),0,0,cos(theta)*((x_3-x_1)**3)/(6*E*I_yy)-T*r*(x_3-x_I)/(G*J)], #w(x_3) - p(x_3)*r = -d_3*sin(theta)
                   [x_I*sin(theta),sin(theta),-x_I*cos(theta),-cos(theta),z_hat*sin(theta)-r*cos(theta),-((x_I-x_1)**3)*sin(theta)/(6*E*I_zz)-1/(G*J)*(z_hat-r)*(x_I-x_1)*(z_hat*sin(theta)-r*cos(theta)),((x_I-x_1)**3)*cos(theta)/(6*E*I_yy),0,0,0,0,0],
                   [x_2,1,0,0,z_hat,-((x_2-x_1)**3)/(6*E*I_zz)-z_hat*(z_hat-r)*(x_2-x_1)/(G*J),0,0,0,0,0,-((x_2-x_I)**3)/(6*E*I_zz)*sin(theta)+z_hat*T*(x_2-x_I)/(G*J)], #v(x_2) + p(x_2)*z = 0
                   [x_1,1,0,0,z_hat,0,0,0,0,0,0,0], #v(x_1) + p(x_1)*z = d_1*cos(theta)
                   [x_3,1,0,0,z_hat,-((x_3-x_1)**3)/(6*E*I_zz)-z_hat*(z_hat-r)*(x_3-x_1)/(G*J),0,-((x_3-x_2)**3)/(6*E*I_zz)-z_hat*(z_hat-r)*(x_3-x_2)/(G*J),0,0,0,-((x_3-x_I)**3)*sin(theta)/(6*E*I_zz)+T*z_hat*(x_3-x_I)/(G*J)], #v(x_3) + p(x_3)*z = d_3*cos(theta)
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #S_y(l) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #S_z(l) = 0
                   [0,0,0,0,0,0,x_1-l_a,0,x_2-l_a,0,x_3-l_a,(x_I-l_a)*sin(theta)], #M_y(l) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,(l_a-x_I)*cos(theta)], #M_z(l) = 0
                   [0,0,0,0,0,-(z_hat-r)*(l_a-x_1),0,-(z_hat-r)*(l_a-x_2),0,-(z_hat-r)*(l_a-x_3),0,T]], dtype='float') #T(l) = 0

    b = np.matrix([[-d_1*sin(theta)+TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [+TripleIntegralZSC(x_2,z_hat)/(G*J)],
                   [-d_3*sin(theta)+P*(cos(theta)*((x_3-x_II)**3)/(2*E*I_yy)-T*r*(x_3-x_II))+TripleIntegralZSC(x_3,z_hat)],
                   [-FiveIntegral(x_I)*sin(theta)/(E*I_zz)+TripleIntegralZSC(l_a, z_hat)*(z_hat*sin(theta)-r*cos(theta))],
                   [+FiveIntegral(x_2)/(E*I_zz)-TripleIntegralZSC(x_2,z_hat)/(G*J)],
                   [d_1*cos(theta)+FiveIntegral(x_1)/(E*I_zz)-TripleIntegralZSC(x_1,z_hat)/(G*J)],
                   [d_3*cos(theta)+FiveIntegral(x_3)/(E*I_zz)-P*(sin(theta)*((x_3-x_II)**3)/(2*E*I_zz)-(x_3-x_II)/(G*J))-TripleIntegralZSC(x_3,z_hat)/(G*J)],
                   [-P*cos(theta)],
                   [-ThreeIntegral(l_a)+P*sin(theta)],
                   [P*(x_II-l_a)*cos(theta)],
                   [-DoubleIntegral(l_a,)+P*(l_a-x_II)*sin(theta)],
                   [T*P-DoubleIntegralZSC(l_a, z_hat)]], dtype='float')
    '''
    x = np.linalg.solve(A,b)
    header = 'C1, C2, C3, C4, C5, F_1y, F_1z, F_2y, F_2z, F_3y, F_3z, P_I'
    np.savetxt("reactionForces.dat", x, delimiter=",", header = header)
    print(np.allclose(np.dot(A, x), b))
    return x
