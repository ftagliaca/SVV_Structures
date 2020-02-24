from math import sqrt, cos, sin, tan
import numpy as np
from tools import macaulay, integrate2D

def normalStress(y, z, Aileron, M_z, M_x):
    '''
    Input:

    y = y-coordinate
    z = z-coordinate
    Aileron = aileron class containing geometrical and material properties
    M_z = moment around z-axis
    M_x = moment around x-axis

    Output:

    sigma_z = normal stress along z axis, in Pa
    '''
    sigma_z = (M_z*Aileron.Iyy*y + M_x*Aileron.Izz*z)/(Aileron.Izz*Aileron.Iyy)
    return sigma_z

def vonMises(sigma, tau):
    '''
    Input:

    sigma = list with len = 3 containing [sigma_xx, sigma_yy, sigma_zz]
    tau = list with len = 3 containing [tau_xy, tau_xz, tau_yz]

    Output:

    sigma_vm = Von Mises stress
    '''
    sigma_vm = sqrt(((sigma[0]-sigma[1])**2 + (sigma[1]-sigma[2])**2 + (sigma[0]-sigma[2])**2)/2\
    3*(tau[0]**2+tau[1]**2+tau[2]**2))
    return sigma_vm

def solveInternal(alr, q):
    '''
    Input:

    alr = aileron class, containing all the geometrical properties
    q = distributed load q (function)

    Output:

    X = numpy matrix of size (11,1) containing all the unknown as follows
        C1, C2, C3, C4, F_1y, F_1z, F_2y, F_2z, F_3y, F_3z, P_I
    '''

    A = np.matrix([[0,0,alr.x_1,1,0,0,0,0,0,0,0],
                   [0,0,alr.x_2,1,0,((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Iyy),0,0,0,0,cos(alr.theta)*((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Iyy)],
                   [0,0,alr.x_3,1,0,((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Iyy),0,((alr.x_3-alr.x_2)**3)/(6*alr.E*alr.Iyy),0,0,cos(alr.theta)*((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Iyy)],
                   [-alr.x_I,-1,alr.x_I/tan(alr.theta),1/tan(alr.theta),((alr.x_I-alr.x_1)**3)/(6*alr.E*alr.Izz),-((alr.x_I-alr.x_1)**3)/(6*alr.E*alr.Iyy*tan(alr.theta)),0,0,0,0,0],
                   [alr.x_2,1,0,0,-((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Izz),0,0,0,0,0,-((alr.x_2-alr.x_I)**3)/(6*alr.E*alr.Izz)*sin(alr.theta)],
                   [alr.x_1,1,0,0,0,0,0,0,0,0,0],
                   [alr.x_3,1,0,0,-((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Izz),0,-((alr.x_3-alr.x_2)**3)/(6*alr.E*alr.Izz),0,0,0,-((alr.x_3-alr.x_I)**3)*sin(alr.theta)/(6*alr.E*alr.Izz)],
                   [0,0,0,0,0,-1,0,-1,0,-1,-cos(alr.theta)],
                   [0,0,0,0,1,0,1,0,1,0,sin(alr.theta)]])

    d = np.matrix([[-alr.d_1*sin(alr.theta)],
                   [0],
                   [-alr.d_3*sin(alr.theta)+alr.P*cos(alr.theta)*((alr.x_3-alr.x_II)**3)/(2*alr.E*alr.Iyy)],
                   [integrate2D(q,0,alr.x_I,-alr.C_a,0,10,10,p=4)/(alr.E*alr.Izz)],
                   [-integrate2D(q,0,alr.x_2,-alr.C_a,0,10,10,p=4)/(alr.E*alr.Izz)],
                   [alr.d_1*cos(alr.theta)-integrate2D(q,0,alr.x_1,-alr.C_a,0,10,10,p=4)/(alr.E*alr.Izz)],
                   [alr.d_3*cos(alr.theta)-integrate2D(q,0,alr.x_3,-alr.C_a,0,10,10,p=4)/(alr.E*alr.Izz)-alr.P*sin(alr.theta)*((alr.x_3-alr.x_II)**3)/(2*alr.E*alr.Izz)],
                   [-alr.P*cos(alr.theta)],
                   [integrate2D(q,0,alr.x_I,-alr.C_a,0,10,10,p=1)+alr.P*sin(alr.theta)]])
