import numpy as np
from aileronProperties import Aileron

def macaulay(x,x1):
    return x-x1 if (x-x1)>0 else 0

def integrate1(f, a, b, n = 100, p = 1):
    X = np.linspace(a, b, num = n)
    h = (b-a)/(n-1)
    if p == 1:
        result = np.sum(f(X))
        result -= f(a)
        result -= f(b)
        result *= h
    else:
        result = 0
        for i, _ in enumerate(X):
            if i != 0:
                result += (integrate1(f, 0, X[i], n = n, p = p-1)+integrate1(f, 0, X[i-1], n = n, p = p-1))*h/2

    return result

def integrate(f, a, b, n = 100, p = 1):
    X = np.linspace(a, b, num = n)
    result = 0
    for i, x in enumerate(X):
        if i != 0:
            if p == 1:
                result += (f(X[i])+f(X[i-1]))*(X[i]-X[i-1])/2
            else:
                result += (integrate(f, 0, X[i], n = n, p = p-1)+integrate(f, 0, X[i-1], n = n, p = p-1))*(X[i]-X[i-1])/2


    return result


def integrate2D(f, a, b, c, d, nx, ny, p = 1):
    #http://hplgit.github.io/prog4comp/doc/pub/p4c-sphinx-Python/._pylight004.html#reusing-code-for-one-dimensional-integrals
    print("Integrating, please wait...")
    def g(x):
        return integrate1(lambda y: f(x, y), c, d, n = ny)

    return integrate1(g, a, b, n = nx, p = p)


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
    z_hat = -0.215
    P = alr.P

    T = cos(theta)*r-sin(theta)*z_hat

    def dtau(z, x):
        return q(z,x)*(z-z_hat)

    A = np.matrix([[0,0,x_1,1,-r,0,0,0,0,0,0,0], #w(x_1) - p(x_1)*r = -d_1*sin(theta)
                   [0,0,x_2,1,-r,0,((x_2-x_1)**3)/(6*E*I_yy),0,0,0,0,cos(theta)*((x_2-x_1)**3)/(6*E*I_yy)-T*r*(x_2-x_I)/(G*J)], #w(x_2) - p(x_2)*r = 0
                   [0,0,x_3,1,-r,0,((x_3-x_1)**3)/(6*E*I_yy),0,((x_3-x_2)**3)/(6*E*I_yy),0,0,cos(theta)*((x_3-x_1)**3)/(6*E*I_yy)-T*r*(x_3-x_I)/(G*J)], #w(x_3) - p(x_3)*r = -d_3*sin(theta)
                   [x_I*sin(theta),sin(theta),-x_I*cos(theta),-cos(theta),z_hat*sin(theta)-r*cos(theta),-((x_I-x_1)**3)*sin(theta)/(6*E*I_zz),((x_I-x_1)**3)*cos(theta)/(6*E*I_yy),0,0,0,0,0],
                   [x_2,1,0,0,z_hat,-((x_2-x_1)**3)/(6*E*I_zz),0,0,0,0,0,-((x_2-x_I)**3)/(6*E*I_zz)*sin(theta)+z_hat*T*(x_2-x_I)/(G*J)], #v(x_2) + p(x_2)*z = 0
                   [x_1,1,0,0,z_hat,0,0,0,0,0,0,0], #v(x_1) + p(x_1)*z = d_1*cos(theta)
                   [x_3,1,0,0,z_hat,-((x_3-x_1)**3)/(6*E*I_zz),0,-((x_3-x_2)**3)/(6*E*I_zz),0,0,0,-((x_3-x_I)**3)*sin(theta)/(6*E*I_zz)+T*z_hat*(x_3-x_I)/(G*J)], #v(x_3) + p(x_3)*z = d_3*cos(theta)
                   [0,0,0,0,0,0,-1,0,-1,0,-1,-cos(theta)], #S_y(l) = 0
                   [0,0,0,0,0,1,0,1,0,1,0,sin(theta)], #S_z(l) = 0
                   [0,0,0,0,0,0,x_1-l_a,0,x_2-l_a,0,x_3-l_a,(x_I-l_a)*sin(theta)], #M_y(l) = 0
                   [0,0,0,0,0,l_a-x_1,0,l_a-x_2,0,l_a-x_3,0,(l_a-x_I)*cos(theta)], #M_z(l) = 0
                   [0,0,0,0,0,0,0,0,0,0,0,T]], dtype='float') #T(l) = 0

    b = np.matrix([[-d_1*sin(theta)-integrate2D(dtau,-C_a,0,0,x_1,10,10,p = 2)/(G*J)],
                   [-integrate2D(dtau,-C_a,0,0,x_2,10,10,p = 2)/(G*J)],
                   [-d_3*sin(theta)+P*(cos(theta)*((x_3-x_II)**3)/(2*E*I_yy)-T*r*(x_3-x_II))-integrate2D(dtau,-C_a,0,0,x_3,10,10,p=2)],
                   [integrate2D(q,-C_a,0,0,x_I,10,10,p=4)*sin(theta)/(E*I_zz)-integrate2D(dtau,-C_a, 0, 0, l_a, 10, 10, p=2)*(z_hat*sin(theta)-r*cos(theta))],
                   [-integrate2D(q,-C_a,0,0,x_2,10,10,p=4)/(E*I_zz)+integrate2D(dtau,-C_a,0,0,x_2,10,10,p=2)/(G*J)],
                   [d_1*cos(theta)-integrate2D(q,-C_a,0,0,x_1,10,10,p=4)/(E*I_zz)+integrate2D(dtau,-C_a,0,0,x_1,10, 10, p=2)/(G*J)],
                   [d_3*cos(theta)-integrate2D(q,-C_a,0,0,x_3,10,10,p=4)/(E*I_zz)-P*(sin(theta)*((x_3-x_II)**3)/(2*E*I_zz)-(x_3-x_II)/(G*J))+integrate2D(dtau,-C_a,0,0,x_3,10, 10, p=2)/(G*J)],
                   [-P*cos(theta)],
                   [integrate2D(q,-C_a,0,0,x_I,10,10,p=1)+P*sin(theta)],
                   [P*(x_II-l_a)*cos(theta)],
                   [integrate2D(q,-C_a,0,0,l_a,10,10,p=2)],
                   [T*P+integrate2D(dtau,-C_a, 0, 0, l_a, 10, 10, p=1)]], dtype='float')
    '''
    A = np.matrix([[0,0,alr.x_1,1,0,0,0,0,0,0,0],
                   [0,0,alr.x_2,1,0,((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Iyy),0,0,0,0,cos(alr.theta)*((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Iyy)],
                   [0,0,alr.x_3,1,0,((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Iyy),0,((alr.x_3-alr.x_2)**3)/(6*alr.E*alr.Iyy),0,0,cos(alr.theta)*((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Iyy)],
                   [-alr.x_I,-1,alr.x_I/tan(alr.theta),1/tan(alr.theta),((alr.x_I-alr.x_1)**3)/(6*alr.E*alr.Izz),-((alr.x_I-alr.x_1)**3)/(6*alr.E*alr.Iyy*tan(alr.theta)),0,0,0,0,0],
                   [alr.x_2,1,0,0,-((alr.x_2-alr.x_1)**3)/(6*alr.E*alr.Izz),0,0,0,0,0,-((alr.x_2-alr.x_I)**3)/(6*alr.E*alr.Izz)*sin(alr.theta)],
                   [alr.x_1,1,0,0,0,0,0,0,0,0,0],
                   [alr.x_3,1,0,0,-((alr.x_3-alr.x_1)**3)/(6*alr.E*alr.Izz),0,-((alr.x_3-alr.x_2)**3)/(6*alr.E*alr.Izz),0,0,0,-((alr.x_3-alr.x_I)**3)*sin(alr.theta)/(6*alr.E*alr.Izz)],
                   [0,0,0,0,0,-1,0,-1,0,-1,-cos(alr.theta)],
                   [0,0,0,0,1,0,1,0,1,0,sin(alr.theta)],
                   [0,0,0,0,0,alr.x_1-alr.l_a,0,alr.x_2-alr.l_a,0,alr.x_3-alr.l_a,(alr.x_I-alr.l_a)*sin(alr.theta)],
                   [0,0,0,0,alr.l_a-alr.x_1,0,alr.l_a-alr.x_2,0,alr.l_a-alr.x_3,0,(alr.l_a-alr.x_I)*cos(alr.theta)]], dtype='float')

    b = np.matrix([[-alr.d_1*sin(alr.theta)],
                   [0],
                   [-alr.d_3*sin(alr.theta)+alr.P*cos(alr.theta)*((alr.x_3-alr.x_II)**3)/(2*alr.E*alr.Iyy)],
                   [integrate2D(q,-alr.C_a,0,0,alr.x_I,10,10,p=4)/(alr.E*alr.Izz)],
                   [-integrate2D(q,-alr.C_a,0,0,alr.x_2,10,10,p=4)/(alr.E*alr.Izz)],
                   [alr.d_1*cos(alr.theta)-integrate2D(q,-alr.C_a,0,0,alr.x_1,10,10,p=4)/(alr.E*alr.Izz)],
                   [alr.d_3*cos(alr.theta)-integrate2D(q,-alr.C_a,0,0,alr.x_3,10,10,p=4)/(alr.E*alr.Izz)-alr.P*sin(alr.theta)*((alr.x_3-alr.x_II)**3)/(2*alr.E*alr.Izz)],
                   [-alr.P*cos(alr.theta)],
                   [integrate2D(q,-alr.C_a,0,0,alr.x_I,10,10,p=1)+alr.P*sin(alr.theta)],
                   [alr.P*(alr.x_II-alr.l_a)*cos(alr.theta)],
                   [integrate2D(q,-alr.C_a,0,0,alr.l_a,10,10,p=2)]], dtype='float')
    '''
    x = np.linalg.solve(A,b)
    header = 'C1, C2, C3, C4, C5, F_1y, F_1z, F_2y, F_2z, F_3y, F_3z, P_I'
    np.savetxt("reactionForces.dat", x, delimiter=",", header = header)
    print(np.allclose(np.dot(A, x), b))
    return x
