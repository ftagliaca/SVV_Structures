import numpy as np
from aileronProperties import Aileron

A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)

def MapAeroLoading(filename):
    return np.genfromtxt(filename, delimiter=",")

AeroLoading = MapAeroLoading("data/aerodynamicloada320.dat")

def find_idx(x, z):
    #input: x, z; location of a point on the aileron surface
    #outputs: index_x, index_z: the indices of the location in the weight matrix (p11)

    def find_z(x):
        for idx_z in range(0, 82):
            if A320.z_i(idx_z + 1, i1 = 1) < z <= A320.z_i(idx_z, i1 = 1):
                return idx_z
        if z > A320.z_i(0, i1 = 1):
            return 0
        else:
            return 81
    def find_x(loc_x):
        for idx_x in range(0, 42):
            if A320.x_i(idx_x, i1 = 1) <= loc_x < A320.x_i(idx_x + 1, i1 = 1):
                return idx_x
        if x < A320.x_i(0, i1 = 1):
            return 0
        else:
            return 41
    return find_x(x), find_z(z)

def LinearInterpolatePos(Q1, Q2, x_0, x_1, x):
    return Q1 + (Q2-Q1)/(x_1 - x_0)*(x-x_0)


def integrate_1d(x, y, x_f):
    if x_f > x[-1]:   #if x_f is outside of the range covered by input, we return the total integral of what we can integrate over
        x_f = x[-1]
    if x_f < x[0]:    #if x_f is lower than the lowest x_value, return 0
        return 0
    if len(x) < 2:
        return 0
    total = 0
    i = 1
    while x[i] < x_f:
        total += (x[i] - x[i-1])*(y[i]+y[i-1])/2
        i += 1
    i -= 1
    if x_f > x[i]:
        total += (x_f - x[i])*(LinearInterpolatePos(y[i], y[i+1], x[i], x[i+1], x_f) + y[i])/2
    return total

def integrate_1d_list(x, y, x_f):
    int_list = []
    x_new = []
    i = 1
    if x_f > x[-1]:
        x_f = x[-1]
    if x_f < x[0]:
        return [0], [x_f]
    if len(x) < 2:
        return [0], [x_f]
    int_list.append(0)
    x_new.append(x[0])
    while x[i] < x_f:
        int_list.append(integrate_1d(x, y, x[i]))
        x_new.append(x[i])
        i += 1
    i -= 1
    if x_f > x[i]:
        int_list.append(integrate_1d(x, y, x_f))
        x_new.append(x_f)

    return int_list, x_new


def integrate_1d_tau(x, y, x_f, x_sc):
    #inputs: x; an array containing all the x-locationx of the points. y; an array containing all the y values of the points. x_i; the value of x to where we integrate.
    #outputs int_x; the total integrated value until x_i.

    #assert len(x) == len(y) #both arrays must have the same size

    if x_f > x[-1]:   #if x_f is outside of the range covered by input, we return the total integral of what we can integrate over
        x_f = x[-1]
    if x_f < x[0]:    #if x_f is lower than the lowest x_value, return 0
        return 0
    if len(x) < 2:
        return 0
    total = 0
    i = 1
    while x[i] < x_f:
        total += (x[i] - x[i-1])*(y[i]+y[i-1])/2 * ((x[i] + x[i-1])*0.5 - x_sc)
        #print(((x[i] + x[i-1])*0.5 - x_sc))
        #print(x[i-1], x[i], total)
        i += 1
    i -= 1
    if x_f > x[i]:
        total += (x_f - x[i])*(LinearInterpolatePos(y[i], y[i+1], x[i], x[i+1], x_f) + y[i])/2* ((x[i] + x_f)*0.5 - x_sc)
        #print((x[i] + x_f)*0.5 - x_sc)
        #print(x[i], x_f, total)
    return total


def integrate_1d_list_tau(x, y, x_f, x_sc):
    #inputs: x, y; lists containing the locations and values of all data points. x_f; the maximum location until which we integrate
    #outputs: x_new; a list containing all the data locatations up to x_f, and x_f if that is not already in the list. int_list; a list containing the integrated values at each location in x_new
    int_list = []
    x_new = []
    i = 1
    if x_f > x[-1]:
        x_f = x[-1]
    if x_f < x[0]:
        return [0], [x_f]
    if len(x) < 2:
        return [0], [x_f]
    int_list.append(0)
    x_new.append(x[0])
    while x[i] < x_f:
        int_list.append(integrate_1d_tau(x, y, x[i], x_sc))
        x_new.append(x[i])
        i += 1
    i -= 1
    if x_f > x[i]:
        int_list.append(integrate_1d_tau(x, y, x_f, x_sc))
        x_new.append(x_f)

    return int_list, x_new




#Making two arrays with all the x, z locationx
def make_x_z():
    x = []
    z = []
    for i in range(0, 41):
        x.append(A320.x_i(i, i1=1))
    for i in range(0, 81):
        z.append(A320.z_i(i, i1=1))
    z = z[::-1] #so it goes from -C to 0 instead of 0 to -C
    return x, z

x, z = make_x_z()
#Define w_bar
def make_w_bar(AeroLoading = AeroLoading):
    x, z = make_x_z()
    w_bar = []
    for i in range(len(x)):
        w_bar = [integrate_1d(z, AeroLoading[:,i], z[-1])] + w_bar
    return w_bar
w_bar = make_w_bar()

def make_tau(x_sc, AeroLoading = AeroLoading):
    x, z = make_x_z()
    tau = []
    for i in range(len(x)):
        tau = [integrate_1d_tau(z, AeroLoading[:,i], z[-1], x_sc)] + tau
    return tau



def DoubleIntegral(x_f):
    x = make_x_z()[0]
    w_bar = make_w_bar()

    return integrate_1d(x, w_bar, x_f)

def ThreeIntegral(x_f):
    x = make_x_z()[0]
    # Integration 1
    w_bar = make_w_bar()

    # Integration 2
    int_list_2, x_list_2 = integrate_1d_list(x, w_bar, x_f)

    # Integration 3
    return integrate_1d(int_list_2, x_list_2, x_f)


def FiveIntegral(x_f):
    x = make_x_z()[0]
    w_bar = make_w_bar()

    int_list_2, x_list_2 = integrate_1d_list(x, w_bar, x_f)

    int_list_3, x_list_3 = integrate_1d_list(x_list_2, int_list_2, x_f)

    int_list_4, x_list_4 = integrate_1d_list(x_list_3, int_list_3, x_f)

    return integrate_1d(x_list_4, int_list_4, x_f)

def DoubleIntegralZSC(x_f, z_sc):
    x = make_x_z()[0]
    tau = make_tau(z_sc)
    return integrate_1d(x, tau, x_f)

def TripleIntegralZSC(x_f, z_sc):
    x = make_x_z()[0]

    tau = make_tau(z_sc)

    int_list_2, x_list_2 = integrate_1d_list(x, tau, x_f)

    return integrate_1d(x_list_2, int_list_2, x_f)
