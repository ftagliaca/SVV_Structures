import numpy as np
from aileronProperties import Aileron
from aero_loads import AerodynamicLoad

A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
AELoading = AerodynamicLoad(A320, "data/aerodynamicloada320.dat")

def MapAeroLoading(filename):
    return np.genfromtxt(filename, delimiter=",")

z,x,AeroLoading = AELoading.interpolate_predefined_grid()

def LinearInterpolatePos(Q1, Q2, x_0, x_1, x):
    #print("Oh come on!")
    #print(x)
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


#Define w_bar
def make_w_bar(AeroLoading = AeroLoading):
    w_bar = []
    for i in range(len(x)):
        w_bar = [integrate_1d(z, AeroLoading[:,i], z[-1])] + w_bar
    return w_bar

#Define tau
def make_tau(x_sc, AeroLoading = AeroLoading):
    tau = []
    for i in range(len(x)):
        tau = [integrate_1d_tau(z, AeroLoading[:,i], z[-1], x_sc)] + tau
    return tau

w_bar = make_w_bar()

def Integral(x_f, p = 2):
    ret_list_2 = w_bar
    x_list_2 = x
    for _ in range(p-2):
        ret_list_2, x_list_2 = integrate_1d_list(x, w_bar, x_f)

    return integrate_1d(x_list_2, ret_list_2, x_f)


def IntegralShear():
    pass

def DoubleIntegralZSC(x_f, z_sc):
    tau = make_tau(z_sc)
    return integrate_1d(x, tau, x_f)

def TripleIntegralZSC(x_f, z_sc):
    tau = make_tau(z_sc)

    int_list_2, x_list_2 = integrate_1d_list(x, tau, x_f)

    return integrate_1d(x_list_2, int_list_2, x_f)
