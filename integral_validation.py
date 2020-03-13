from integrals import integrate_1d, integrate_1d_list
import numpy as np
from math import sin, cos, pi

def f1(x):
    return x**2
def f2(x):
    return x**3-2*x**2+5*x-90
def f3(x,y):
    return x**2 + y**2
def f4(x,y):
    return np.sin(x) + np.cos(y)

x_1 = np.linspace(0,10,101)
y_1 = f1(x_1)

x_2 = np.linspace(0,20,70)
y_2 = f2(x_2)

x_3 = np.linspace(0,10,70)
z_3 = np.linspace(0,5,70)
y_3 = []
for i in range(len(x_3)):
    y_3.append(f3(x_3[i],z_3))
y_3 = np.column_stack(y_3)

x_4 = np.linspace(0,pi,70)
z_4 = np.linspace(0,pi/2,70)
y_4 = []
for i in range(len(x_4)):
    y_4.append(f4(x_4[i],z_4))
y_4 = np.column_stack(y_4)

x_s = np.linspace(0,pi,101)
y_s = np.sin(x_s)

#Compute and print integration error for different functions
ef_1 = - 2*(10-0)**3/(12*100**2)
e_1  = integrate_1d(x_1, y_1, 10)-1000/3
ep_1 = abs(e_1/(5000 / 3))

ef_2 = - (6*20-4)*(20-0)**3/(12*70**2)
e_2  = integrate_1d(x_2, y_2, 20)-101600/3
ep_2 = abs(e_2/(101600/3))

ef_3 = - 4*(10-0)**3*(5-0)**3/(12*70**2*70**2)
w_3 = []
for i in range(len(x_3)):
    w_3.append(integrate_1d(z_3, y_3[:,i], z_3[-1]))
e_3  = integrate_1d(x_3, w_3, x_3[-1]) - 6250/3
ep_3 = abs(e_3/(6250/3))

w_4 = []
for i in range(len(x_4)):
    w_4.append(integrate_1d(z_4, y_4[:,i], z_4[-1]))
e_4  = integrate_1d(x_4, w_4, x_4[-1]) - 2*pi
ep_4 = abs(e_4/(2*pi))

ef_s = - (pi-0)**3/(12*101**2)
e_s  = integrate_1d(x_s, y_s, pi)-2
ep_s = abs(e_s/2)

print("f(x) = x^2, a = 0, b = 10, Dx = {4} \nerror: = {0} or {1}% \ncomputed error = {2} ({3})".format(e_1,round(ep_1*100,4),ef_1,abs(ef_1)>abs(e_1), x_1[1]-x_1[0]))
print("f(x) = x^3-2*x^2+5*x-90, a = 0, b = 20, Dx = {4} \nerror: = {0} or {1}% \ncomputed error = {2} ({3})".format(e_2,round(ep_2*100,4),ef_2,abs(ef_2)>abs(e_2), x_2[1]-x_2[0]))
print("f(x) = sin(x), a = 0, b = pi, Dx = {4} \nerror: = {0} or {1}% \ncomputed error = {2} ({3})".format(e_s,round(ep_s*100,4),ef_s,abs(ef_s)>abs(e_s), x_s[1]-x_s[0]))
print("f(x) = x^2 + y^2, a = 0, b = 10,c = 0, d = 5, Dx = {4}, Dy = {5} \nerror: = {0} or {1}% \ncomputed error = {2} ({3})".format(e_3,round(ep_3*100,4),ef_3,abs(ef_3)>abs(e_3), x_3[1]-x_3[0],z_3[1]-z_3[0]))
print("f(x) = sin(x) + cos(y), a = 0, b = pi,c = 0, d = pi/2, Dx = {4}, Dy = {5} \nerror: = {0} or {1}% \ncomputed error = {2} ({3})".format(e_4,round(ep_4*100,4),1,abs(1)>abs(e_4), x_4[1]-x_4[0],z_4[1]-z_4[0]))
