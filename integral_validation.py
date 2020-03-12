from integrals import integrate_1d
import numpy as np
from math import sin, cos

def f2(x):
    return x**2
def fsin(x):
    return np.sin(x)

x = np.linspace(0,100,101)
y1 = f2(x)
y2 = fsin(x)

print((integrate_1d(x, y1, 10)-10**3/3)/(10**3 / 3)*100)
