import numpy as np
from functools import partial
import time
from math import sin

def integrate(f, a, b, n = 100, p = 1):
    X = np.linspace(a, b, num = n+1)
    h = (b-a)/n
    if p == 1:
        result = (np.sum(f(X))-0.5*f(a)-0.5*f(b))*h
    else:
        #result = (np.sum(integrate(f, 0, X, n = n, p = p-1))-0.5*integrate(f, 0, a, n = n, p = p-1)-0.5*integrate(f, 0, b, n = n, p = p-1))*h
        result = 0
        for i, x in enumerate(X):
            if i != 0:
                result += (integrate(f, 0, X[i], n = n, p = p-1)+integrate(f, 0, X[i-1], n = n, p = p-1))*h/2
    return result


def integrate1(f, a, b, n = 100, p = 1):
    X = np.linspace(a, b, num = n)
    result = 0
    for i, x in enumerate(X):
        if i != 0:
            if p == 1:
                result += (f(X[i])+f(X[i-1]))*(X[i]-X[i-1])/2
            else:
                result += (integrate1(f, 0, X[i], n = n, p = p-1)+integrate1(f, 0, X[i-1], n = n, p = p-1))*(X[i]-X[i-1])/2


    return result

def integrate2D(f, a, b, c, d, nx, ny, p = 1):
    #http://hplgit.github.io/prog4comp/doc/pub/p4c-sphinx-Python/._pylight004.html#reusing-code-for-one-dimensional-integrals
    print("Integrating, please wait...")
    def g(x):
        return integrate(lambda y: f(x, y), c, d, n = ny)

    return integrate(g, a, b, n = nx, p = p)

def integrate2D1(f, a, b, c, d, nx, ny, p = 1):
    #http://hplgit.github.io/prog4comp/doc/pub/p4c-sphinx-Python/._pylight004.html#reusing-code-for-one-dimensional-integrals
    print("Integrating, please wait...")
    def g(x):
        return integrate1(lambda y: f(x, y), c, d, n = ny)

    return integrate1(g, a, b, n = nx, p = p)

f = lambda x: np.sin(x)
g = lambda x,y: x+y
t0 = time.time()
print(integrate(f, 0, 10, n = 100, p = 3))
t1 = time.time()
d_t1 = t1-t0
print("Time took = ", d_t1, "s")
t0 = time.time()
print(integrate1(f, 0, 10, n = 100, p = 3))
t1 = time.time()
d_t2 = t1-t0
print("Time took = ", d_t2, "")
print("The first program was ",d_t2/d_t1, "times faster than the secon one")
print("2D integration")
t0 = time.time()
print(integrate2D(g, 0, 10, 0, 10, 10, 10, p = 3))
t1 = time.time()
d_t1 = t1-t0
print("Time took = ", d_t1, "s")
t0 = time.time()
print(integrate2D1(g, 0, 10, 0, 10, 10, 10, p = 3))
t1 = time.time()
d_t2 = t1-t0
print("Time took = ", d_t2, "")
print("The first program was ",d_t2/d_t1, "times faster than the secon one")
