import numpy as np

def macaulay(x,x1):
    return x-x1 if (x-x1)>0 else 0

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
    def g(x):
        return integrate(lambda y: f(x, y), c, d, n = ny)

    return integrate(g, a, b, n = nx, p = p)
