
import sympy as sp
from sympy import symbols,simplify
from scipy.interpolate import lagrange
def lagrange_pol(x,y,p):
    u=symbols('u')
    n=len(x)
    b=0
    for i in range(n):
        prod = 1
        for k in range(n):
            if i!=k:
                prod*=(u-x[k])/(x[i]-x[k])
        b+=y[i]*prod
    poly=sp.lambdify(u,b)
    ply = lagrange(x, y)
    return poly(p),simplify(b),ply(x)



