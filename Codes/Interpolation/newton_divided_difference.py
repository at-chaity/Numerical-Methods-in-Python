import numpy as np
import sympy as sp
from sympy import symbols,simplify
import math

def dividif_pol(x,y,p):
    u=symbols('u')
    n=len(x)
    h=(x[n-1]-x[0])/(n-1)
    dy=np.zeros([n,n+2])
    for i in range(n):
        dy[i][0]=dy[i][0]+x[i]
        dy[i][1]=dy[i][1]+y[i]
    a=n-1
    j=2
    while j<(n+2):
        i=0
        while i<a:
            dy[i][j]=dy[i+1][j-1]-dy[i][j-1]
            i+=1
            a-=1
        j+=1
    def u_cal(u,n):
        prod=1
        for i in range(1,n):
            prod*=(u-i+1)
        return prod
    for j in range(1,n+2):
        sum=dy[0][j]
    for i in range(1,n):
        sum+=(u_cal(u,i)*dy[0][i])/math.factorial(i-1)
    poly=sp.lambdify(u,sum)
    u=(p-x[0])/h
    return poly(u),simplify(sum)












