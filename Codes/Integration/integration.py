import numpy as np
import sympy as sp
from sympy import symbols

x=np.array(list(map(float,input("Enter x: ").split())))
y=np.array(list(map(float,input("Enter f(x): ").split())))
n=len(x)
h=(x[n-1]-x[0])/n

#trapezoidal rule
ti=(h/2)*(y[n-1]+y[0])
#simpson's 1/3 rule
si13=(h/3)*(y[n-1]+(4*y[int(np.ceil(n/2))])+y[0])
#simpson's 3/8 rule
si38=(3*h/8)*(y[n-1]+(3*y[int(np.floor(n/2))])+(3*y[int(np.floor(n/2)+1)])+y[0])

print("Result by Trapezoidal rule= ",ti)
print("Result by simpson's 1/3 rule= ",si13)
print("Result by simpson's 3/8 rule= ",si38)

