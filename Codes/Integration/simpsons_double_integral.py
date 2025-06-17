#This is for double integral over a non-rectangular region
#works well when the integrand is of two variables
#for loop excludes the end point so it had to increase

##f: function to be integrated
#a: Lower bound of outer integral (float)
#b: Upper bound of outer integral (float)
#c: Lower bound of inner integral (function)
#d: Upper bound of inner integral (function
#m: Number of intervals for outer integral (must be even integer)
#n: Number of intervals for inner integral (must be even integer)
#main_val: Returns the result with scipy.integrate
#error_s: Returns the error of simpson's method

import numpy as np
import sympy as sp
from scipy.integrate import dblquad
def simpson_doublei(fun_in,c1,d1,a1,b1,m,n,main_val=True,return_error=False):
    #string to numerical function
    x=sp.symbols('x')
    y=sp.symbols('y')
    symbolic_f=sp.sympify(fun_in)
    f=sp.lambdify([x,y],symbolic_f,'numpy')
    symbolic_c=sp.sympify(c1)
    c=sp.lambdify(x,symbolic_c,'numpy')
    symbolic_d=sp.sympify(d1)
    d=sp.lambdify(x,symbolic_d,'numpy')
    #a and b in float
    def dec_in(v_str):
        return float(sp.sympify(v_str))
    a=dec_in(a1)
    b=dec_in(b1)
    if c==d and a==b:
        print('Integrated value is 0')
    elif m%2==0 and n%2==0:
        [result, error] = dblquad(lambda y,x: f(x, y), a, b, c, d) #result by scipy
        h = (b - a) / n
        j1 = 0
        j2 = 0
        j3 = 0
        s_in=0
        #simpson's method
        for i in range(0, n + 1):
            x = a + i * h
            hx = (d(x) - c(x)) / m
            k1 = f(x, c(x)) + f(x, d(x))
            k2 = 0
            k3 = 0
            for j in range(1, m):
                y = c(x) + j * hx
                q = f(x, y)
                if j % 2 == 0:
                    k2 += q
                else:
                    k3 += q
            L = (hx / 3) * (k1 + 2 * k2 + 4 * k3) #inner integration
            if i == 0 or i == n:
                j1 += L
            elif i % 2 == 0:
                j2 += L
            else:
                j3 += L
        s_in+=(h/3)*(j1+2*j2+4*j3) #outer integration
        if main_val:
            if return_error:
                return s_in,result,abs(result-s_in)
            else:
                return s_in,result
        else:
            return s_in
    else:
        print("m and n must be even")


