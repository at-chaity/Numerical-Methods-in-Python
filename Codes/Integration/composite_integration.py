#This function gives result of integration using Simpson's method, Trapezoidal method and scipy.integrate.quad method

#f: function to be integrated
#a: Lower bound of integral (float)
#b: Upper bound of integral (float)
#n: Number of intervals for simpson's method (must be even integer)
#n1:Number of intervals for Trapezoidal method (integer)
#main_val: Returns the result with scipy.integrate
#error_s: Returns the error of simpson's method
#error_t: Returns the error of Trapezoidal method


import numpy as np
import scipy.integrate as integrate
import sympy as sp
def simpson_trapezoidal(fun_in,a1,b1,n,n1,main_val=False,error_s=False,error_t=False):
    x=sp.symbols('x')
    symbolic_f=sp.sympify(fun_in) #string to symbolic function
    f_num=sp.lambdify(x,symbolic_f,'numpy') #symbolic function to numeric function
    #for pi and exponential value of a and b
    def dec_in(v_str):
        v_str = v_str.strip().lower()
        if v_str in ['pi']:
            return np.pi
        elif v_str in ['e']:
            return np.e
        else:
            return float(v_str)

    a = dec_in(a1)
    b = dec_in(b1)
    h=(b-a)/n
    #simpson composite rule
    xi_0=f_num(a)+f_num(b)
    xi_1=0
    xi_2=0
    i=1
    [integral_1, er] = integrate.quad(f_num, a, b) #solution by scipy
    while i<=n-1:
        x=a+i*h
        if i%2==0:
            xi_2+=f_num(x)
        else:
            xi_1+=f_num(x)
        i+=1
    xi=(h/3)*(xi_0+(2*np.sum(xi_2))+(4*np.sum(xi_1)))
    error_1 = abs(integral_1 - xi)
    #trapezoidal composite rule
    h1=(b-a)/n1
    xt_1=0
    i1=1
    while i1<=n1-1:
        x1=a+i1*h1
        xt_1+=f_num(x1)
        i1+=1
    xt=(h/2)*(xi_0+2*xt_1)
    error_2=abs(integral_1-xt)
    if main_val:
        if error_s and error_t:
            return xi,xt,integral_1,error_1,error_2
        else:
            return xi,xt,integral_1
    else:
        return xi,xt