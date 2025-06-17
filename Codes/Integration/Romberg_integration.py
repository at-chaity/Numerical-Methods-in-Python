# Romberg integration method
#f: function to be integrated
#a: Lower bound of integral (float)
#b: Upper bound of integral (float)
#n: Number of intervals

import numpy as np
import sympy as sp
import scipy.integrate as integrate

def romberg(fun_in,a1,b1,n,main_val=False,return_error=False):
    x=sp.symbols('x')
    #turning string into functions and values
    symbolic_f=sp.sympify(fun_in)
    f_num=sp.lambdify(x,symbolic_f,'numpy')
    def dec_in(v_str):
        v_str=v_str.strip().lower()
        if v_str in ['pi']:
            return np.pi
        elif v_str in ['e']:
            return np.e
        else:
            return float(v_str)
    a=dec_in(a1)
    b=dec_in(b1)
    [integral_1, err] = integrate.quad(f_num, a, b) #result by scipy
    R = np.zeros((2, n))
    h = b - a
    R[0, 0] = (h / 2) * (f_num(a) + f_num(b))
    # Romberg iteration
    for i in range(1, n):
        h /= 2
        sum_f = sum(f_num(a + ( k + 0.5) * h) for k in range(0, 2**(i-1)))
        R[1, 0] = 0.5 * R[0,0] + h * sum_f
        for j in range(1, i+1):
            R[1, j] = R[1, j - 1] + ((R[1, j - 1] - R[0, j - 1]) / (4**j - 1))
        for j in range(0, i+1):
            R[0,j]=R[1,j]
    if main_val:
        if return_error:
            return R[1,n-1],integral_1,abs(integral_1 - R[1, n-1])
        else:
            return R[1,n-1],integral_1
    else:
        return R[1,n-1]




