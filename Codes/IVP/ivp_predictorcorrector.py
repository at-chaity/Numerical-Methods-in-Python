#Predictor-corrector method.
#Uses Runge Kutta method for first 3 approximations
#uses Adams-Bashforth method of order 4 for later predictions
#uses Adams-Moulton method of order 3 for corrections

#f_in: The function f(t,y), y'= ")
#a1  : Lower bound of t, a<=t<=b
#b1  : Upper bound of t, a<=t<=b
#w1  : The initial condition, y(a)=y0
#n   : The number of steps
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

def predictor_corrector(f_in,a1,b1,w1,n,return_main=False,return_error=False):
    t, y = sp.symbols('t y')
    f_symbolic = sp.sympify(f_in)
    f = sp.lambdify([t, y], f_symbolic, 'numpy')

    def dec_in(v_str):
        return float(sp.sympify(v_str))

    a = dec_in(a1)
    b = dec_in(b1)
    w = dec_in(w1)
    h=(b-a)/n
    t0=np.linspace(a,b,n+1) #issues with solve_ivp
    w0=np.zeros([n+1])
    t0[0]=a
    w0[0]=w
#Runge-Kutta method for first 3 approximations
    for i in range(1,4):
        k1=h*f(t0[i-1],w0[i-1])
        k2=h*f(t0[i-1]+(h/2),w0[i-1]+(k1/2))
        k3=h*f(t0[i-1]+(h/2),w0[i-1]+(k2/2))
        k4=h*f(t0[i-1]+h,w0[i-1]+k3)
        w0[i]=w0[i-1]+((k1+(2*k2)+(2*k3)+k4)/6)

#predictor-corrector
    for i in range(4,n+1):
        f1=f(t0[i-4],w0[i-4])
        f2=f(t0[i-3],w0[i-3])
        f3=f(t0[i-2],w0[i-2])
        f4=f(t0[i-1],w0[i-1])
        p_w=w0[i-1]+((h/24)*((55*f4)-(59*f3)+(37*f2)-(9*f1))) #predictor 4th order Adam Bashforth
        f5=f(t0[i],p_w)
        w0[i]=w0[i-1]+((h/24)*((9*f5)+(19*f4)-(5*f3)+f2)) #corrector 3rd order Adam Muolton

#scipy exact solution
    sol=solve_ivp(f,[a,b],[w],t_eval=t0) #solves in RK54
    exact_sol=sol.y[0]
    errs=np.abs(exact_sol-w0)
    if return_main:
        if return_error:
            return w0,exact_sol,errs
        else:
            return w0,exact_sol
    else:
        return w0



