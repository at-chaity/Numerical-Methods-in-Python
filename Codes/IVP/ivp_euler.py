#eulers method gives accurate results when h is small
#sympy doesn't produce good results because the exact solution is complex or the f rhs is implicit

#f_in: The function f(t,y), y'= ")
#a1  : Lower bound of t, a<=t<=b
#b1  : Upper bound of t, a<=t<=b
#w1  : The initial condition, y(a)=y0
#n   : The number of steps
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

def euler(f_in,a1,b1,w1,n,return_main=False,return_error=False):
    t,y=sp.symbols('t y')
    f_symbolic=sp.sympify(f_in)
    f=sp.lambdify([t,y],f_symbolic,'numpy')
    def dec_in(v_str):
        return float(sp.sympify(v_str))
    a=dec_in(a1)
    b=dec_in(b1)
    w=dec_in(w1)
    h=(b-a)/n
    t0=np.linspace(a,b,n+1) #issues with solve_ivp
    w0=np.zeros([n+1])
    t0[0]=a
    w0[0]=w
    # Euler's method
    for i in range(1,n+1):
        w0[i]=w0[i-1]+h*f(t0[i-1],w0[i-1])

#exact solution
#ode=sp.Eq(sp.Derivative(y,t),f_symbolic)
#gen_sol=sp.dsolve(ode,y)
#exact_rhs=gen_sol.rhs #taking the expression of rhs only
#C1=sp.symbols('C1')
#exact_rhs_with_c=exact_rhs.subs(sp.symbols('C1'),C1) #substituting C1 for solve
#ics_eq=exact_rhs_with_c.subs(t,a)-w #initial condition y(a)-y0,substitute t=a
#c_val=sp.solve(ics_eq,C1)[0] #solving y(a)-y0=0
#exact_sol=exact_rhs.subs(sp.symbols('C1'),c_val) #substitute the value of C1 in exact solution
#exact_fun=sp.lambdify(t,exact_sol,'numpy') #parsing into numerical function of t only

#exact_vals=exact_fun(t0)
#print(exact_vals)
#errs=np.abs(exact_vals-w0)
#print(errs)
#sometimes returns value Zero which is not float but string

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



