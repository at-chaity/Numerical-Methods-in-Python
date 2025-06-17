# Runge- Kutta method of order  4

#f_in: The function f(t,y), y'= ")
#a1  : Lower bound of t, a<=t<=b
#b1  : Upper bound of t, a<=t<=b
#w1  : The initial condition, y(a)=y0
#n   : The number of steps
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
def rk4(f_in,a1,b1,w1,n,return_main=False,return_error=False):
    t, y = sp.symbols('t y')
    f_symbolic = sp.sympify(f_in)
    f = sp.lambdify([t, y], f_symbolic, 'numpy')

    def dec_in(v_str):
        return float(sp.sympify(v_str))

    a = dec_in(a1)
    b = dec_in(b1)
    w = dec_in(w1)
#RK4 method
    h=(b-a)/n
    t0=np.linspace(a,b,n+1) #issues with solve_ivp
    w0=np.zeros([n+1])
    t0[0]=a
    w0[0]=w
    for i in range(1,n+1):
        k1=h*f(t0[i-1],w0[i-1])
        k2=h*f(t0[i-1]+(h/2),w0[i-1]+(k1/2))
        k3=h*f(t0[i-1]+(h/2),w0[i-1]+(k2/2))
        k4=h*f(t0[i-1]+h,w0[i-1]+k3)
        w0[i]=w0[i-1]+((k1+(2*k2)+(2*k3)+k4)/6)

#scipy exact solution
    sol=solve_ivp(f,[a,b],[w],t_eval=t0) #solves in RK54
    exact_sol=sol.y[0]
    errs=np.abs(exact_sol-w0)
#print("Exact solution= ",exact_sol)

    if return_main:
        if return_error:
            return w0,exact_sol,errs
        else:
            return w0,exact_sol
    else:
        return w0

#visualization
#plt.figure(figsize=(6,6))
#plt.plot(t0,exact_sol,label="Exact solution by solve_ivp",color='red',linewidth=2,marker='o')
#plt.plot(t0,w0,label="RK4 method",color='blue',linestyle='--',marker='x')
#plt.title("Runge-Kutta method vs Exact solution")
#plt.xlabel("t")
#plt.ylabel("y")
#plt.grid(True)
#plt.legend()
#plt.tight_layout()
#plt.show()
