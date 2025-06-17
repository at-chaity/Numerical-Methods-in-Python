#Shooting method for linear boundary value problems
#p(x)y''+q(x)y'+r(x)=y'';a<=x<=b,y(a)=c,y(b)=d
#p1: p(x) function
#q1: q(x) function
#r1: r(x) function
#a1: (float) a<=x<=b
#b1: (float) a<=x<=b
#c1: Boundary condition 1,value of y(a)=c
#d1: Boundary condition 1,value of y(b)=d
#n : The number of subintervals

import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp
def lin_shooting(p1,q1,r1,a1,b1,c1,d1,n):
    x=sp.symbols('x')
    symbolic_p=sp.sympify(p1)
    p=sp.lambdify(x,symbolic_p,'numpy')
    symbolic_q=sp.sympify(q1)
    q=sp.lambdify(x,symbolic_q,'numpy')
    symbolic_r=sp.sympify(r1)
    r=sp.lambdify(x,symbolic_r,'numpy')

#a,b,alpha,beta in float
    def dec_in(v_str):
        return float(sp.sympify(v_str))
    a=dec_in(a1)
    b=dec_in(b1)
    alpha=dec_in(c1)
    beta=dec_in(d1)
    #exact solution
    def fun(x, y):
        y1, y2 = y[0], y[1]
        dy1 = y2
        dy2 = p(x) * y2 + q(x) * y1 + r(x)
        return np.vstack((dy1, dy2))

    # Boundary conditions: y(1) = 1, y(2) = 2
    def bc(ya, yb):
        return np.array([ya[0] - alpha, yb[0] - beta])

    x2 = np.linspace(a,b,n+1)
    y_guess = np.zeros((2, x2.size))  # [y, y']
    y_guess[0] = np.linspace(alpha, beta, x2.size)  # crude guess for y

    # Solve BVP
    sol = solve_bvp(fun, bc, x2, y_guess)
    #shooting method
    h=(b-a)/n
    u1 = np.zeros(n + 1)
    u2 = np.zeros(n + 1)
    v1 = np.zeros(n + 1)
    v2 = np.zeros(n + 1)
    w1 = np.zeros(n + 1)
    w2 = np.zeros(n + 1)
    u1[0]=alpha
    v2[0]=1
    for i in range(1,n+1):
        x1=a+(i*h)
        p1=p(x1)
        p2=p(x1+(h/2))
        p3=p(x1+h)
        q1=q(x1)
        q2=q(x1+(h/2))
        q3=q(x1+h)
        r1=r(x1)
        r2=r(x1+(h/2))
        r3=r(x1+h)
        k11=h*u2[i-1]
        k12=h*((p1*u2[i-1])+(q1*u1[i-1])+r1)
        k21=h*(u2[i-1]+(k12/2))
        k22=h*((p2*(u2[i-1]+(k12/2)))+(q2*(u1[i-1]+(k11/2)))+r2)
        k31=h*(u2[i-1]+(k22/2))
        k32=h*((p2*(u2[i-1]+(k22/2)))+(q2*(u1[i-1]+(k21/2)))+r2)
        k41=h*(u2[i-1]+k32)
        k42=h*((p3*(u2[i-1]+k32))+(q3*(u1[i-1]+k31))+r3)
        u1[i]=u1[i-1]+(k11+k21+k31+k41)/6
        u2[i]=u2[i-1]+(k12+k21+k32+k42)/6
        k_11 = h * v2[i-1]
        k_12 = h * ((p1 * v2[i-1]) + (q1 * v1[i-1]))
        k_21 = h * (v2[i-1] + (k_12 / 2))
        k_22 = h * ((p2 * (v2[i-1] + (k_12 / 2))) + (q2 * (v1[i-1] + (k_11 / 2))))
        k_31 = h * (v2[i-1] + (k_22 / 2))
        k_32 = h * ((p2 * (v2[i-1] + (k_22 / 2))) + (q2 * (v1[i-1] + (k_21 / 2))))
        k_41 = h * (v2[i-1] + k_32)
        k_42 = h * ((p3 * (v2[i-1] + k_32)) + (q3 * (v1[i-1] + k_31)))
        v1[i ] = v1[i-1] + (k_11 + k_21 + k_31 + k_41) / 6
        v2[i ] = v2[i-1] + (k_12 + k_21 + k_32 + k_42) / 6
    w1[0]=alpha
    w2[0]=(beta-u1[n])/v1[n]
    for i in range(1,n+1):
        w1[i]=u1[i]+w2[0]*v1[i]
        w2[i]=u2[i]+w2[0]*v2[i]
    return x2,w1,w2,sol.y[0],sol.y[1]