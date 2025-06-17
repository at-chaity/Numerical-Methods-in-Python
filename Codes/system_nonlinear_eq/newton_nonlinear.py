import numpy as np
import sympy as sp
from sympy import symbols,sympify,cos,sin,tan,exp,sqrt,log

def newton_solver(n,x,f_sym,f_np,x0,tol):
    N=1000
    err=[]
    for k in range(N):
        #evaluate Fs
        f_val=np.array([f(*x0) for f in f_np],dtype=float)
        #evaluate jacobian
        j_sym=[[sp.diff(f,xi) for xi in x]for f in f_sym]
        j_np=[[sp.lambdify(x,df,modules='numpy') for df in row] for row in j_sym]
        n1=len(j_np)
        j_eval=np.zeros((n1,n1))
        for i in range(n1):
            for j in range(n1):
                j_eval[i,j] +=j_np[i][j](*x0)
        #solve JY=-F
        try:
            y=np.linalg.solve(j_eval,-f_val)
        except np.linalg.LinAlgError:
            print('Jacobian is singular')
            return
        x1=x0+y
        error=np.linalg.norm(y,ord=2)
        err.append(error)
        if error<tol:
            return x1,err
        else:
            x0=x1
    else:
        print('No convergence')






