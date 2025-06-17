#Period method for eigenvalue

#a1: The coefficient matrix A (e.g.[[a,b],[c,d]])
#x1: The initial guess of the solutions,x0=[x1,x2]
#tol: (float) Tolerance for convergence

import numpy as np
def power(a1,x1,tol):
    a=np.array(a1,dtype=float)
    xo=np.array(x1,dtype=float)
    n=len(xo)
    N=1000
    eigen_val=0
#power method
    for k in range(1,N):
        x=np.matmul(a,xo)
        eigen_val_new=np.dot(x,xo)/np.dot(xo,xo) #Rayleigh quotient for better accuracy
        eigen_vec=x/np.linalg.norm(x)
        if abs(eigen_val_new-eigen_val)<tol:
            return eigen_val_new,eigen_vec
        else:
            xo=eigen_vec
            eigen_val=eigen_val_new
    else:
        return eigen_val,x1






