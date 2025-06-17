# Successive Over Relaxation(SOR) method for solving linear system of equations

#a1: The coefficient matrix A (e.g.[[a,b],[c,d]])
#b1: The rhs vector B (e.g.[p,q])
#x1: The initial guess of the solutions,x0=[x1,x2]
#w:  Relaxation parameter (w>1 for over relaxation)
#tol:(float) Tolerance for convergence

import numpy as np
def sor(a1,b1,x1,w,tol,return_maxiter=False):
    a=np.array(a1,dtype=float)
    b=np.array(b1,dtype=float)
    xo=np.array(x1,dtype=float)
    if a.shape[0]!=a.shape[1] or a.shape[0]!=b.shape[0]:
        print("Matrix A and vector B must have same row count")
    else:
        N=1000
        n=len(b)
        x=np.copy(xo)
        for k in range(N):
            for i in range(n):
                 sum1=0
                 sum2=0
                 for j in range(i):
                     sum1 +=np.sum(a[i][j]*x[j])
                 for j in range(i+1,n):
                     sum2 += np.sum(a[i][j] * xo[j])
                 x[i]=(1-w)*xo[i]+w*(b[i]-sum2-sum1)/a[i][i]
            if np.linalg.norm(x-xo,ord=np.inf)<tol:
                if return_maxiter:
                    return x, k + 1
                else:
                    return x
            else:
                xo = x.copy()
        else:
            return xo
