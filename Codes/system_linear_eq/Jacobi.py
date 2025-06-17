#Jacobi Iteration method for solving system of Linear equations
#works well with larger systems than gaussian elimination

#a1: The coefficient matrix A (e.g.[[a,b],[c,d]])
#b1: The rhs vector B (e.g.[p,q])
#x1: The initial guess of the solutions,x0=[x1,x2]
#tol: (float) Tolerance for convergence
import numpy as np
def jacobi(a1,b1,x1,tol,return_maxiter=False):
    a=np.array(a1,dtype=float)
    b=np.array(b1,dtype=float)
    xo=np.array(x1,dtype=float)
    if a.shape[0]!=a.shape[1] or a.shape[0]!=b.shape[0]:
        print("Matrix A and vector B must have same row count")
    else:
        N=1000
        n=len(b)
        x=np.zeros(n)
        for k in range(N):
            for i in range(n):
                sum=0
                for j in range(n):
                    if j!=i:
                        sum+=np.sum(a[i][j]*xo[j])
                x[i]=(b[i]-sum)/a[i][i]
            if np.linalg.norm(x-xo,ord=np.inf)<tol:
                if return_maxiter:
                    return x,k+1
                else:
                    return x
            else:
                xo=x.copy()
        else:
            return xo
