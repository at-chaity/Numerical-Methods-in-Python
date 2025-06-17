#Solution of system of Linear equations with LU decomposition
#a1: The coefficient matrix A,(e.g.[[a,b],[c,d]])= ")
#b1: The rhs vector B,(e.g.[p,q])= ")

import numpy as np
def lu_decom(a1,b1):
    a=np.array(a1)
    b=np.array(b1)
    if a.shape[0]!=a.shape[1] or a.shape[0]!=b.shape[0]:
        print("Matrix A and vector B must have same row count")
    else:
        n=len(b)
        L=np.zeros([n,n])
        U=L.copy()
    #Doolittle's algorithm (assuming A doesn't need pivoting)
        for i in range(n):
            for j in range(i,n):
                U[i][j]=a[i][j]-sum(L[i][k]*U[k][j] for k in range(i)) #compute U
            L[i][i]=1
            for j in range(i+1,n):
                L[j][i]=(a[j][i]-sum(L[j][k]*U[k][i] for k in range(i)))/U[i][i]
    #print(U)
    #print(L)
    #forward substitution for LY=B
        y=np.zeros(n)
        for i in range(n):
            y[i]=b[i]-np.dot(L[i,:i],y[:i])
    #backward substitution for UX=Y
        x=np.zeros(n)
        for i in range(n-1,-1,-1):
            x[i]=(y[i]-np.dot(U[i,i+1:],x[i+1:]))/U[i][i]

        return x

