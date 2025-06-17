#Solution of system of Linear equations with Gaussian elimination
#a1: The coefficient matrix A,(e.g.[[a,b],[c,d]])= ")
#b1: The rhs vector B,(e.g.[p,q])= ")

import numpy as np
def gauss_elimination(a1,b1):
    a=np.array(a1)
    b=np.array(b1)
    if a.shape[0]!=a.shape[1] or a.shape[0]!=b.shape[0]:
        print("Matrix A and vector B must have same row count")
    else:
        ab=np.hstack([a,b.reshape(-1,1)]) #augmented matrix
        n=len(b)
        #forward elimination
        for i in range(n):
            pivot_row=i
            max_val=abs(ab[i][i])
            for j in range(i+1,n):
                if np.abs(ab[j][i])>max_val: #checking if other elements are bigger for numerical stability
                 max_val=abs(ab[j][i])
                 pivot_row=j
            if pivot_row!=i:
                ab[[i,pivot_row]]=ab[[pivot_row,i]] #swapping rows
            pivot=ab[i,i]
            if pivot == 0:
                print("No unique solution exists (zero pivot encountered)")
                exit()
            ab[i]=ab[i]/pivot # 1 is pivot position
            for k in range(i+1,n):
                factor=ab[k,i]
                ab[k]=ab[k]-factor*ab[i] # making 0 in lower rows
        #Backward substitution
        x=np.zeros(n)
        for i in range(n-1,-1,-1):
            x[i]=ab[i][-1]
            for j in range(i+1,n):
                x[i]=x[i]-ab[i,j]*x[j]
        return x











