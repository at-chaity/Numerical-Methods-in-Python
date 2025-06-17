import numpy as np
#input
a1=eval(input("Enter coefficient matrix A,(e.g.[[a,b],[c,d]])= "))
#convverting input into arrays
a=np.array(a1,dtype=float)
m,n=a.shape #size of a
ata=a.T@a #turning a into square matrix for eigenvalues
eigen_vals,V=np.linalg.eig(ata) #determine eigenvalue
#indexing eigenvalues in descending order
idx=np.argsort(eigen_vals)[::-1]
eigen_vals=eigen_vals[idx]
V=V[:,idx]

singular_values=np.sqrt(eigen_vals) #s matrix
#construct sigma matrix
sigma=np.zeros_like(a)
np.fill_diagonal(sigma,singular_values)
#construct u matrix, u=av/sigma
U=np.zeros((m,m))
for i in range(len(singular_values)):
    U[:,i]=a@V[:,i]/singular_values[i]
U,_=np.linalg.qr(U)
a_recon=U @ sigma @ V.T

print("U= ",U)
print("Sigma= ",sigma)
print("V= ",V)
print("a_recon= ",a_recon)