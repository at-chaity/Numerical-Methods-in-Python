import numpy as np
def spline_pol(x,a,p):
    n=len(x)-1
    h=np.zeros(n+1)
    alpha=h.copy()
    l=h.copy()
    u=h.copy()
    z=h.copy()
    c=np.zeros(n+1)
    d=np.zeros(n)
    b=d.copy()
    h=np.diff(x)
    for i in range(1,n):
        alpha[i]=((3/h[i])*(a[i+1]-a[i]))-((3/h[i-1])*(a[i]-a[i-1]))
    l[0]=1
    for i in range(1,n):
        l[i]=2*(x[i+1]-x[i-1])-(h[i-1]*u[i-1])
        u[i]=h[i]/l[i]
        z[i]=(alpha[i]-(h[i-1]*u[i-1]))/l[i]
    l[n]=1
    for j in range(n-1,-1,-1):
        c[j]=z[j]-(u[j]*c[j+1])
        b[j]=((a[j+1]-a[j])/h[j])-(h[j]*(c[j+1]+2*c[j])/3)
        d[j]=(c[j+1]-c[j])/(3*h[j])

    #Evaluate the cubic spline at a given point
    result = []
    for i in range(len(x)-1):
        if x[i] <= p <= x[i+1]:
            dx = p - x[i]
            y = a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
            result.append(y)
            break
    return np.array(result),a,b,c,d

