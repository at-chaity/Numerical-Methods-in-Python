#   Bisection method to find the root of a function in [a,b]
#f : Function
#a(float): Left end point of the Interval
#b(float): Right end point of the Interval
#tol(float): Tolerance for error
#N0(int): Maximum number of iterations
from sympy import symbols,sympify,lambdify
def bisection(user_input,a,b,tol,n0,return_iters=False,return_history=False):
    x = symbols("x")
    f_fun = sympify(user_input)
    f = lambdify(x, f_fun)
    fa=[]
    fp=[]
    if f(a)*f(b)>0: #f(a)*f(b) must be equal to 0
        return 'Error in interval'
    else:
        for i in range(n0):
            fa_val=f(a)
            fa.append(fa_val) #filling the value of f(a)
            p=a+(b-a)/2
            fp_val=f(p)
            fp.append(fp_val)
            if fp_val==0 or ((b-a)/2)<=tol:
                if return_history:
                    return p,i+1,fp
                if return_iters:
                    return p, i+1
                return p
                #print(f"Root is= {p:.5f}")
                #print("Number of iterations= ",i)
            else:
                if fa_val*fp_val>0:
                    a=p
                    fa[i]=fp[i]
                else:
                    b=p
        if return_history:
            return a,n0,fp
        if return_iters:
            return p,n0
        return p
