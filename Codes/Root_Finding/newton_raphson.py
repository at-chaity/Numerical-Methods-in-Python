# Newton-Raphson method fo root finding. It doesn't need an interval [a,b] where f(a)*f(b)=0
#f : Function
#p0: Initial guess of the root
#tol(float): Tolerance for error
#N0(int): Maximum number of iterations

from sympy import symbols, diff, lambdify, sympify

def newton_raphson(user_input,p0,tol,N0,return_iters=False,return_history=False):
    x=symbols("x")
    f_fun=sympify(user_input)
    f_fun_prime=diff(f_fun,x)
    f=lambdify(x,f_fun)
    f_prime=lambdify(x,f_fun_prime)
    history=[p0]
    for i in range(N0):
        p=p0-(f(p0)/f_prime(p0))
        history.append(p)
        if abs(p-p0)<=tol:
            if return_history:
                return p,i+1,history
            if return_iters:
                return p,i+1
            return p

        else:
            p0=p
    if return_history:
        return p0, N0, history
    elif return_iters:
        return p0, N0
    else:
        return p0