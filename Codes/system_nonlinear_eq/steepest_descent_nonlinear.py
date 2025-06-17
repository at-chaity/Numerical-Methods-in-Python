import numpy as np
import sympy as sp
from sympy import symbols, sympify, cos, sin, tan, exp, sqrt, log

def steepest_descent_solver(n, x_sym, f_sym, f_np,x0,tol):
    max_iter = 1000
    # Define g(x) = sum(fi^2)
    g_sym = sum(f**2 for f in f_sym)
    g_np = sp.lambdify(x_sym, g_sym, modules='numpy')

    # Jacobian of F
    j_sym = [[sp.diff(f, xi) for xi in x_sym] for f in f_sym]
    j_np = [[sp.lambdify(x_sym, df, modules='numpy') for df in row] for row in j_sym]
    err=[]
    for k in range(max_iter):
        # Evaluate F(x0)
        f_val = np.array([f(*x0) for f in f_np], dtype=float)

        # Evaluate Jacobian at x0
        J = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                J[i, j] = j_np[i][j](*x0)

        # Compute descent direction: -∇g = -2 * J.T @ F
        z = 2 * J.T @ f_val
        norm_z = np.linalg.norm(z)
        if norm_z == 0:
            print("Zero gradient. May be already at minimum.")
            break
        z = z / norm_z  # normalize

        g1 = g_np(*x0)

        # Line search using polynomial interpolation
        alpha1 = 0
        alpha3 = 1
        g3 = g_np(*list(x0 - alpha3 * z))
        while g3 >= g1:
            alpha3 /= 2
            if alpha3 < 1e-8:
                break
            g3 = g_np(*list(x0 - alpha3 * z))

        alpha2 = alpha3 / 2
        g2 = g_np(*list(x0 - alpha2 * z))

        h1 = (g2 - g1) / alpha2
        h2 = (g3 - g2) / (alpha3 - alpha2)
        h3 = (h2 - h1) / alpha3

        if h3 == 0:
            alpha0 = alpha3 / 2
        else:
            alpha0 = 0.5 * (alpha2 - h1 / h3)

        x1 = x0 - alpha0 * z

        error = np.linalg.norm(x1 - x0)
        err.append(error)
        if error < tol:
            return x1,err
        else:
            x0 = x1
    else:
        print("No convergence.")
