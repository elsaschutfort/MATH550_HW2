import numpy as np
import matplotlib.pyplot as plt
from poisson_ib import solve_poisson_ib, J_operator
import time

g_ib = lambda x, y: np.cos(2*x) * np.sin(2*y)
g_box = lambda x, y: 0*x
f = lambda x, y: 0*x

N = 64
alphas = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0])

dx = dy = (2*np.pi) / N
th = np.linspace(0, 2*np.pi, 400, endpoint=False)
xc, yc = np.pi + np.cos(th), np.pi + np.sin(th)
"""
for alpha in alphas:
    start = time.time()
    X, Y, U, q, xib, yib, info = solve_poisson_ib(0, 2*np.pi, 0, 2*np.pi, N, N, f, g_box, g_ib, alpha=alpha)
    end = time.time()
    Xi, Yi = X[1:-1, 1:-1], Y[1:-1, 1:-1]
    u_on_circle = J_operator(U[1:-1, 1:-1], xc, yc, Xi, Yi, dx, dy)
    mismatch = np.max(np.abs(u_on_circle - np.cos(2*xc)*np.sin(2*yc)))
    print(f"alpha = {alpha:.2f}, max|u_on_circle - g_ib| = {mismatch:.2e}, time = {end - start:.4f} s")

"""
u_exact = lambda x, y: np.cos(2*x) * np.sin(2*y)
f  = lambda x, y: -8 * np.cos(2*x) * np.sin(2*y)

errInf, errL2 = [], []

for alpha in alphas:
    start = time.time()
    X, Y, U, q, xib, yib, info = solve_poisson_ib(0, 2*np.pi, 0, 2*np.pi, N, N, f, u_exact, u_exact, alpha=alpha)
    end = time.time()
    U_exact = u_exact(X, Y)
    E = U - U_exact
    errInf = np.max(np.abs(E)) / np.max(np.abs(U_exact))

    print(f"alpha = {alpha:.2f}, max|u_on_circle - g_ib| = {errInf:.2e}, time = {end - start:.4f} s")
