import numpy as np
import matplotlib.pyplot as plt
from poisson_ib import solve_poisson_ib, solve_poisson_ib_optimized
import time

u_exact = lambda x, y: np.cos(2*x) * np.sin(2*y)
f  = lambda x, y: -8 * np.cos(2*x) * np.sin(2*y)

NxVals = np.array([16, 32, 64, 128])
times_optimized, time_normal = [], []
for N in NxVals:
    start = time.perf_counter()
    X, Y, U, q, xib, yib, info = solve_poisson_ib_optimized(0, 2*np.pi, 0, 2*np.pi, N, N, f, u_exact, u_exact, prescribed_tol=1e-10)
    times_optimized.append(time.perf_counter() - start)
    U_exact = u_exact(X, Y)
    E = U - U_exact

for N in NxVals:
    start = time.perf_counter()
    X, Y, U, q, xib, yib, info = solve_poisson_ib(0, 2*np.pi, 0, 2*np.pi, N, N, f, u_exact, u_exact)
    time_normal.append(time.perf_counter() - start)
    U_exact = u_exact(X, Y)
    E = U - U_exact

print("Times for optimized version:", times_optimized)
print("Times for normal version:", time_normal)