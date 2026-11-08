import numpy as np
import matplotlib.pyplot as plt
from poisson_ib import solve_poisson_ib, solve_poisson_ib_optimized
import time

u_exact = lambda x, y: np.cos(2*x) * np.sin(2*y)
f  = lambda x, y: -8 * np.cos(2*x) * np.sin(2*y)

NxVals = np.array([16, 32, 64, 128, 256])
prescribed_tol = 1e-10
n_repeats = 2
times_optimized, time_normal = [], []

def timed_solve(solver, N, **kwarfs):
    best_time = np.inf
    for _ in range(n_repeats):
        start_time = time.time()
        output = solver(0, 2*np.pi, 0, 2*np.pi, N, N, f, u_exact, u_exact, **kwarfs)
        end_time = time.time()
        best_time = min(best_time, end_time - start_time)
    return output, best_time

time_normal, time_optimized = [], []
diff_solutions = []

for N in NxVals:
    (X, Y, U_n, q_n, xib, yib, info_n), t_n = timed_solve(solve_poisson_ib, N)
    (X, Y, U_o, q_o, xib, yib, info_o), t_o = timed_solve(solve_poisson_ib_optimized, N,
                                                         prescribed_tol=prescribed_tol)
    U_exact = u_exact(X, Y)

    time_normal.append(t_n)
    time_optimized.append(t_o)
    diff_solutions.append(np.max(np.abs(U_o - U_n)))

time_normal, time_optimized = np.array(time_normal), np.array(time_optimized)
speedup = time_normal / time_optimized

print(f"\n{'N':>5} {'t normal (s)':>13} {'t optimized (s)':>16} {'speed-up':>9} "
      f"{'max|U_opt - U_normal|':>22}")
for i, N in enumerate(NxVals):
    print(f"{N:5d} {time_normal[i]:13.4f} {time_optimized[i]:16.4f} {speedup[i]:9.1f} "
          f"{diff_solutions[i]:22.2e}")


fig2, ax2 = plt.subplots(figsize=(7, 5))
slope_normal = np.polyfit(np.log(NxVals), np.log(time_normal), 1)[0]
slope_opt = np.polyfit(np.log(NxVals), np.log(time_optimized), 1)[0]
ax2.loglog(NxVals, time_normal, 'o-', color='blue', linewidth=3, markersize=8,
           label=r'Normal, slope $\approx$ %0.2f' % slope_normal)
ax2.loglog(NxVals, time_optimized, 's-', color='orange', linewidth=3, markersize=8,
           label=r'Optimized, slope $\approx$ %0.2f' % slope_opt)
for N, t, s in zip(NxVals, time_optimized, speedup):
    ax2.annotate(f'{s:.1f}x', (N, t), textcoords='offset points', xytext=(0, -18),
                 ha='center', fontsize=10)
ax2.grid(True, which='both')
ax2.set_xlabel('N')
ax2.set_ylabel('Wall-clock time (s)')
ax2.set_title(f'Run Time (prescribed_tol = {prescribed_tol:g})')
ax2.legend(loc='upper left')

plt.tight_layout()
plt.show()