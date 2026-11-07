import numpy as np
import matplotlib.pyplot as plt
from poisson_ib import solve_poisson_ib

u_exact = lambda x, y: np.cos(2*x) * np.sin(2*y)
f  = lambda x, y: -8 * np.cos(2*x) * np.sin(2*y)

NxVals = np.array([16, 32, 64, 128])
errInf, errL2 = [], []

for N in NxVals:
    X, Y, U, q, xib, yib, info = solve_poisson_ib(0, 2*np.pi, 0, 2*np.pi, N, N, f, u_exact, u_exact)
    U_exact = u_exact(X, Y)
    E = U - U_exact
    errInf.append(np.max(np.abs(E)) / np.max(np.abs(U_exact)))
    errL2.append(np.linalg.norm(E) / np.linalg.norm(U_exact))

X_last, Y_last, U_last, xib_last, yib_last = X, Y, U, xib, yib
errInf = np.array(errInf)
errL2 = np.array(errL2)

fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(projection='3d', computed_zorder=False)
surf = ax1.plot_surface(X_last, Y_last, U_last, cmap='viridis',
                        edgecolor='none', antialiased=True, rcount=129, ccount=129)
xc_ = np.append(xib_last, xib_last[0]); yc_ = np.append(yib_last, yib_last[0])
ax1.plot(xc_, yc_, u_exact(xc_, yc_), 'r-', linewidth=2, zorder=10,
         label='Immersed boundary')
fig1.colorbar(surf, ax=ax1, shrink=0.6, pad=0.08, label='u(x,y)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('u(x,y)')
ax1.view_init(elev=30, azim=-45)
ax1.set_title('Numerical Solution to Poisson Equation with Immersed Boundary')
ax1.legend(loc='upper left')
plt.tight_layout()
plt.show()

slope_inf = np.polyfit(np.log(NxVals), np.log(errInf), 1)[0]
slope_l2 = np.polyfit(np.log(NxVals), np.log(errL2), 1)[0]

fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.loglog(NxVals, errInf, 'o-', color='blue', linewidth=3, markersize=8,
           label=r'$L_\infty$ Error Slope $\approx$ %0.2f' % slope_inf)
ax2.loglog(NxVals, errL2, 's-', color='orange', linewidth=3, markersize=8,
           label=r'$L_2$ Error Slope $\approx$ %0.2f' % slope_l2)
ax2.loglog(NxVals, errInf[0] * (NxVals / NxVals[0]) ** (-2.0), '--', color='black',
           linewidth=2.5, label=r'Slope $O(N^{-2})$')
ax2.loglog(NxVals, errInf[0] * (NxVals / NxVals[0]) ** (-1.5), ':', color='black',
           linewidth=2.5, label=r'Slope $O(N^{-1.5})$')
ax2.grid(True, which='both')
ax2.set_xlabel('N')
ax2.set_ylabel('Relative Error')
ax2.set_title('Convergence of 2D Poisson Solution with Immersed Boundary')
ax2.legend(loc='lower left')
plt.tight_layout()
plt.show()

fig3, ax3 = plt.subplots(figsize=(7, 6))
pcm = ax3.pcolormesh(X_last, Y_last, U_last, cmap='viridis', shading='auto')
ax3.plot(np.append(xib_last, xib_last[0]), np.append(yib_last, yib_last[0]),
         'r-', linewidth=2, label='Immersed boundary')
fig3.colorbar(pcm, ax=ax3, label='u(x,y)')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_aspect('equal')
ax3.set_title('Numerical Solution to Poisson Equation with Immersed Boundary')
ax3.legend(loc='upper right')
plt.tight_layout()
plt.show()