import numpy as np
import matplotlib.pyplot as plt
from poisson_ib import solve_poisson_ib

g_ib = lambda x, y: np.cos(2*x) * np.sin(2*y)
g_box = lambda x, y: 0*x
f = lambda x, y: 0*x

N = 128

X, Y, U, q, xib, yib, info = solve_poisson_ib(0, 2*np.pi, 0, 2*np.pi, N, N, f, g_box, g_ib)


fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(projection='3d', computed_zorder=False)
surf = ax1.plot_surface(X, Y, U, cmap='viridis',
                        edgecolor='none', antialiased=True, rcount=129, ccount=129)
xc_ = np.append(xib, xib[0]); yc_ = np.append(yib, yib[0])
ax1.plot(xc_, yc_, g_ib(xc_, yc_), 'r-', linewidth=2, zorder=10,
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
