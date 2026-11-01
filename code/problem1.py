import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

tpi = 2.0 * np.pi

plt.rcParams.update({'font.size': 18})
def idx(i, j, Ny):
    return i * Ny + j

def solve_poisson_2d(a, b, c, d, Nx, Ny, f_fn, u_exact):
    dx = (b - a) / Nx
    dy = (d - c) / Ny
    x = np.linspace(a, b, Nx + 1)
    y = np.linspace(c, d, Ny + 1)
    nx, ny = Nx + 1, Ny + 1
    X, Y = np.meshgrid(x, y, indexing='ij')

    Dx = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(nx, nx)) / dx**2
    Dy = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(ny, ny)) / dy**2

    I_x = sp.eye(nx)
    I_y = sp.eye(ny)
    L = sp.kron(Dx, I_y) + sp.kron(I_x, Dy)
    L = L.tolil() 

    F = -f_fn(X, Y).flatten()

    for i in range(nx):
        for j in range(ny):
            k = idx(i, j, ny)
            if i == 0 or i == nx - 1 or j == 0 or j == ny - 1:
                L[k, :] = 0.0
                L[k, k] = 1.0
                F[k] = u_exact(x[i], y[j])

    L = L.tocsr()
    phi = spla.spsolve(L, F)
    return phi.reshape((nx, ny)), X, Y

a, b = 0.0, tpi
c, d = 0.0, tpi
u_fn = lambda x, y: np.sin(2*y) * np.cos(2*x)
f_fn = lambda x, y: -8.0 * np.sin(2*y) * np.cos(2*x)

errInf, errL2 = [], []

NxVals = np.array([10, 15, 20, 35, 40, 50, 70, 100])
NyVals = NxVals.copy()

phi_last, X_last, Y_last = None, None, None

for Nx, Ny in zip(NxVals, NyVals):
    phi, X, Y = solve_poisson_2d(a, b, c, d, Nx, Ny, f_fn, u_fn)
    u_exact = u_fn(X, Y)

    errInf.append(np.max(np.abs(u_exact - phi)) / np.max(np.abs(u_exact)))
    errL2.append(np.linalg.norm(u_exact - phi) / np.linalg.norm(u_exact))

    phi_last, X_last, Y_last = phi, X, Y 
errInf = np.array(errInf)
errL2 = np.array(errL2)

fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(projection='3d')
surf = ax1.plot_surface(X_last, Y_last, phi_last, cmap='viridis',
                         edgecolor='none', antialiased=True)
fig1.colorbar(surf, ax=ax1, shrink=0.6, pad=0.08, label='u(x,y)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('u(x,y)')
ax1.set_title('Numerical Solution to Poisson Equation')
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
ax2.grid(True, which='both')
ax2.set_xlabel('N')
ax2.set_ylabel('Relative Error')
ax2.set_title('Convergence of 2D Poisson Numerical Solution')
ax2.legend(loc='lower left')
plt.tight_layout()
plt.show()
