import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

def make_ib_points(Nib, xc=np.pi, yc=np.pi, r=1.0):
    theta = np.linspace(0.0, 2*np.pi, Nib, endpoint=False)
    return xc + r * np.cos(theta), yc + r * np.sin(theta)
 
 
def nib_for_spacing(dx, r=1.0, alpha=1.0):
    return max(3, int(round(2.0 * np.pi * r / (alpha * dx))))

def spreading_operator(q, xib, yib, X, Y, dx):
    S = np.zeros(X.shape)
    for k in range(len(q)):
        S += q[k] * delta_a(X - xib[k], Y - yib[k], dx)
    return S

def J_operator(u, xib, yib, X, Y, dx, dy):
    J = np.zeros(len(xib))
    for k in range(len(xib)):
        J[k] = np.sum(u * delta_a(X - xib[k], Y - yib[k], dx)) * dx * dy
    return J

def delta_a(x, y, dx):
    return 1 / (2 * np.pi * (1.2*dx) ** 2 ) * np.exp(-((x**2 + y**2) / (2 * (1.2 * dx)**2)))

def Amult(unknowns, L, xib, yib, X, Y, dx, dy):
    n = X.size
    u = unknowns[:n]
    q = unknowns[n:]
    top = L @ u - spreading_operator(q, xib, yib, X, Y, dx).ravel()
    bottom = J_operator(u.reshape(X.shape), xib, yib, X, Y, dx, dy)
    return np.concatenate([top, bottom])

def laplacian(Nx, Ny, dx, dy):
    nx, ny = Nx - 1, Ny - 1 
    Ix, Iy = sp.identity(nx), sp.identity(ny)
    return (sp.kron(D2(nx, dx), Iy) + sp.kron(Ix, D2(ny, dy))).tocsr()

def D2(m, h):
        return sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(m, m)) / h**2


def solve_poisson_ib(a, b, c, d, Nx, Ny, f_fn, g_box, g_ib, tol=1e-10):
    dx = (b - a) / Nx
    dy = (d - c) / Ny
    x = np.linspace(a, b, Nx + 1)
    y = np.linspace(c, d, Ny + 1)
    X, Y = np.meshgrid(x, y, indexing='ij')

    Xi, Yi = X[1:-1, 1:-1], Y[1:-1, 1:-1]
    xi, yi = x[1:-1], y[1:-1]
    n = Xi.size

    xib, yib = make_ib_points(nib_for_spacing(dx))
    Nib = len(xib)
    L = laplacian(Nx, Ny, dx, dy)

    F = f_fn(Xi, Yi)
    F[0, :]  -= g_box(a, yi) / dx**2     # next to left wall   x = a
    F[-1, :] -= g_box(b, yi) / dx**2     # next to right wall  x = b
    F[:, 0]  -= g_box(xi, c) / dy**2     # next to bottom wall y = c
    F[:, -1] -= g_box(xi, d) / dy**2     # next to top wall    y = d

    rhs = np.concatenate([F.ravel(), g_ib(xib, yib)])

    A_op = spla.LinearOperator((n + Nib, n + Nib),
               matvec=lambda v: Amult(v, L, xib, yib, Xi, Yi, dx, dy))
    sol, info = spla.gmres(A_op, rhs, rtol=tol, atol=0.0, restart=1000, maxiter=1)

    U = g_box(X, Y).copy()
    U[1:-1, 1:-1] = sol[:n].reshape(Xi.shape)
    q = sol[n:]
    return X, Y, U, q, xib, yib, info