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

def local_window(xk, yk, x, y, dx, dy, tol):
    s = 1.2 * dx
    R = s * np.sqrt(2 * np.log(1 / (2 * np.pi * s**2 * tol)))
    i0 = max(int(np.ceil((xk - R - x[0]) / dx)), 0)
    i1 = min(int(np.floor((xk + R - x[0]) / dx)) + 1, len(x))
    j0 = max(int(np.ceil((yk - R - y[0]) / dy)), 0)
    j1 = min(int(np.floor((yk + R - y[0]) / dy)) + 1, len(y))
    return i0, i1, j0, j1

def spreading_operator_optimized(q, xib, yib, X, Y, dx, dy, prescribed_tol=1e-10):
    x, y = X[:, 0], Y[0, :]
    S = np.zeros(X.shape)
    for k in range(len(q)):
        i0, i1, j0, j1 = local_window(xib[k], yib[k], x, y, dx, dy, prescribed_tol)
        S[i0:i1, j0:j1] += q[k] * delta_a(X[i0:i1, j0:j1] - xib[k],
                                          Y[i0:i1, j0:j1] - yib[k], dx)
    return S

def J_operator_optimized(u, xib, yib, X, Y, dx, dy, prescribed_tol=1e-10):
    x, y = X[:, 0], Y[0, :]
    J = np.zeros(len(xib))
    for k in range(len(xib)):
        i0, i1, j0, j1 = local_window(xib[k], yib[k], x, y, dx, dy, prescribed_tol)
        J[k] = np.sum(u[i0:i1, j0:j1] * delta_a(X[i0:i1, j0:j1] - xib[k],
                                                Y[i0:i1, j0:j1] - yib[k], dx)) * dx * dy
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

def Amult_optimized(unknowns, L, xib, yib, X, Y, dx, dy, prescribed_tol = 0):
    n = X.size
    u = unknowns[:n]
    q = unknowns[n:]
    top = L @ u - spreading_operator_optimized(q, xib, yib, X, Y, dx, dy, prescribed_tol).ravel()
    bottom = J_operator_optimized(u.reshape(X.shape), xib, yib, X, Y, dx, dy, prescribed_tol)
    return np.concatenate([top, bottom])

def laplacian(Nx, Ny, dx, dy):
    nx, ny = Nx - 1, Ny - 1 
    Ix, Iy = sp.identity(nx), sp.identity(ny)
    return (sp.kron(D2(nx, dx), Iy) + sp.kron(Ix, D2(ny, dy))).tocsr()

def D2(m, h):
        return sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(m, m)) / h**2


def solve_poisson_ib(a, b, c, d, Nx, Ny, f_fn, g_box, g_ib, tol=1e-10, alpha=1.0):
    dx = (b - a) / Nx
    dy = (d - c) / Ny
    x = np.linspace(a, b, Nx + 1)
    y = np.linspace(c, d, Ny + 1)
    X, Y = np.meshgrid(x, y, indexing='ij')

    Xi, Yi = X[1:-1, 1:-1], Y[1:-1, 1:-1]
    xi, yi = x[1:-1], y[1:-1]
    n = Xi.size

    xib, yib = make_ib_points(nib_for_spacing(dx, alpha=alpha))
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

def solve_poisson_ib_optimized(a, b, c, d, Nx, Ny, f_fn, g_box, g_ib, prescribed_tol = 0, tol=1e-10):
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
               matvec=lambda v: Amult_optimized(v, L, xib, yib, Xi, Yi, dx, dy, prescribed_tol))
    sol, info = spla.gmres(A_op, rhs, rtol=tol, atol=0.0, restart=1000, maxiter=1)

    U = g_box(X, Y).copy()
    U[1:-1, 1:-1] = sol[:n].reshape(Xi.shape)
    q = sol[n:]
    return X, Y, U, q, xib, yib, info