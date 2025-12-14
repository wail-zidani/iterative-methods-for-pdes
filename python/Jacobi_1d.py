"""
Jacobi solver for 1D Poisson problem:
    -u''(x) = f(x) on (0,1), with u(0)=u(1)=0

Finite differences (2nd order) on an interior grid of size N.
Jacobi iterative method.

Author: Wail Zidani
M1 MACS – Nantes Université & Centrale Nantes
"""

import numpy as np


def build_poisson_1d(N, f):
    h = 1.0 / (N + 1)
    x = np.linspace(h, 1.0 - h, N)

    main = 2.0 * np.ones(N)
    off = -1.0 * np.ones(N - 1)
    A = (1.0 / h**2) * (np.diag(main) + np.diag(off, 1) + np.diag(off, -1))

    b = f(x)
    return A, b, x


def jacobi(A, b, tol=1e-10, max_iter=100000):
    n = len(b)
    x = np.zeros(n)

    D = np.diag(A)
    R = A - np.diag(D)

    b_norm = np.linalg.norm(b)
    if b_norm == 0:
        b_norm = 1.0

    for k in range(max_iter):
        x_new = (b - R @ x) / D
        r = b - A @ x_new
        rel_res = np.linalg.norm(r) / b_norm

        if rel_res < tol:
            return x_new, k + 1, rel_res

        x = x_new

    return x, max_iter, rel_res


def main():
    N = 200
    f = lambda x: (np.pi**2) * np.sin(np.pi * x)
    u_exact = lambda x: np.sin(np.pi * x)

    A, b, x = build_poisson_1d(N, f)
    u_num, iters, res = jacobi(A, b)

    error = np.max(np.abs(u_num - u_exact(x)))

    print("Jacobi method – 1D Poisson")
    print(f"N = {N}")
    print(f"Iterations = {iters}")
    print(f"Relative residual = {res:.2e}")
    print(f"Max error = {error:.2e}")


if _name_ == "_main_":
    main()
