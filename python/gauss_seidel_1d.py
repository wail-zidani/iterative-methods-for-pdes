import numpy as np


def build_poisson_1d(N: int):
    """
    -u''(x) = f(x) on (0,1), u(0)=u(1)=0
    Finite differences on interior points.
    Returns A (NxN), b (N,), grid x (N,), and h.
    """
    h = 1.0 / (N + 1)
    x = np.linspace(h, 1.0 - h, N)

    A = np.zeros((N, N))
    diag = 2.0 / (h * h)
    off = -1.0 / (h * h)

    for i in range(N):
        A[i, i] = diag
        if i > 0:
            A[i, i - 1] = off
        if i < N - 1:
            A[i, i + 1] = off

    return A, x, h


def jacobi(A, b, x0=None, tol=1e-10, max_iter=50000):
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()

    D = np.diag(A)
    R = A - np.diag(D)

    b_norm = np.linalg.norm(b)
    if b_norm == 0:
        b_norm = 1.0

    for k in range(max_iter):
        x_new = (b - R @ x) / D
        r = b - A @ x_new
        rel_res = np.linalg.norm(r) / b_norm
        x = x_new
        if rel_res < tol:
            return x, k + 1, rel_res

    r = b - A @ x
    rel_res = np.linalg.norm(r) / b_norm
    return x, max_iter, rel_res


def gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=50000):
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()

    b_norm = np.linalg.norm(b)
    if b_norm == 0:
        b_norm = 1.0

    for k in range(max_iter):
        for i in range(n):
            # sum_{j<i} A_ij x_j  +  sum_{j>i} A_ij x_j
            s1 = A[i, :i] @ x[:i]
            s2 = A[i, i+1:] @ x[i+1:]
            x[i] = (b[i] - s1 - s2) / A[i, i]

        r = b - A @ x
        rel_res = np.linalg.norm(r) / b_norm
        if rel_res < tol:
            return x, k + 1, rel_res

    r = b - A @ x
    rel_res = np.linalg.norm(r) / b_norm
    return x, max_iter, rel_res


def main():
    N = 200

    # Example with known exact solution u(x)=sin(pi x) => -u'' = pi^2 sin(pi x)
    A, x, h = build_poisson_1d(N)
    f = (np.pi ** 2) * np.sin(np.pi * x)
    b = f.copy()

    u_exact = np.sin(np.pi * x)

    u_j, it_j, res_j = jacobi(A, b, tol=1e-10, max_iter=200000)
    err_j = np.max(np.abs(u_j - u_exact))

    u_gs, it_gs, res_gs = gauss_seidel(A, b, tol=1e-10, max_iter=200000)
    err_gs = np.max(np.abs(u_gs - u_exact))

    print("1D Poisson (FD) | u_exact = sin(pi x)")
    print(f"N = {N}, h = {h:.3e}")
    print("-" * 55)
    print(f"Jacobi      : iterations = {it_j:6d} | rel_res = {res_j:.2e} | max_err = {err_j:.2e}")
    print(f"Gauss-Seidel: iterations = {it_gs:6d} | rel_res = {res_gs:.2e} | max_err = {err_gs:.2e}")


if _name_ == "_main_":
    main()
