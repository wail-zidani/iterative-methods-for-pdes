# Iterative methods for PDEs

Numerical resolution of partial differential equations using *finite differences*
and *iterative methods*.

This repository contains academic implementations developed during the
*M1 MACS (Modélisation, Analyse et Calcul Scientifique)*  
Nantes Université & Centrale Nantes.

---

## Contents

### 1D Poisson equation
\[
•⁠  ⁠u''(x) = f(x), \quad x \in (0,1), \quad u(0)=u(1)=0
\]

Discretization by finite differences.

### Implemented methods
•⁠  ⁠*Jacobi method*
•⁠  ⁠*Gauss–Seidel method*

Each solver computes:
•⁠  ⁠convergence iterations
•⁠  ⁠relative residual
•⁠  ⁠numerical error vs exact solution

---

## Languages & tools
•⁠  ⁠Python (NumPy)
•⁠  ⁠Numerical linear algebra
•⁠  ⁠Iterative solvers for PDEs

---

## Academic context
This work is part of my training in **Applied Mathematics, Numerical Analysis
and Scientific Computing (HPC-oriented)**.
## Problem
We consider the 1D Poisson equation
-u''(x) = f(x),  x ∈ (0,1)
u(0) = u(1) = 0

Exact solution used for validation:
u(x) = sin(πx)

## Implemented methods
•⁠  ⁠Jacobi iterative method
•⁠  ⁠Gauss–Seidel iterative method

## Languages
•⁠  ⁠Python (NumPy)
•⁠  ⁠Fortran 90 (scientific computing / HPC)

## How to run

### Python
```bash
python python/jacobi_1d.py
python python/gauss_seidel_1d.py
gfortran -O2 fortran/jacobi_1d_poisson.f90 -o jacobi
./jacobi

gfortran -O2 fortran/gauss_seidel_1d_poisson.f90 -o gs
./gs
