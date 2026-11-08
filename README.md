# MATH 550: Numerical PDEs - Homework 2

## Project Structure

```
MATH550-HW0/
├── README.md
├── code/
│   ├── problem1.py
|   ├── problem2.py
|   ├── problem3.py
|   ├── problem4.py
|   ├── alpha_comparison.py
│   └── poisson_ib.py
└── writeup/
    ├── README.md
    └── figures/
        ├── ConvergenceForPE.png
        ├── NumericalSolutionToPE.png
        ├── Problem2_ConvergencePlot.png
        ├── Problem2_ImmersedBoundary.png
        ├── Problem2_ImmersedBoundaryHeatMap.png
        ├── Problem3_PotentialPE.png
        ├── Problem4_RunTimeHighTol.png
        └── Problem4_RunTimeLowTol.png
```

## Files

| Path | Description |
|------|-------------|
| `writeup/README.md` | Homework submission |
| `code/problem1.py` | Solves the Poisson equation without an immersed boundary using the 5-point finite-difference Laplacian with Dirichlet boundary conditions, and plots the convergence for the test solution $u = \cos (2x) \sin (2y)$ |
| `code/problem2.py` | Solves the Poisson equation with the circular immersed boundary. |
| `code/problem3.py` | Solves an almost minimal surface problem. |
| `code/problem4.py` | Compares run time and accuracy of the original and optimized spreading/interpolation operators |
| `code/alpha_comparison.py` | Experiments with spacing for Problem 2 and Problem 3. |
| `code/poisson_ib.py` | Core solver. |
| `writeup/figures/ConvergenceForPE.png` | Problem 1: relative $L_\infty$ and $L_2$ error vs. $N$ with an $O(N^{-2})$ reference slope |
| `writeup/figures/NumericalSolutionToPE.png` | Problem 1: 3D surface plot of the numerical solution |
| `writeup/figures/Problem2_ConvergencePlot.png` | Problem 2: relative $L_\infty$ and $L_2$ error vs. $N$ with an $O(N^{-2})$ reference slope |
| `writeup/figures/Problem2_ImmersedBoundary.png` | Problem 2: 3D surface plot of the solution with the immersed boundary curve drawn in red |
| `writeup/figures/Problem2_ImmersedBoundaryHeatMap.png` | Problem 2: top-view of the solution with the immersed boundary circle |
| `writeup/figures/Problem3_PotentialPE.png` | Problem 3: 3D surface plot of the almost minimal surface with the immersed boundary |
| `writeup/figures/Problem4_RunTimeHighTol.png` |  Problem 4: run time vs. $N$ for the original and optimized solvers at prescribed tolerance $10^{-4}$. |
| `writeup/figures/Problem4_RunTimeLowTol.png` |  Problem 4: run time vs. $N$ for the original and optimized solvers at prescribed tolerance $10^{-10}$. |

## Running the Code

```bash 
pip install numpy scipy matplotlib
pip install cmocean
```