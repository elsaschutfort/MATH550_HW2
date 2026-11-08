# HW 2: Working with Immersed Boundaries

## Problem 1: Poisson Equation

To find $f$ we can take the second derivative of $u(x, y)$ with respect to x and y. 

$$\frac{\partial^2 u}{\partial x^2} = -4\cos{(2x)}\sin{(2y)} \qquad \qquad \frac{\partial^2 u}{\partial y^2} = -4\cos{(2x)}\sin{(2y)}$$

Thus 

$$f = -8\cos{(2x)}\sin{(2y)}$$

When we use the finite difference method to solve the Poisson equation we get 

| ![Descriptive alt text](figures/NumericalSolutionToPE.png) |
| :-: |
| *Figure 1: Numerical Solution to Poisson Equation for u(x,y) = sin(2y)cos(2y)* |

And the solution converges wuth second order in the grid size $\lambda x$

| ![Descriptive alt text](figures/ConvergenceForPE.png) |
| :-: |
| *Figure 2: Relative error of the 2D Poisson solution with Dirchelet boundary conditions versus N.* |

## Problem 2: Poisson Equation with an Immersed Boundary

| ![Descriptive alt text](figures/Problem2_ImmersedBoundary.png) |
| :-: |
| *Figure 3: Numerical Solution to Poisson Equation for u(x,y) = sin(2y)cos(2y) with Immersed Boundary. * |

| ![Descriptive alt text](figures/Problem2_ConvergencePlot.png) |
| :-: |
| *Figure 4: Convergence of the Poisson solution with an immersed boundary. Relative $L_{\infty}$ and $L_2$ errors against N, with $O^N^{-2}$ and $O^N^{-1.5}$ reference slopes.* |

| ![Descriptive alt text](figures/Problem2_ImmersedBoundaryHeatMap.png) |
| :-: |
| *Figure 5: Top voew of the numerical solution $u(x,y) = \lamba^2 u = -8\cos{2x}\sin{2y}$ with the immersed boundary (red circle radius 1 centered at $(\pi, \pi)$)), on a 128x128 grid. The solution matches the test solution.* |


*Once your code works, experiment with this spacing. What happens when the IB points are very far apart compared to the finite difference mesh size? What happens when they are very close together? Look at both the accuracy (Problem 3 will give a better idea of this) and the speed of your code?*

## Problem 3: An almost minimal surface

## Problem 4: Charlie's rule