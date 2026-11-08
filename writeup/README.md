# HW 2: Working with Immersed Boundaries

## Problem 1: Poisson Equation

We solve $\nabla^2 u = f(x,y)$ on $[0, 2\pi] \times [0, 2\pi]$ with Dirichlet boundary conditions, choosing $f$ and the boundary conditions so that the test solution is 

$$u(x, y) = \cos (2x) \sin (2y).$$

Differentiating the test solution twice in each direction.

$$\frac{\partial^2 u}{\partial x^2} = -4\cos{(2x)}\sin{(2y)} \qquad \qquad \frac{\partial^2 u}{\partial y^2} = -4\cos{(2x)}\sin{(2y)}$$

so

$$f(x, y) = \nabla^2 u = -8\cos{(2x)}\sin{(2y)}.$$

For the boundary conditions we can evalute the test solution on each edge of the box gives 

$$u(0, y) = u(2\pi, y) = \sin(2y), \qquad u(x, 0) = u(x, 2\pi) = 0.$$

### Results

| ![Descriptive alt text](figures/NumericalSolutionToPE.png) |
| :-: |
| *Figure 1: Numerical Solution to Poisson Equation for u(x,y) = sin(2y)cos(2y)* |

And the solution converges wuth second order in the grid size $\lambda x$

| ![Descriptive alt text](figures/ConvergenceForPE.png) |
| :-: |
| *Figure 2: Relative error of the 2D Poisson solution with Dirchelet boundary conditions versus N.* |

## Problem 2: Poisson Equation with an Immersed Boundary

We placed an unknown charge $q_k$ at each of $N_{ib}$ immeresed boundary points on the circle (spaced $\approx \Delta x$ apart) adn spread them onto teh grid with a Guassian kernel $\delta_a$ of width $\sigma = 1.2 \Delta x$, giving $\nabla^2 u = f + S(q)$. The charges are determined by requiring that $u$, interpolated back to the IB points with the same kernel, equals $\cos(2x^{ib})\sin(2y^{ib})$; the resulting saddle-point system for $\begin{bmatrix} u(x,y) \\ q \end{bmatrix}$ is solved.

### Results

| ![Descriptive alt text](figures/Problem2_ImmersedBoundary.png) |
| :-: |
| *Figure 3: Numerical Solution to Poisson Equation for u(x,y) = sin(2y)cos(2y) with Immersed Boundary. * |

| ![Descriptive alt text](figures/Problem2_ConvergencePlot.png) |
| :-: |
| Figure 4: Convergence of the Poisson solution with an immersed boundary. Relative $L_{\infty}$ and $L_2$ errors against N, with ${O}^{N^{-2}}$ and $O^{N^{-1.5}}$ reference slopes. |

| ![Descriptive alt text](figures/Problem2_ImmersedBoundaryHeatMap.png) |
| :-: |
| Figure 5: Top voew of the numerical solution $u(x,y) = \lambda^2 u = -8\cos{2x}\sin{2y}$ with the immersed boundary (red circle radius 1 centered at $(\pi, \pi)$), on a 128x128 grid. The solution matches the test solution. |


**Once your code works, experiment with this spacing. What happens when the IB points are very far apart compared to the finite difference mesh size? What happens when they are very close together? Look at both the accuracy (Problem 3 will give a better idea of this) and the speed of your code?**

We control the spacing through $\alpha$, where the distance between neighboring immersed boundary points is $\Delta s \approx \alpha \Delta x$ (so $N_{ib} = \frac{2 \pi}{\alpha \Delta x^2}$ ). I ran $\alpha \in \{0.25, 0.5, 1, 2, 4, 8\}$ at $N = 64$ for both Problems 2 and 3, and recorded run time, and accuracy. Since Problem 3 has no exact solution, we can measure accuracy by comparing against the boundary data. 

For Problem 2

| $\alpha$ | Error | Time (s) |
| --- | --- | --- |
| 0.25 | 5.33e-02 | 2.0540 |
| 0.5 | 5.33e-02 | 1.1108 |
| 1 | 5.33e-02 | 0.6534 |
| 2 | 5.33e-02 | 0.4132 |
| 4 | 5.41e-02 | 0.2443 |
| 8 | 5.76e-02 | 0.1207 |


For Problem 3

| $\alpha$ | Error | Time (s) |
| --- | --- | --- |
| 0.25 | 8.26e-10 | 12.50 |
| 0.5 | 1.42e-09 | 7.58 |
| 1 | 1.40e-09 | 5.11 |
| 2 | 9.62e-07 | 1.12 |
| 4 | 1.22e-02 | 0.27 |
| 8 | 2.50e-01 | 0.15 |

**IB points far apart ($\alpha > 1$):** When the IB points are very far apart compared to the finite difference mesh size the Guassian around the neighboring points don't overlap. The condition $u = boundary$ is enforced near each immersed boundary point, but in the gaps u is in plac, so it leaks through. Rather than acting like a wire the ciecle behaves like a few seperate points. Additionally, since there are fewer points the spreading operator and Jaconian operator must do fewer iterations and thus are less expensive. Thus, when the IB points are very far apart compared to the finite difference mesh size then the solution is inaccurate and runs fast. 

**IB points close togehter ($\alpha < 1$):** When the IB poitns are very close compared to the finite difference mesh size then there could be multiple points per grid cell (e.g. with alpha = 0.25 and N = 64 then there is about 4 points per grid cell) which makes the neighboring Guassians almost identoical, so the extra poitns add no new information and thus the accuracy of the solution does not improve. As we've discussed before the spreading operator and Jacobian operator iterate through all points and since there is more points this means the solver will take longer to run. 

## Problem 3: An almost minimal surface

### Setup

Problem 3 resuses the Problem 2 code with $f - 0$ and homogeneous Dirichlet data on the box. The only change is in the right-hand side, the top block is now zero, so teh circle data in the botton block is the only thing driving the solution. Unlike Problem 2, there is no exact solution to compare against.

### Results

| ![Descriptive alt text](figures/Problem3_PotentialPE.png) |
| :-: |
| Figure 6: Solution to Problem 3 on 128 x 128 grid, with $u=0$ on the boundaries of the box and $u = \cos{(2x)}\sin{(2y)}$ prescribed on te immersed boundary  (red cuve: circle of radius 1 centered at $(\pi, \pi)$ ). The surface is flat away from the circle and is pulled up and down to match the boundaries along it, forming an "almost minimal surface." |

Along the circle, the boundary data has two peaks and two dips reaching about $\pm 1$. Away from the circle, $u$ decays smoothly to zero at the walls, so the solution resembles a soap film stretched across a bent wire frame. Across the circle the sope of the surface changes abruptly. This is profuced by the IB charges $q$. 

## Problem 4: Charlie's rule

For Charlie's rule I decided to attempt to speed up the code for the spreading and Jacobian operators by setting $\delta_a$ to zero when $\delta_a$ is within zero with a prescribed tolerance. To test if change sped up the code without losing accuracy timed how long it took to run the solver for Problem 2 with the changes and without for $N = [16, 32, 64, 128]$. 

