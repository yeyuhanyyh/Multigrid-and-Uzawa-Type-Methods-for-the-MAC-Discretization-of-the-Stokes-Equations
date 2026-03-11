# Multigrid and Uzawa-Type Methods for the MAC Discretization of the Stokes Equations

This repository contains MATLAB implementations of three iterative solvers for the MAC (marker-and-cell) discretization of the steady Stokes equations on a staggered grid.

The project studies three related questions:

1. How effectively does a $V$-cycle multigrid method with the distributive Gauss-Seidel (DGS) smoother solve the full saddle-point system?
2. How does the Uzawa iteration perform when the velocity subproblem is solved by conjugate gradient (CG)?
3. How much can the Inexact Uzawa method be accelerated by using a multigrid-preconditioned conjugate gradient (PCG) solver for the velocity block?

This repository accompanies a final project for **Numerical Linear Algebra**.

* * *

## Repository Structure

```text
.
├── q1_multigrid_dgs.m
├── q2_uzawa_cg.m
├── q3_inexact_uzawa_pcg.m
└── README.md
```

Each MATLAB file is self-contained and corresponds to one part of the project:

- `q1_multigrid_dgs.m`: Problem 1, $V$-cycle multigrid with DGS smoothing;
- `q2_uzawa_cg.m`: Problem 2, Uzawa iteration with CG solves of the velocity subproblem;
- `q3_inexact_uzawa_pcg.m`: Problem 3, Inexact Uzawa iteration with a multigrid-preconditioned CG solve.

* * *

## Continuous Model

We solve the steady Stokes equations on the unit square $\Omega = (0,1)^2$:

```math
-\Delta \mathbf{u} + \nabla p = \mathbf{f}, \qquad \nabla \cdot \mathbf{u} = 0 \quad \text{in } \Omega,
```

where $\mathbf{u} = (u,v)^\top$ is the velocity and $p$ is the pressure.

The boundary conditions are mixed Dirichlet-Neumann conditions:

For the horizontal velocity component $u$,

```math
\frac{\partial u}{\partial n} = b \quad \text{on } y=0, \qquad
\frac{\partial u}{\partial n} = t \quad \text{on } y=1, \qquad
u = 0 \quad \text{on } x=0,1.
```

For the vertical velocity component $v$,

```math
\frac{\partial v}{\partial n} = l \quad \text{on } x=0, \qquad
\frac{\partial v}{\partial n} = r \quad \text{on } x=1, \qquad
v = 0 \quad \text{on } y=0,1.
```

* * *

## Manufactured Solution

To verify the implementation, the code uses the exact solution

```math
u(x,y) = \bigl(1-\cos(2\pi x)\bigr)\sin(2\pi y),
```

```math
v(x,y) = -\bigl(1-\cos(2\pi y)\bigr)\sin(2\pi x),
```

```math
p(x,y) = \frac{x^3}{3} - \frac{1}{12}.
```

Substituting this exact solution into the Stokes equations gives the forcing terms

```math
f_1(x,y) = -4\pi^2\bigl(2\cos(2\pi x)-1\bigr)\sin(2\pi y) + x^2,
```

```math
f_2(x,y) = 4\pi^2\bigl(2\cos(2\pi y)-1\bigr)\sin(2\pi x).
```

The corresponding Neumann boundary data are generated automatically in the MATLAB codes.

* * *

## MAC Discretization

We use the standard MAC staggered grid with mesh size $h = 1/N$.

- The horizontal velocity is stored on vertical cell faces:
  $u_{i,j-\frac12} \approx u(x_i, y_{j-\frac12})$.
- The vertical velocity is stored on horizontal cell faces:
  $v_{i-\frac12,j} \approx v(x_{i-\frac12}, y_j)$.
- The pressure is stored at cell centers:
  $p_{i-\frac12,j-\frac12} \approx p(x_{i-\frac12}, y_{j-\frac12})$.

In the MATLAB implementation:

- `u` has size `(N+1) x N`,
- `v` has size `N x (N+1)`,
- `p` has size `N x N`.

The discrete divergence at a cell center is

```math
(\nabla_h \cdot \mathbf{u})_{i-\frac12,j-\frac12}
=
\frac{u_{i+1,j}-u_{i,j}}{h}
+
\frac{v_{i,j+1}-v_{i,j}}{h}.
```

The discrete Stokes system has the saddle-point form

```math
\begin{pmatrix}
A & B \\
B^\top & 0
\end{pmatrix}
\begin{pmatrix}
U \\
P
\end{pmatrix}
=
\begin{pmatrix}
F \\
0
\end{pmatrix},
```

where:

- $U$ is the vector of all velocity unknowns,
- $P$ is the vector of all pressure unknowns,
- $A$ is the discrete velocity Laplacian block,
- $B$ is the discrete gradient operator.

Since pressure is only determined up to an additive constant, all three codes normalize the pressure by subtracting its mean after each outer iteration.

* * *

## Error Metric

The velocity error reported by the codes is the discrete $L^2$-type error

```math
e_N
=
h\left(
\sum_{j=1}^{N}\sum_{i=1}^{N-1}
\left|u_{i,j-\frac12}-u(x_i,y_{j-\frac12})\right|^2
+
\sum_{j=1}^{N-1}\sum_{i=1}^{N}
\left|v_{i-\frac12,j}-v(x_{i-\frac12},y_j)\right|^2
\right)^{1/2}.
```

This is the quantity used in the project report to verify second-order accuracy.

* * *

## Algorithms

### 1. V-cycle Multigrid with DGS Smoothing

The first solver applies a recursive $V$-cycle multigrid method to the full Stokes system.

The smoother is the **distributive Gauss-Seidel (DGS)** method, which consists of two parts:

1. a Gauss-Seidel relaxation step for the momentum equations;
2. a local divergence-correction step to reduce the discrete incompressibility defect.

At a high level, one $V$-cycle consists of:

- $\nu_1$ pre-smoothing steps,
- residual restriction to the next coarser grid,
- recursive coarse-grid correction,
- prolongation of the correction,
- $\nu_2$ post-smoothing steps.

The code is matrix-free: all operator applications are performed by finite-difference stencils rather than by assembling large sparse matrices.

### 2. Uzawa Iteration with CG

The second solver uses the classical Uzawa iteration:

```math
A U^{k+1} = F - B P^k,
```

```math
P^{k+1} = P^k - \alpha \, \nabla_h \cdot U^{k+1}.
```

At each outer iteration, the velocity subproblem is solved by matrix-free conjugate gradient (CG).

This method typically requires very few outer iterations, but the inner CG solve can become expensive if a very strict tolerance is used.

### 3. Inexact Uzawa with Multigrid-Preconditioned CG

The third solver replaces the exact or highly accurate velocity solve in Uzawa iteration by a preconditioned conjugate gradient method.

At each outer iteration:

1. form the shifted velocity right-hand side $F - B P^k$;
2. approximately solve the velocity block by PCG;
3. use a multigrid $V$-cycle as the preconditioner;
4. update pressure by the discrete divergence.

The multigrid preconditioner uses **symmetric Gauss-Seidel** smoothing, which is suitable for PCG.

This method is designed to keep the number of outer Uzawa steps small while greatly reducing the cost of the velocity solve.

* * *

## Numerical Experiments

### Problem 1: Multigrid Baseline

`q1_multigrid_dgs.m` uses the default parameters

- `N = 1024`
- `coarsestN = 4`
- `nu1 = 2`
- `nu2 = 2`
- `tol = 1e-8`

The main quantities of interest are:

- the number of outer $V$-cycles,
- the total runtime,
- the final residual,
- the velocity error.

This experiment is meant to demonstrate:

- second-order accuracy of the MAC discretization,
- near grid-independent convergence of multigrid,
- the effectiveness of DGS smoothing.

### Problem 2: Uzawa + CG

`q2_uzawa_cg.m` uses the default parameters

- `N = 256`
- `alpha = 1.0`
- `cgTol = 1e-10`
- `tol = 1e-8`

The key outputs are:

- the number of outer Uzawa iterations,
- the number of inner CG iterations per outer step,
- the final residual,
- the velocity error.

This experiment shows that Uzawa may converge in very few outer steps, but the total computational cost depends strongly on the inner CG solve.

### Problem 3: Inexact Uzawa + PCG

`q3_inexact_uzawa_pcg.m` uses the default parameters

- `N = 1024`
- `alpha = 1.0`
- `nu1 = 2`
- `nu2 = 2`
- `coarsestN = 2`
- `pcgTol = 1e-3`
- `tol = 1e-8`
- `firstOuterPCGSteps = 2`

The key outputs are:

- the number of outer Inexact Uzawa iterations,
- the number of inner PCG iterations,
- the final residual,
- the velocity error.

This experiment demonstrates that a multigrid preconditioner can make the velocity solve much cheaper while preserving the efficiency of the outer iteration.

* * *

## How to Run

Open MATLAB in the repository folder and run any of the following commands.

### Problem 1

```matlab
results1 = q1_multigrid_dgs;
```

or with explicit parameters,

```matlab
results1 = q1_multigrid_dgs(1024, 4, 2, 2, 1e-8, 50);
```

### Problem 2

```matlab
results2 = q2_uzawa_cg;
```

or with explicit parameters,

```matlab
results2 = q2_uzawa_cg(256, 1.0, 1e-10, 1e-8, 50);
```

### Problem 3

```matlab
results3 = q3_inexact_uzawa_pcg;
```

or with explicit parameters,

```matlab
results3 = q3_inexact_uzawa_pcg(1024, 1.0, 2, 2, 2, 1e-3, 1e-8, 50, 2);
```

* * *

## Returned Output

Each function returns a MATLAB `struct` containing:

- the problem size,
- the solver parameters,
- the number of outer iterations,
- the total runtime,
- the final residual,
- the final velocity error,
- the final arrays `u`, `v`, and `p`,
- the residual history,
- the error history.

For Problems 2 and 3, the returned struct also contains the inner iteration counts.

Typical fields include:

- `N`
- `iterations`
- `time`
- `finalResidual`
- `finalError`
- `u`
- `v`
- `p`
- `residualHistory`
- `errorHistory`

and, where appropriate,

- `innerIterations`

* * *

## Notes

- All three MATLAB programs are written in matrix-free form.
- No large sparse matrix is assembled explicitly.
- Pressure is normalized to zero mean after each outer iteration.
- The multigrid solvers assume that `N` is a power-of-two multiple of `coarsestN`.
- The code is written to be compact, readable, and easy to run directly from MATLAB.

* * *

## Project Summary

A concise summary of this repository is:

Implemented matrix-free MATLAB solvers for staggered-grid Stokes equations, including $V$-cycle multigrid with distributive Gauss-Seidel smoothing, Uzawa iteration with conjugate gradient velocity solves, and Inexact Uzawa iteration with multigrid-preconditioned conjugate gradient.
