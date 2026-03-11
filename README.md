# Multigrid and Uzawa-Type Methods for the MAC Discretization of the Stokes Equations

This repository contains three self-contained MATLAB programs for the final project in Numerical Linear Algebra.

The code solves the MAC (marker-and-cell) discretization of the steady Stokes equations on a staggered grid by three different iterative frameworks:

1. V-cycle multigrid with the distributive Gauss-Seidel (DGS) smoother;
2. Uzawa iteration with a conjugate-gradient (CG) solve for the velocity subproblem;
3. Inexact Uzawa iteration with a multigrid-preconditioned conjugate-gradient (PCG) solve for the velocity subproblem.

The repository is intentionally kept small: each problem is implemented in a single MATLAB file.

Files
- `q1_multigrid_dgs.m`
- `q2_uzawa_cg.m`
- `q3_inexact_uzawa_pcg.m`

--------------------------------------------------
1. Model problem
--------------------------------------------------

We consider the steady Stokes equations on the unit square

    Omega = (0,1)^2

with velocity `u = (u,v)^T`, pressure `p`, and forcing `f = (f,g)^T`:

    -Delta u + grad p = f    in Omega
    div u = 0                in Omega

The boundary conditions are mixed Dirichlet-Neumann conditions:

For the horizontal velocity component:
    du/dn = b   on y = 0
    du/dn = t   on y = 1
    u = 0       on x = 0, 1

For the vertical velocity component:
    dv/dn = l   on x = 0
    dv/dn = r   on x = 1
    v = 0       on y = 0, 1

To verify the implementation, the code uses the manufactured exact solution

    u(x,y) = (1 - cos(2*pi*x)) * sin(2*pi*y)
    v(x,y) = -(1 - cos(2*pi*y)) * sin(2*pi*x)
    p(x,y) = x^3 / 3 - 1/12

The forcing terms are obtained by substituting this exact solution into the Stokes equations:

    f1(x,y) = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y) + x^2
    f2(x,y) =  4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x)

This allows the velocity error to be measured directly.

--------------------------------------------------
2. MAC staggered-grid discretization
--------------------------------------------------

The code uses the standard MAC staggered grid.

For a mesh size `h = 1/N`:

- the horizontal velocity `u` is stored on vertical cell faces,
- the vertical velocity `v` is stored on horizontal cell faces,
- the pressure `p` is stored at cell centers.

More precisely:

    u_{i,j-1/2} ~ u(x_i, y_{j-1/2})
    v_{i-1/2,j} ~ v(x_{i-1/2}, y_j)
    p_{i-1/2,j-1/2} ~ p(x_{i-1/2}, y_{j-1/2})

In the MATLAB implementation:

- `u` has size `(N+1) x N`
- `v` has size `N x (N+1)`
- `p` has size `N x N`

The discrete divergence at a cell center is

    div_h(u,v) =
        ( u(i+1,j) - u(i,j) ) / h
      + ( v(i,j+1) - v(i,j) ) / h

The discrete velocity error reported by the code is

    e_N = h * sqrt( sum |u_h - u_exact|^2 + sum |v_h - v_exact|^2 )

where the sums run over all staggered-grid velocity unknowns.

The discrete Stokes system has the saddle-point form

    [ A   B ] [ U ] = [ F ]
    [ B^T 0 ] [ P ]   [ 0 ]

where:
- `U` is the vector of all velocity unknowns,
- `P` is the vector of all pressure unknowns,
- `A` is the discrete velocity Laplacian block,
- `B` is the discrete pressure-gradient operator.

Because pressure is determined only up to an additive constant, the code normalizes pressure after each outer iteration by subtracting its mean.

--------------------------------------------------
3. What each MATLAB file does
--------------------------------------------------

3.1 `q1_multigrid_dgs.m`

This file solves the full MAC Stokes system by a V-cycle multigrid method with the distributive Gauss-Seidel (DGS) smoother.

Main ingredients:
- matrix-free stencil application;
- DGS momentum relaxation;
- local divergence correction;
- restriction and prolongation between grids;
- recursive V-cycle down to the coarsest grid.

The DGS smoother consists of two parts:

(1) momentum relaxation:
    approximately solve the velocity equations with pressure fixed;

(2) divergence correction:
    locally modify neighboring face velocities so that the discrete divergence defect is reduced or removed.

Default parameters:
- `N = 1024`
- `coarsestN = 4`
- `nu1 = 2`
- `nu2 = 2`
- `tol = 1e-8`
- `maxCycles = 50`

How to run:
    results = q1_multigrid_dgs;

or with explicit parameters:
    results = q1_multigrid_dgs(1024, 4, 2, 2, 1e-8, 50);

Expected behavior:
- a small, almost grid-independent number of outer V-cycles;
- second-order velocity accuracy;
- much better efficiency than a basic smoother such as Jacobi.

--------------------------------------------------
3.2 `q2_uzawa_cg.m`

This file solves the Stokes system by the Uzawa iteration.

The outer iteration is

    A U^{k+1} = F - B P^k
    P^{k+1}   = P^k - alpha * div_h(U^{k+1})

In the code, the velocity subproblem is solved by matrix-free conjugate gradient (CG).

Default parameters:
- `N = 256`
- `alpha = 1.0`
- `cgTol = 1e-10`
- `tol = 1e-8`
- `maxOuter = 50`

How to run:
    results = q2_uzawa_cg;

or with explicit parameters:
    results = q2_uzawa_cg(256, 1.0, 1e-10, 1e-8, 50);

Expected behavior:
- very few outer Uzawa iterations;
- but potentially many inner CG iterations;
- therefore the method can be accurate but still expensive if the inner tolerance is too strict.

Important note:
This single-file version is written to be clear and self-contained.
It keeps the main Uzawa + CG structure of the original project, but the inner stopping rule is written in a standard residual-based form.

--------------------------------------------------
3.3 `q3_inexact_uzawa_pcg.m`

This file solves the Stokes system by the Inexact Uzawa iteration.

The outer iteration again updates pressure through the discrete divergence, but now the velocity subproblem is solved only approximately by PCG.

At each outer step:
- form the shifted velocity right-hand side `F - B P^k`,
- solve the velocity equation approximately by PCG,
- use a multigrid V-cycle as the preconditioner,
- update pressure.

The multigrid preconditioner uses symmetric Gauss-Seidel smoothing, which is appropriate for PCG.

Default parameters:
- `N = 1024`
- `alpha = 1.0`
- `nu1 = 2`
- `nu2 = 2`
- `coarsestN = 2`
- `pcgTol = 1e-3`
- `tol = 1e-8`
- `maxOuter = 50`
- `firstOuterPCGSteps = 2`

How to run:
    results = q3_inexact_uzawa_pcg;

or with explicit parameters:
    results = q3_inexact_uzawa_pcg(1024, 1.0, 2, 2, 2, 1e-3, 1e-8, 50, 2);

Expected behavior:
- very few outer iterations;
- only a few PCG steps per outer iteration;
- overall efficiency comparable to the multigrid solver in Problem 1.

Important note:
This single-file version is a cleaned, standard inexact-Uzawa implementation.
It is designed for readability and reproducibility, and its inner stopping rule is simpler than the original multi-file project code.

--------------------------------------------------
4. Returned output
--------------------------------------------------

Each function returns a MATLAB `struct` containing:
- the problem size,
- the solver parameters,
- the number of outer iterations,
- the total runtime,
- the final residual,
- the final velocity error,
- the final arrays `u`, `v`, and `p`,
- the residual history,
- the error history,
- and, for Q2 and Q3, the inner iteration counts.

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

For Q2 and Q3 there is also:
- `innerIterations`

--------------------------------------------------
5. Interpreting the numerical results
--------------------------------------------------

For Problem 1 (multigrid with DGS):
- the most important indicators are the number of V-cycles and the final velocity error;
- the number of cycles should remain small as `N` increases;
- this is the signature of an efficient multigrid method.

For Problem 2 (Uzawa + CG):
- the number of outer iterations may be very small;
- however, the total runtime can still be large because the inner CG solve may be expensive.

For Problem 3 (Inexact Uzawa + PCG):
- the goal is to reduce the cost of the velocity solve while keeping the outer Uzawa iteration effective;
- in a good regime, only a few PCG iterations are needed per outer step.

When comparing timings, keep in mind:
- runtime depends strongly on hardware,
- runtime also depends on MATLAB version,
- iteration counts and error levels are usually more stable than wall-clock time.

--------------------------------------------------
6. Notes on implementation
--------------------------------------------------

- All linear operators are applied in matrix-free form.
- No large sparse matrix is assembled explicitly.
- Pressure is normalized to zero mean after each outer iteration.
- The multigrid codes assume that `N` is a power-of-two multiple of `coarsestN`.
- The implementation is intentionally written in a direct and readable style rather than aggressively vectorized MATLAB style.

--------------------------------------------------
7. Minimal repository layout
--------------------------------------------------

A minimal clean repository can contain only:

- `q1_multigrid_dgs.m`
- `q2_uzawa_cg.m`
- `q3_inexact_uzawa_pcg.m`
- `README.md`

This is enough to reproduce the three parts of the project in a compact and understandable form.

--------------------------------------------------
8. Suggested citation / project description
--------------------------------------------------

Suggested one-sentence project summary:

Implemented matrix-free MATLAB solvers for staggered-grid Stokes equations, including V-cycle multigrid with distributive Gauss-Seidel smoothing, Uzawa iteration with conjugate gradient velocity solves, and Inexact Uzawa iteration with multigrid-preconditioned conjugate gradient.
