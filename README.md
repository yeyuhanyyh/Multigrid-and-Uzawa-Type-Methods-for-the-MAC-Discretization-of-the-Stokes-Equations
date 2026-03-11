# Stokes Solvers on a MAC Grid

This repository contains three self-contained MATLAB programs for the final project in Numerical Linear Algebra.

Files
- `q1_multigrid_dgs.m`
- `q2_uzawa_cg.m`
- `q3_inexact_uzawa_pcg.m`

Each file is standalone:
- it builds the manufactured-solution test problem,
- runs the corresponding solver,
- prints a summary,
- and returns a MATLAB `struct` with the numerical results.

Model problem
We solve the steady Stokes equations on the unit square `(0,1)^2`:

`-Δu + ∇p = f`
`∇·u = 0`

with mixed boundary conditions:
- homogeneous Dirichlet conditions on the normal velocity boundaries,
- Neumann conditions on the tangential velocity boundaries.

The manufactured solution used in the code is

`u(x,y) = (1 - cos(2πx)) sin(2πy)`
`v(x,y) = -(1 - cos(2πy)) sin(2πx)`
`p(x,y) = x^3/3 - 1/12`

The forcing terms are generated from this exact solution, so the velocity error can be measured directly.

Grid layout
The code uses the standard MAC staggered grid:
- `u` is stored on vertical cell faces and has size `(N+1) x N`;
- `v` is stored on horizontal cell faces and has size `N x (N+1)`;
- `p` is stored at cell centers and has size `N x N`.

Boundary values are stored explicitly in the outer rows or columns of `u` and `v`.

What each file does

1. `q1_multigrid_dgs.m`
   Solves the full Stokes saddle-point system by a recursive V-cycle multigrid method with a distributive Gauss-Seidel (DGS) smoother.

   Default parameters:
   - `N = 1024`
   - `coarsestN = 4`
   - `nu1 = 2`
   - `nu2 = 2`
   - `tol = 1e-8`

2. `q2_uzawa_cg.m`
   Solves the Stokes system by the Uzawa iteration.
   At each outer step, the velocity subproblem is solved by matrix-free conjugate gradient (CG).

   Default parameters:
   - `N = 256`
   - `alpha = 1.0`
   - `cgTol = 1e-10`
   - `tol = 1e-8`

3. `q3_inexact_uzawa_pcg.m`
   Solves the Stokes system by the Inexact Uzawa iteration.
   At each outer step, the velocity subproblem is solved by preconditioned conjugate gradient (PCG), where the preconditioner is a multigrid V-cycle with symmetric Gauss-Seidel smoothing.

   Default parameters:
   - `N = 1024`
   - `alpha = 1.0`
   - `nu1 = 2`
   - `nu2 = 2`
   - `coarsestN = 2`
   - `pcgTol = 1e-3`
   - `tol = 1e-8`
   - `firstOuterPCGSteps = 2`

How to run
Open MATLAB in the repository folder and call any of the three functions:

`results1 = q1_multigrid_dgs;`
`results2 = q2_uzawa_cg;`
`results3 = q3_inexact_uzawa_pcg;`

You can also override the default parameters. For example:

`results1 = q1_multigrid_dgs(1024, 4, 2, 2, 1e-8, 50);`
`results2 = q2_uzawa_cg(256, 1.0, 1e-10, 1e-8, 50);`
`results3 = q3_inexact_uzawa_pcg(1024, 1.0, 2, 2, 2, 1e-3, 1e-8, 50, 2);`

Returned output
Each function returns a struct with fields such as:
- the grid size,
- the solver parameters,
- the total runtime,
- the number of outer iterations,
- the final residual,
- the final velocity error,
- the final arrays `u`, `v`, `p`,
- and the residual/error history.

Implementation notes
- All operators are applied in matrix-free form by finite-difference stencils.
- Pressure is normalized to have zero mean after each outer iteration. This removes the constant-pressure ambiguity.
- The multigrid codes assume that `N` is a power-of-two multiple of `coarsestN`.
- The code is written to stay close to the original project implementation, but the file structure and comments have been cleaned up and rewritten in English.

Suggested repository layout
If you want a very small and clean repository, keep only these files:

- `q1_multigrid_dgs.m`
- `q2_uzawa_cg.m`
- `q3_inexact_uzawa_pcg.m`
- `README.md`

That is enough to reproduce the three parts of the project without dozens of helper files.
