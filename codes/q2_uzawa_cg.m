function results = q2_uzawa_cg(N, alpha, cgTol, tol, maxOuter)
%Q2_UZAWA_CG Solve the MAC-discretized Stokes system by Uzawa + CG.
%
%   results = q2_uzawa_cg()
%   results = q2_uzawa_cg(N, alpha, cgTol, tol, maxOuter)
%
%   Default parameters:
%       N        = 256
%       alpha    = 1.0
%       cgTol    = 1e-10
%       tol      = 1e-8
%       maxOuter = 50
%
%   The code is self-contained and uses matrix-free stencil operations.

    if nargin < 1 || isempty(N),        N = 256;    end
    if nargin < 2 || isempty(alpha),    alpha = 1.0; end
    if nargin < 3 || isempty(cgTol),    cgTol = 1e-10; end
    if nargin < 4 || isempty(tol),      tol = 1e-8; end
    if nargin < 5 || isempty(maxOuter), maxOuter = 50; end

    [u, v, p, D, f, g, uExact, vExact] = setup_problem(N);

    [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
    r0 = residual_norm(ru, rv, rp, N);
    re = r0;
    er = velocity_error(u, v, uExact, vExact, N);

    residualHistory = zeros(maxOuter + 1, 1);
    errorHistory    = zeros(maxOuter + 1, 1);
    innerHistory    = zeros(maxOuter, 1);
    residualHistory(1) = re;
    errorHistory(1)    = er;

    elapsed = 0.0;
    outerIter = 0;

    while re > tol * r0 && outerIter < maxOuter
        tic;
        [u, v, cgIter] = cg_velocity_solve(u, v, p, f, g, N, cgTol);
        p = uzawa_pressure_update(u, v, p, N, alpha);
        p = p - mean(p(:));
        elapsed = elapsed + toc;

        [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
        re = residual_norm(ru, rv, rp, N);
        er = velocity_error(u, v, uExact, vExact, N);

        outerIter = outerIter + 1;
        innerHistory(outerIter)    = cgIter;
        residualHistory(outerIter + 1) = re;
        errorHistory(outerIter + 1)    = er;

        fprintf('Outer iteration %d, CG iterations = %d, residual = %.6e\n', ...
            outerIter, cgIter, re);
    end

    fprintf('\n');
    fprintf('Q2: Uzawa iteration with CG velocity solves\n');
    fprintf('N               = %d\n', N);
    fprintf('alpha           = %.6f\n', alpha);
    fprintf('CG tolerance    = %.2e\n', cgTol);
    fprintf('iterations      = %d\n', outerIter);
    fprintf('time (seconds)  = %.6f\n', elapsed);
    fprintf('final residual  = %.6e\n', re);
    fprintf('final error     = %.6e\n', er);

    results = struct();
    results.N = N;
    results.alpha = alpha;
    results.cgTol = cgTol;
    results.iterations = outerIter;
    results.time = elapsed;
    results.finalResidual = re;
    results.finalError = er;
    results.u = u;
    results.v = v;
    results.p = p;
    results.innerIterations = innerHistory(1:outerIter);
    results.residualHistory = residualHistory(1:outerIter+1);
    results.errorHistory = errorHistory(1:outerIter+1);
end

function [u, v, p, D, f, g, uExact, vExact] = setup_problem(N)
%SETUP_PROBLEM Build the manufactured-solution test problem on a MAC grid.

    h = 1 / N;

    u = zeros(N + 1, N);
    v = zeros(N, N + 1);
    p = zeros(N, N);
    D = zeros(N, N);

    f = zeros(N + 1, N);
    g = zeros(N, N + 1);

    uExact = zeros(N + 1, N);
    vExact = zeros(N, N + 1);

    bottomNeumann = zeros(N + 1, 1);
    topNeumann    = zeros(N + 1, 1);
    leftNeumann   = zeros(N + 1, 1);
    rightNeumann  = zeros(N + 1, 1);

    for i = 2:N
        x = (i - 1) * h;
        bottomNeumann(i) = -2 * pi * (1 - cos(2 * pi * x));
        topNeumann(i)    =  2 * pi * (1 - cos(2 * pi * x));
    end

    for j = 2:N
        y = (j - 1) * h;
        leftNeumann(j)  =  2 * pi * (1 - cos(2 * pi * y));
        rightNeumann(j) = -2 * pi * (1 - cos(2 * pi * y));
    end

    for i = 1:N+1
        for j = 1:N
            x = (i - 1) * h;
            y = (j - 0.5) * h;
            f(i, j) = -4 * pi^2 * (2 * cos(2 * pi * x) - 1) * sin(2 * pi * y) + x^2;
            uExact(i, j) = (1 - cos(2 * pi * x)) * sin(2 * pi * y);
        end
    end

    for i = 1:N
        for j = 1:N+1
            x = (i - 0.5) * h;
            y = (j - 1) * h;
            g(i, j) = 4 * pi^2 * (2 * cos(2 * pi * y) - 1) * sin(2 * pi * x);
            vExact(i, j) = -(1 - cos(2 * pi * y)) * sin(2 * pi * x);
        end
    end

    for i = 2:N
        f(i, 1) = f(i, 1) + bottomNeumann(i) / h;
        f(i, N) = f(i, N) + topNeumann(i) / h;
    end

    for j = 2:N
        g(1, j) = g(1, j) + leftNeumann(j) / h;
        g(N, j) = g(N, j) + rightNeumann(j) / h;
    end
end

function [u, v, cgIter] = cg_velocity_solve(u, v, p, f, g, N, cgTol)
%CG_VELOCITY_SOLVE Solve A*[u;v] = [f;g] - grad(p) by CG.

    maxInner = N * N;

    [rhsU, rhsV] = subtract_pressure_gradient(p, f, g, N);
    [ru, rv] = velocity_residual_with_pressure(u, v, p, f, g, N);

    bNormSq = velocity_dot(rhsU, rhsV, rhsU, rhsV);
    rNormSq = velocity_dot(ru, rv, ru, rv);

    cgIter = 0;
    if bNormSq == 0 || rNormSq == 0
        return;
    end

    pu = zeros(size(u));
    pv = zeros(size(v));
    rhoOld = 1.0;

    while cgIter < maxInner && rNormSq > (cgTol^2) * bNormSq
        cgIter = cgIter + 1;

        if cgIter == 1
            pu = ru;
            pv = rv;
        else
            beta = rNormSq / rhoOld;
            pu = ru + beta * pu;
            pv = rv + beta * pv;
        end

        [wu, wv] = apply_velocity_operator(pu, pv, N);
        alpha = rNormSq / velocity_dot(pu, pv, wu, wv);

        u = u + alpha * pu;
        v = v + alpha * pv;

        ru = ru - alpha * wu;
        rv = rv - alpha * wv;

        rhoOld = rNormSq;
        rNormSq = velocity_dot(ru, rv, ru, rv);
    end
end

function p = uzawa_pressure_update(u, v, p, N, alpha)
%UZAWA_PRESSURE_UPDATE Standard pressure update in Uzawa iteration.

    h = 1 / N;
    for i = 1:N
        for j = 1:N
            div = (u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j)) / h;
            p(i, j) = p(i, j) - alpha * div;
        end
    end
end

function [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N)
%STOKES_RESIDUAL Full residual of the Stokes system.

    [ru, rv] = velocity_residual_with_pressure(u, v, p, f, g, N);

    h = 1 / N;
    rp = zeros(N, N);
    for i = 1:N
        for j = 1:N
            rp(i, j) = D(i, j) - (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h;
        end
    end
end

function [ru, rv] = velocity_residual_with_pressure(u, v, p, f, g, N)
%VELOCITY_RESIDUAL_WITH_PRESSURE Residual of A*[u;v] + grad(p) = [f;g].

    [rhsU, rhsV] = subtract_pressure_gradient(p, f, g, N);
    [Au, Av] = apply_velocity_operator(u, v, N);

    ru = rhsU - Au;
    rv = rhsV - Av;
end

function [shiftedF, shiftedG] = subtract_pressure_gradient(p, f, g, N)
%SUBTRACT_PRESSURE_GRADIENT Compute [f;g] - grad(p) on the MAC grid.

    h = 1 / N;

    shiftedF = zeros(N + 1, N);
    shiftedG = zeros(N, N + 1);

    shiftedF(1, :)   = 0.0;
    shiftedF(N + 1, :) = 0.0;
    shiftedG(:, 1)   = 0.0;
    shiftedG(:, N + 1) = 0.0;

    for i = 2:N
        for j = 1:N
            shiftedF(i, j) = f(i, j) - (p(i, j) - p(i - 1, j)) / h;
        end
    end

    for j = 2:N
        for i = 1:N
            shiftedG(i, j) = g(i, j) - (p(i, j) - p(i, j - 1)) / h;
        end
    end
end

function [Au, Av] = apply_velocity_operator(u, v, N)
%APPLY_VELOCITY_OPERATOR Apply the discrete velocity Laplacian block A.

    h = 1 / N;
    Au = zeros(N + 1, N);
    Av = zeros(N, N + 1);

    Au(1, :)   = 0.0;
    Au(N + 1, :) = 0.0;

    for i = 2:N
        j = 1;
        Au(i, j) = (3 * u(i, j) - u(i - 1, j) - u(i + 1, j) - u(i, j + 1)) / (h * h);
        for j = 2:N-1
            Au(i, j) = (4 * u(i, j) - u(i - 1, j) - u(i + 1, j) - u(i, j - 1) - u(i, j + 1)) / (h * h);
        end
        j = N;
        Au(i, j) = (3 * u(i, j) - u(i - 1, j) - u(i + 1, j) - u(i, j - 1)) / (h * h);
    end

    Av(:, 1)   = 0.0;
    Av(:, N + 1) = 0.0;

    for j = 2:N
        i = 1;
        Av(i, j) = (3 * v(i, j) - v(i + 1, j) - v(i, j - 1) - v(i, j + 1)) / (h * h);
        for i = 2:N-1
            Av(i, j) = (4 * v(i, j) - v(i - 1, j) - v(i + 1, j) - v(i, j - 1) - v(i, j + 1)) / (h * h);
        end
        i = N;
        Av(i, j) = (3 * v(i, j) - v(i - 1, j) - v(i, j - 1) - v(i, j + 1)) / (h * h);
    end
end

function value = velocity_dot(u1, v1, u2, v2)
%VELOCITY_DOT Euclidean inner product on the velocity unknowns.

    value = sum(u1(:) .* u2(:)) + sum(v1(:) .* v2(:));
end

function re = residual_norm(ru, rv, rp, N)
%RESIDUAL_NORM Euclidean norm of the full residual, scaled by 1/N.

    re = sqrt(sum(ru(:).^2) + sum(rv(:).^2) + sum(rp(:).^2)) / N;
end

function er = velocity_error(u, v, uExact, vExact, N)
%VELOCITY_ERROR Euclidean velocity error, scaled by 1/N.

    er = sqrt(sum((u(:) - uExact(:)).^2) + sum((v(:) - vExact(:)).^2)) / N;
end