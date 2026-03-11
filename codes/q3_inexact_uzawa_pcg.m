function results = q3_inexact_uzawa_pcg(N, alpha, nu1, nu2, coarsestN, pcgTol, tol, maxOuter, firstOuterPCGSteps)
%Q3_INEXACT_UZAWA_PCG Solve the MAC-discretized Stokes system by
%inexact Uzawa iteration with a multigrid-preconditioned CG inner solve.
%
%   results = q3_inexact_uzawa_pcg()
%   results = q3_inexact_uzawa_pcg(N, alpha, nu1, nu2, coarsestN, ...
%                                  pcgTol, tol, maxOuter, firstOuterPCGSteps)
%
%   Default parameters:
%       N                 = 1024
%       alpha             = 1.0
%       nu1               = 2
%       nu2               = 2
%       coarsestN         = 2
%       pcgTol            = 1e-3
%       tol               = 1e-8
%       maxOuter          = 50
%       firstOuterPCGSteps= 2
%
%   The code is self-contained and uses matrix-free stencil operations.

    if nargin < 1 || isempty(N),                 N = 1024;    end
    if nargin < 2 || isempty(alpha),             alpha = 1.0; end
    if nargin < 3 || isempty(nu1),               nu1 = 2;     end
    if nargin < 4 || isempty(nu2),               nu2 = 2;     end
    if nargin < 5 || isempty(coarsestN),         coarsestN = 2; end
    if nargin < 6 || isempty(pcgTol),            pcgTol = 1e-3; end
    if nargin < 7 || isempty(tol),               tol = 1e-8;  end
    if nargin < 8 || isempty(maxOuter),          maxOuter = 50; end
    if nargin < 9 || isempty(firstOuterPCGSteps), firstOuterPCGSteps = 2; end

    ratio = N / coarsestN;
    if N < coarsestN || mod(N, coarsestN) ~= 0 || abs(log2(ratio) - round(log2(ratio))) > 1e-12
        error('N must be a power-of-two multiple of coarsestN.');
    end

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
        [shiftedF, shiftedG] = subtract_pressure_gradient(p, f, g, N);

        if outerIter == 0
            fixedSteps = firstOuterPCGSteps;
        else
            fixedSteps = [];
        end

        [u, v, pcgIter] = pcg_velocity_solve( ...
            u, v, shiftedF, shiftedG, N, nu1, nu2, coarsestN, pcgTol, fixedSteps);

        p = uzawa_pressure_update(u, v, p, N, alpha);
        p = p - mean(p(:));

        elapsed = elapsed + toc;

        [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
        re = residual_norm(ru, rv, rp, N);
        er = velocity_error(u, v, uExact, vExact, N);

        outerIter = outerIter + 1;
        innerHistory(outerIter) = pcgIter;
        residualHistory(outerIter + 1) = re;
        errorHistory(outerIter + 1)    = er;

        fprintf('Outer iteration %d, PCG iterations = %d, residual = %.6e\n', ...
            outerIter, pcgIter, re);
    end

    fprintf('\n');
    fprintf('Q3: Inexact Uzawa with multigrid-preconditioned CG\n');
    fprintf('N                   = %d\n', N);
    fprintf('alpha               = %.6f\n', alpha);
    fprintf('nu1, nu2            = %d, %d\n', nu1, nu2);
    fprintf('coarsestN           = %d\n', coarsestN);
    fprintf('PCG tolerance       = %.2e\n', pcgTol);
    fprintf('first outer PCG it. = %d\n', firstOuterPCGSteps);
    fprintf('iterations          = %d\n', outerIter);
    fprintf('time (seconds)      = %.6f\n', elapsed);
    fprintf('final residual      = %.6e\n', re);
    fprintf('final error         = %.6e\n', er);

    results = struct();
    results.N = N;
    results.alpha = alpha;
    results.nu1 = nu1;
    results.nu2 = nu2;
    results.coarsestN = coarsestN;
    results.pcgTol = pcgTol;
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

function [u, v, pcgIter] = pcg_velocity_solve(u, v, f, g, N, nu1, nu2, coarsestN, pcgTol, fixedSteps)
%PCG_VELOCITY_SOLVE Solve A*[u;v] = [f;g] by PCG with a V-cycle preconditioner.

    maxInner = N * N;

    [ru, rv] = velocity_residual_without_pressure(u, v, f, g, N);
    initialResidualSq = velocity_dot(ru, rv, ru, rv);

    pcgIter = 0;
    if initialResidualSq == 0
        return;
    end

    pu = zeros(size(u));
    pv = zeros(size(v));
    rhoOld = 1.0;

    while pcgIter < maxInner
        [zu, zv] = velocity_preconditioner(zeros(N + 1, N), zeros(N, N + 1), ...
                                           ru, rv, nu1, nu2, coarsestN, N);
        rho = velocity_dot(ru, rv, zu, zv);

        if pcgIter == 0
            pu = zu;
            pv = zv;
        else
            beta = rho / rhoOld;
            pu = zu + beta * pu;
            pv = zv + beta * pv;
        end

        [wu, wv] = apply_velocity_operator(pu, pv, N);
        alpha = rho / velocity_dot(pu, pv, wu, wv);

        u = u + alpha * pu;
        v = v + alpha * pv;

        ru = ru - alpha * wu;
        rv = rv - alpha * wv;

        pcgIter = pcgIter + 1;
        rhoOld = rho;

        currentResidualSq = velocity_dot(ru, rv, ru, rv);

        if ~isempty(fixedSteps)
            if pcgIter >= fixedSteps
                break;
            end
        else
            if currentResidualSq <= (pcgTol^2) * initialResidualSq
                break;
            end
        end
    end
end

function [u, v] = velocity_preconditioner(u, v, f, g, nu1, nu2, coarsestN, N)
%VELOCITY_PRECONDITIONER Symmetric multigrid V-cycle for the velocity block.

    if N == coarsestN
        for k = 1:100
            [u, v] = symmetric_gauss_seidel_velocity(u, v, f, g, N);
        end
        return;
    end

    for k = 1:nu1
        [u, v] = symmetric_gauss_seidel_velocity(u, v, f, g, N);
    end

    [ru, rv] = velocity_residual_without_pressure(u, v, f, g, N);
    [ruc, rvc] = restrict_velocity_residual(ru, rv, N);

    [uc, vc] = initialize_velocity_fields(N / 2);
    [uc, vc] = velocity_preconditioner(uc, vc, ruc, rvc, nu1, nu2, coarsestN, N / 2);

    [eu, ev] = prolong_velocity_correction(uc, vc, N);
    u = u + eu;
    v = v + ev;

    for k = 1:nu2
        [u, v] = symmetric_gauss_seidel_velocity(u, v, f, g, N);
    end
end

function [u, v] = initialize_velocity_fields(N)
%INITIALIZE_VELOCITY_FIELDS Zero initial guess for the velocity block.

    u = zeros(N + 1, N);
    v = zeros(N, N + 1);
end

function [u, v] = symmetric_gauss_seidel_velocity(u, v, f, g, N)
%SYMMETRIC_GAUSS_SEIDEL_VELOCITY Symmetric GS smoother for the velocity block.

    h = 1 / N;

    % Forward sweep for u
    u(1, :) = 0.0;
    for i = 2:N
        j = 1;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j + 1)) / 3;
        for j = 2:N-1
            u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1)) / 4;
        end
        j = N;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1)) / 3;
    end
    u(N + 1, :) = 0.0;

    % Forward sweep for v
    v(:, 1) = 0.0;
    for j = 2:N
        i = 1;
        v(i, j) = (h * h * g(i, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1)) / 3;
        for i = 2:N-1
            v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1)) / 4;
        end
        i = N;
        v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i, j - 1) + v(i, j + 1)) / 3;
    end
    v(:, N + 1) = 0.0;

    % Backward sweep for v
    v(:, N + 1) = 0.0;
    for j = N:-1:2
        i = N;
        v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i, j - 1) + v(i, j + 1)) / 3;
        for i = N-1:-1:2
            v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1)) / 4;
        end
        i = 1;
        v(i, j) = (h * h * g(i, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1)) / 3;
    end
    v(:, 1) = 0.0;

    % Backward sweep for u
    u(N + 1, :) = 0.0;
    for i = N:-1:2
        j = N;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1)) / 3;
        for j = N-1:-1:2
            u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1)) / 4;
        end
        j = 1;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j + 1)) / 3;
    end
    u(1, :) = 0.0;
end

function [ruc, rvc] = restrict_velocity_residual(ru, rv, N)
%RESTRICT_VELOCITY_RESIDUAL Restrict a velocity residual to the next coarser grid.

    Nc = N / 2;
    ruc = zeros(Nc + 1, Nc);
    rvc = zeros(Nc, Nc + 1);

    ruc(1, :)   = 0.0;
    ruc(Nc + 1, :) = 0.0;
    for i = 2:Nc
        for j = 1:Nc
            ruc(i, j) = (ru(2*i,   2*j)   + ru(2*i,   2*j-1) + ...
                         ru(2*i-2, 2*j)   + ru(2*i-2, 2*j-1)) / 8 + ...
                        (ru(2*i-1, 2*j)   + ru(2*i-1, 2*j-1)) / 4;
        end
    end

    rvc(:, 1)   = 0.0;
    rvc(:, Nc + 1) = 0.0;
    for j = 2:Nc
        for i = 1:Nc
            rvc(i, j) = (rv(2*i,   2*j)   + rv(2*i-1, 2*j)   + ...
                         rv(2*i,   2*j-2) + rv(2*i-1, 2*j-2)) / 8 + ...
                        (rv(2*i,   2*j-1) + rv(2*i-1, 2*j-1)) / 4;
        end
    end
end

function [eu, ev] = prolong_velocity_correction(euc, evc, N)
%PROLONG_VELOCITY_CORRECTION Prolong a coarse-grid velocity correction.

    Nc = N / 2;
    eu = zeros(N + 1, N);
    ev = zeros(N, N + 1);

    for j = 1:Nc
        for i = 1:Nc
            eu(2*i-1, 2*j-1) = euc(i, j);
            eu(2*i-1, 2*j)   = euc(i, j);
            eu(2*i,   2*j-1) = 0.5 * (euc(i, j) + euc(i + 1, j));
            eu(2*i,   2*j)   = 0.5 * (euc(i, j) + euc(i + 1, j));
        end
    end
    eu(1, :)   = 0.0;
    eu(N + 1, :) = 0.0;

    for i = 1:Nc
        for j = 1:Nc
            ev(2*i-1, 2*j-1) = evc(i, j);
            ev(2*i,   2*j-1) = evc(i, j);
            ev(2*i-1, 2*j)   = 0.5 * (evc(i, j) + evc(i, j + 1));
            ev(2*i,   2*j)   = 0.5 * (evc(i, j) + evc(i, j + 1));
        end
    end
    ev(:, 1)   = 0.0;
    ev(:, N + 1) = 0.0;
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

function [ru, rv] = velocity_residual_without_pressure(u, v, f, g, N)
%VELOCITY_RESIDUAL_WITHOUT_PRESSURE Residual of A*[u;v] = [f;g].

    [Au, Av] = apply_velocity_operator(u, v, N);
    ru = f - Au;
    rv = g - Av;
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