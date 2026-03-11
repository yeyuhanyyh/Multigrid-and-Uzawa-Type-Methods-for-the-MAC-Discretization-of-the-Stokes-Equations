function results = q1_multigrid_dgs(N, coarsestN, nu1, nu2, tol, maxCycles)
%Q1_MULTIGRID_DGS Solve the MAC-discretized Stokes system by a DGS V-cycle.
%
%   results = q1_multigrid_dgs()
%   results = q1_multigrid_dgs(N, coarsestN, nu1, nu2, tol, maxCycles)
%
%   Default parameters:
%       N         = 1024
%       coarsestN = 4
%       nu1       = 2
%       nu2       = 2
%       tol       = 1e-8
%       maxCycles = 50
%
%   The code is self-contained and uses matrix-free stencil operations.

    if nargin < 1 || isempty(N),         N = 1024;    end
    if nargin < 2 || isempty(coarsestN), coarsestN = 4; end
    if nargin < 3 || isempty(nu1),       nu1 = 2;     end
    if nargin < 4 || isempty(nu2),       nu2 = 2;     end
    if nargin < 5 || isempty(tol),       tol = 1e-8;  end
    if nargin < 6 || isempty(maxCycles), maxCycles = 50; end

    ratio = N / coarsestN;
    if N < coarsestN || mod(N, coarsestN) ~= 0 || abs(log2(ratio) - round(log2(ratio))) > 1e-12
        error('N must be a power-of-two multiple of coarsestN.');
    end

    [u, v, p, D, f, g, uExact, vExact] = setup_problem(N);

    [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
    r0 = residual_norm(ru, rv, rp, N);
    re = r0;
    er = velocity_error(u, v, uExact, vExact, N);

    residualHistory = zeros(maxCycles + 1, 1);
    errorHistory    = zeros(maxCycles + 1, 1);
    residualHistory(1) = re;
    errorHistory(1)    = er;

    elapsed = 0.0;
    iter = 0;

    while re > tol * r0 && iter < maxCycles
        tic;
        [u, v, p] = vcycle_stokes(u, v, p, f, g, D, nu1, nu2, coarsestN, N);
        p = p - mean(p(:));
        elapsed = elapsed + toc;

        [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
        re = residual_norm(ru, rv, rp, N);
        er = velocity_error(u, v, uExact, vExact, N);

        iter = iter + 1;
        residualHistory(iter + 1) = re;
        errorHistory(iter + 1)    = er;
    end

    fprintf('\n');
    fprintf('Q1: V-cycle multigrid with DGS smoothing\n');
    fprintf('N               = %d\n', N);
    fprintf('coarsestN       = %d\n', coarsestN);
    fprintf('nu1, nu2        = %d, %d\n', nu1, nu2);
    fprintf('iterations      = %d\n', iter);
    fprintf('time (seconds)  = %.6f\n', elapsed);
    fprintf('final residual  = %.6e\n', re);
    fprintf('final error     = %.6e\n', er);

    results = struct();
    results.N = N;
    results.coarsestN = coarsestN;
    results.nu1 = nu1;
    results.nu2 = nu2;
    results.iterations = iter;
    results.time = elapsed;
    results.finalResidual = re;
    results.finalError = er;
    results.u = u;
    results.v = v;
    results.p = p;
    results.residualHistory = residualHistory(1:iter+1);
    results.errorHistory = errorHistory(1:iter+1);
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

function re = residual_norm(ru, rv, rp, N)
%RESIDUAL_NORM Euclidean norm of the full residual, scaled by 1/N.

    re = sqrt(sum(ru(:).^2) + sum(rv(:).^2) + sum(rp(:).^2)) / N;
end

function er = velocity_error(u, v, uExact, vExact, N)
%VELOCITY_ERROR Euclidean velocity error, scaled by 1/N.

    er = sqrt(sum((u(:) - uExact(:)).^2) + sum((v(:) - vExact(:)).^2)) / N;
end

function [u, v, p] = vcycle_stokes(u, v, p, f, g, D, nu1, nu2, coarsestN, N)
%VCYCLE_STOKES Recursive V-cycle for the full Stokes system.

    if N == coarsestN
        for k = 1:100
            [u, v, p] = dgs_step(u, v, p, f, g, D, N);
        end
        return;
    end

    for k = 1:nu1
        [u, v, p] = dgs_step(u, v, p, f, g, D, N);
    end

    [ru, rv, rp] = stokes_residual(u, v, p, f, g, D, N);
    [ruc, rvc, rpc] = restrict_stokes_residual(ru, rv, rp, N);

    [uc, vc, pc] = initialize_stokes_fields(N / 2);
    [uc, vc, pc] = vcycle_stokes(uc, vc, pc, ruc, rvc, rpc, nu1, nu2, coarsestN, N / 2);

    [eu, ev, ep] = prolong_stokes_correction(uc, vc, pc, N);
    u = u + eu;
    v = v + ev;
    p = p + ep;

    for k = 1:nu2
        [u, v, p] = dgs_step(u, v, p, f, g, D, N);
    end
end

function [u, v, p] = initialize_stokes_fields(N)
%INITIALIZE_STOKES_FIELDS Zero initial guess on an N-by-N pressure grid.

    u = zeros(N + 1, N);
    v = zeros(N, N + 1);
    p = zeros(N, N);
end

function [u, v, p] = dgs_step(u, v, p, f, g, D, N)
%DGS_STEP One DGS sweep: momentum relaxation + divergence correction.

    [u, v] = dgs_momentum_relaxation(u, v, p, f, g, N);
    [u, v, p] = dgs_divergence_correction(u, v, p, D, N);
end

function [u, v] = dgs_momentum_relaxation(u, v, p, f, g, N)
%DGS_MOMENTUM_RELAXATION First step of DGS: Gauss--Seidel on momentum.

    h = 1 / N;

    u(1, :) = 0.0;
    for i = 2:N
        j = 1;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j + 1) ...
            - h * (p(i, j) - p(i - 1, j))) / 3;

        for j = 2:N-1
            u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) ...
                - h * (p(i, j) - p(i - 1, j))) / 4;
        end

        j = N;
        u(i, j) = (h * h * f(i, j) + u(i - 1, j) + u(i + 1, j) + u(i, j - 1) ...
            - h * (p(i, j) - p(i - 1, j))) / 3;
    end
    u(N + 1, :) = 0.0;

    v(:, 1) = 0.0;
    for j = 2:N
        i = 1;
        v(i, j) = (h * h * g(i, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1) ...
            - h * (p(i, j) - p(i, j - 1))) / 3;

        for i = 2:N-1
            v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i + 1, j) + v(i, j - 1) + v(i, j + 1) ...
                - h * (p(i, j) - p(i, j - 1))) / 4;
        end

        i = N;
        v(i, j) = (h * h * g(i, j) + v(i - 1, j) + v(i, j - 1) + v(i, j + 1) ...
            - h * (p(i, j) - p(i, j - 1))) / 3;
    end
    v(:, N + 1) = 0.0;
end

function [u, v, p] = dgs_divergence_correction(u, v, p, D, N)
%DGS_DIVERGENCE_CORRECTION Second step of DGS: local divergence repair.

    h = 1 / N;

    for i = 2:N-1
        for j = 2:N-1
            r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
            d = r * h / 4;
            u(i, j)     = u(i, j)     - d;
            u(i + 1, j) = u(i + 1, j) + d;
            v(i, j)     = v(i, j)     - d;
            v(i, j + 1) = v(i, j + 1) + d;

            p(i, j)     = p(i, j)     + r;
            p(i + 1, j) = p(i + 1, j) - r / 4;
            p(i - 1, j) = p(i - 1, j) - r / 4;
            p(i, j + 1) = p(i, j + 1) - r / 4;
            p(i, j - 1) = p(i, j - 1) - r / 4;
        end
    end

    j = N;
    for i = 2:N-1
        r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
        d = r * h / 3;
        u(i, j)     = u(i, j)     - d;
        u(i + 1, j) = u(i + 1, j) + d;
        v(i, j)     = v(i, j)     - d;

        p(i, j)     = p(i, j)     + r;
        p(i + 1, j) = p(i + 1, j) - r / 3;
        p(i - 1, j) = p(i - 1, j) - r / 3;
        p(i, j - 1) = p(i, j - 1) - r / 3;
    end

    i = N;
    for j = 2:N-1
        r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
        d = r * h / 3;
        u(i, j)     = u(i, j)     - d;
        v(i, j)     = v(i, j)     - d;
        v(i, j + 1) = v(i, j + 1) + d;

        p(i, j)     = p(i, j)     + r;
        p(i - 1, j) = p(i - 1, j) - r / 3;
        p(i, j + 1) = p(i, j + 1) - r / 3;
        p(i, j - 1) = p(i, j - 1) - r / 3;
    end

    i = 1;
    for j = 2:N-1
        r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
        d = r * h / 3;
        u(i + 1, j) = u(i + 1, j) + d;
        v(i, j)     = v(i, j)     - d;
        v(i, j + 1) = v(i, j + 1) + d;

        p(i, j)     = p(i, j)     + r;
        p(i + 1, j) = p(i + 1, j) - r / 3;
        p(i, j + 1) = p(i, j + 1) - r / 3;
        p(i, j - 1) = p(i, j - 1) - r / 3;
    end

    j = 1;
    for i = 2:N-1
        r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
        d = r * h / 3;
        u(i, j)     = u(i, j)     - d;
        u(i + 1, j) = u(i + 1, j) + d;
        v(i, j + 1) = v(i, j + 1) + d;

        p(i, j)     = p(i, j)     + r;
        p(i + 1, j) = p(i + 1, j) - r / 3;
        p(i - 1, j) = p(i - 1, j) - r / 3;
        p(i, j + 1) = p(i, j + 1) - r / 3;
    end

    i = 1; j = 1;
    r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
    d = r * h / 2;
    u(i + 1, j) = u(i + 1, j) + d;
    v(i, j + 1) = v(i, j + 1) + d;
    p(i, j)     = p(i, j)     + r;
    p(i, j + 1) = p(i, j + 1) - r / 2;
    p(i + 1, j) = p(i + 1, j) - r / 2;

    i = 1; j = N;
    r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
    d = r * h / 2;
    u(i + 1, j) = u(i + 1, j) + d;
    v(i, j)     = v(i, j)     - d;
    p(i, j)     = p(i, j)     + r;
    p(i + 1, j) = p(i + 1, j) - r / 2;
    p(i, j - 1) = p(i, j - 1) - r / 2;

    i = N; j = 1;
    r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
    d = r * h / 2;
    u(i, j)     = u(i, j)     - d;
    v(i, j + 1) = v(i, j + 1) + d;
    p(i, j)     = p(i, j)     + r;
    p(i - 1, j) = p(i - 1, j) - r / 2;
    p(i, j + 1) = p(i, j + 1) - r / 2;

    i = N; j = N;
    r = (u(i, j) - u(i + 1, j) + v(i, j) - v(i, j + 1)) / h - D(i, j);
    d = r * h / 2;
    u(i, j)     = u(i, j)     - d;
    v(i, j)     = v(i, j)     - d;
    p(i, j)     = p(i, j)     + r;
    p(i - 1, j) = p(i - 1, j) - r / 2;
    p(i, j - 1) = p(i, j - 1) - r / 2;

    u(1, :)     = 0.0;
    u(N + 1, :) = 0.0;
    v(:, 1)     = 0.0;
    v(:, N + 1) = 0.0;
end

function [ruc, rvc, rpc] = restrict_stokes_residual(ru, rv, rp, N)
%RESTRICT_STOKES_RESIDUAL Restrict a Stokes residual to the next coarser grid.

    Nc = N / 2;
    ruc = zeros(Nc + 1, Nc);
    rvc = zeros(Nc, Nc + 1);
    rpc = zeros(Nc, Nc);

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

    for i = 1:Nc
        for j = 1:Nc
            rpc(i, j) = (rp(2*i,   2*j)   + rp(2*i-1, 2*j) + ...
                         rp(2*i,   2*j-1) + rp(2*i-1, 2*j-1)) / 4;
        end
    end
end

function [eu, ev, ep] = prolong_stokes_correction(euc, evc, epc, N)
%PROLONG_STOKES_CORRECTION Prolong a coarse-grid Stokes correction.

    Nc = N / 2;

    eu = zeros(N + 1, N);
    ev = zeros(N, N + 1);
    ep = zeros(N, N);

    for i = 1:Nc
        for j = 1:Nc
            ep(2*i,   2*j)   = epc(i, j);
            ep(2*i-1, 2*j)   = epc(i, j);
            ep(2*i,   2*j-1) = epc(i, j);
            ep(2*i-1, 2*j-1) = epc(i, j);
        end
    end

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