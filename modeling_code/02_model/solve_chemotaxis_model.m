function out = solve_chemotaxis_model(params)
    %SOLVE_CHEMOTAXIS_MODEL Validated radially symmetric numerical engine.
    %
    % Ligand is the steady solution of
    %   0 = H(rs-r) - K*a + Da*laplacian_r(a).
    % Upper and basal mesoderm share Dm, undergo chemotaxis with coefficients
    % Alpha*Beta and Beta, and transfer upper-to-basal at rate Gamma.

    arguments
        params (1, 1) struct
    end
    params = local_resolve_boundary_params(params);
    local_validate_params(params);

    R = double(params.R);
    Nr = double(params.Nr);
    rs = linspace(0, R, Nr).';
    dr = rs(2) - rs(1);
    Nt = max(1, round(double(params.T) / double(params.TimeStep)));
    dt = double(params.T) / Nt;
    storeSteps = unique(round(linspace(0, Nt, min(double(params.MaxStoredSteps), Nt + 1))));
    ts = storeSteps(:) .* dt;

    L = local_radial_laplacian_dirichlet(rs, dr);
    source = double(rs <= double(params.SourceRadius));
    steadyMatrix = double(params.K) * speye(Nr) - double(params.Da) * L;
    steadyMatrix(end, :) = 0;
    ligandAbsorption = double(params.LigandBoundaryAbsorption);
    if ligandAbsorption >= 1
        steadyMatrix(end, end) = 1;
    else
        ligandBiot = local_boundary_biot(ligandAbsorption);
        steadyMatrix(end, end - 1) = -1;
        steadyMatrix(end, end) = 1 + ligandBiot * dr / R;
    end
    source(end) = 0;
    a = steadyMatrix \ source;
    a(~isfinite(a)) = 0;

    m0 = exp(-0.5 .* ((rs - double(params.MuM)) ./ double(params.SigmaM)).^2);
    mesodermAbsorption = double(params.MesodermBoundaryAbsorption);
    if mesodermAbsorption >= 1
        m0(end) = 0;
    end
    basalFraction = double(params.InitialBasalFraction);
    mBasal = basalFraction .* m0;
    mUpper = (1 - basalFraction) .* m0;
    initialMass = trapz(rs, 2 * pi * rs .* (mUpper + mBasal));
    if ~(isfinite(initialMass) && initialMass > 0)
        error('run_minimal_chemotaxis:BadInitialMass', 'The initial Gaussian has nonpositive annular mass.');
    end
    mUpper = mUpper ./ initialMass;
    mBasal = mBasal ./ initialMass;

    implicitM = decomposition(local_mesoderm_implicit_matrix( L, Nr, dt, double(params.Dm), dr, R, mesodermAbsorption), 'lu');
    transferFraction = 1 - exp(-double(params.Gamma) * dt);
    nStored = numel(storeSteps);
    stored = local_empty_storage(nStored, Nr);
    storePtr = 1;
    [stored, storePtr] = local_store(stored, storePtr, a, mUpper, mBasal, rs, dr, params);

    idx = 2:(Nr - 1);
    rIp = rs(idx) + dr / 2;
    rIm = rs(idx) - dr / 2;
    lapScale = 1 ./ (rs(idx) * dr);

    for step = 1:Nt
        [advectionBasal, ~, boundaryFluxBasal] = local_chemotaxis_terms( rs, dr, idx, rIp, rIm, lapScale, a, mBasal, double(params.Beta));
        [advectionUpper, ~, boundaryFluxUpper] = local_chemotaxis_terms( rs, dr, idx, rIp, rIm, lapScale, a, mUpper, double(params.Alpha) * double(params.Beta));

        ingressed = transferFraction .* mUpper;
        rhsUpper = mUpper - ingressed - dt .* advectionUpper;
        rhsBasal = mBasal + ingressed - dt .* advectionBasal;
        if mesodermAbsorption >= 1
            rhsUpper(end) = 0;
            rhsBasal(end) = 0;
        else
            boundaryFluxScale = 2 * (R - dr / 2) / (R * dr);
            rhsUpper(end) = rhsUpper(end) + dt * boundaryFluxScale * boundaryFluxUpper;
            rhsBasal(end) = rhsBasal(end) + dt * boundaryFluxScale * boundaryFluxBasal;
        end

        mUpper = implicitM \ rhsUpper;
        mBasal = implicitM \ rhsBasal;
        mUpper = max(mUpper, 0);
        mBasal = max(mBasal, 0);

        if storePtr <= nStored && step == storeSteps(storePtr)
            [stored, storePtr] = local_store(stored, storePtr, a, mUpper, mBasal, rs, dr, params);
        end
    end

    resolvedParams = params;
    resolvedParams.TimeStep = dt;
    resolvedParams.NtComputed = Nt;
    resolvedParams.NtStored = nStored;
    resolvedParams.LigandSourceAmplitude = 1;

    out = struct();
    out.schemaVersion = 2;
    out.solverVersion = char(string(params.SolverVersion));
    out.condition = char(string(params.Condition));
    out.params = resolvedParams;
    out.rs = rs;
    out.ts = ts;
    out.a = stored.a;
    out.m_upper = stored.mUpper;
    out.m_basal = stored.mBasal;
    out.m_surface = stored.mUpper;
    out.m_ingressed = stored.mBasal;
    out.m_mobile = stored.mBasal;
    out.m = stored.mUpper + stored.mBasal;
    out.v_upper = stored.vUpper;
    out.v_basal = stored.vBasal;
    out.v_surface_drift = stored.vUpper;
    out.v_drift = stored.vBasal;
    out.totalMass = stored.totalMass;
    out.basalMassFraction = stored.basalMassFraction;
    out.ingressedMassFraction = stored.basalMassFraction;
end

function storage = local_empty_storage(nStored, Nr)
    blank = nan(nStored, Nr);
    storage = struct( 'a', blank, 'mUpper', blank, 'mBasal', blank, 'vUpper', blank, 'vBasal', blank, 'totalMass', nan(nStored, 1), 'basalMassFraction', nan(nStored, 1));
end

function [storage, nextPtr] = local_store(storage, ptr, a, mUpper, mBasal, rs, dr, params)
    [~, vBasal] = local_chemotaxis_terms_from_grid(a, mBasal, rs, dr, double(params.Beta));
    [~, vUpper] = local_chemotaxis_terms_from_grid( a, mUpper, rs, dr, double(params.Alpha) * double(params.Beta));
    total = mUpper + mBasal;
    totalMass = trapz(rs, 2 * pi * rs .* total);
    basalMass = trapz(rs, 2 * pi * rs .* mBasal);

    storage.a(ptr, :) = a;
    storage.mUpper(ptr, :) = mUpper;
    storage.mBasal(ptr, :) = mBasal;
    storage.vUpper(ptr, :) = vUpper;
    storage.vBasal(ptr, :) = vBasal;
    storage.totalMass(ptr) = totalMass;
    if totalMass > 0
        storage.basalMassFraction(ptr) = basalMass / totalMass;
    end
    nextPtr = ptr + 1;
end

function [advection, velocity] = local_chemotaxis_terms_from_grid(a, density, rs, dr, beta)
    idx = 2:(numel(rs) - 1);
    rIp = rs(idx) + dr / 2;
    rIm = rs(idx) - dr / 2;
    lapScale = 1 ./ (rs(idx) * dr);
    [advection, velocity] = local_chemotaxis_terms( rs, dr, idx, rIp, rIm, lapScale, a, density, beta);
end

function [advection, velocity, boundaryInnerFlux] = local_chemotaxis_terms( rs, dr, idx, rIp, rIm, lapScale, a, density, beta)
    Nr = numel(rs);
    advection = zeros(Nr, 1);
    velocity = zeros(Nr, 1);
    boundaryInnerFlux = 0;
    if beta <= 0 || isempty(idx)
        return;
    end

    density = max(density(:), 0);
    gradientOut = (a(idx + 1) - a(idx)) / dr;
    gradientIn = (a(idx) - a(idx - 1)) / dr;
    upwindOut = density(idx) .* (gradientOut >= 0) + density(idx + 1) .* (gradientOut < 0);
    upwindIn = density(idx - 1) .* (gradientIn > 0) + density(idx) .* (gradientIn <= 0);
    fluxOut = beta .* upwindOut .* gradientOut;
    fluxIn = beta .* upwindIn .* gradientIn;
    boundaryInnerFlux = fluxOut(end);
    advection(idx) = (rIp .* fluxOut - rIm .* fluxIn) .* lapScale;
    velocity(idx) = beta .* 0.5 .* (gradientOut + gradientIn);
    velocity(end) = beta .* (a(end) - a(end - 1)) ./ dr;
    advection(~isfinite(advection)) = 0;
    velocity(~isfinite(velocity)) = 0;
end

function matrix = local_mesoderm_implicit_matrix( L, Nr, dt, diffusivity, dr, R, absorption)
    matrix = speye(Nr) - dt * diffusivity * L;
    if absorption >= 1
        matrix(end, :) = 0;
        matrix(end, end) = 1;
        return;
    end
    biot = local_boundary_biot(absorption);
    permeability = diffusivity * biot / R;
    innerDiffusion = 2 * diffusivity * (R - dr / 2) / (R * dr^2);
    outerLoss = 2 * permeability / dr;
    matrix(end, :) = 0;
    matrix(end, end - 1) = -dt * innerDiffusion;
    matrix(end, end) = 1 + dt * (innerDiffusion + outerLoss);
end

function biot = local_boundary_biot(absorption)
    biot = absorption ./ (1 - absorption);
end

function params = local_resolve_boundary_params(params)
    if ~isfield(params, 'LigandBoundaryAbsorption')
        params.LigandBoundaryAbsorption = 1;
    end
    if ~isfield(params, 'MesodermBoundaryAbsorption')
        params.MesodermBoundaryAbsorption = 1;
    end
    if ~isfield(params, 'InitialBasalFraction')
        params.InitialBasalFraction = 0;
    end
end

function L = local_radial_laplacian_dirichlet(rs, dr)
    Nr = numel(rs);
    idx = 2:(Nr - 1);
    rIp = rs(idx) + dr / 2;
    rIm = rs(idx) - dr / 2;
    mainDiag = zeros(Nr, 1);
    upperDiag = zeros(Nr - 1, 1);
    lowerDiag = zeros(Nr - 1, 1);
    mainDiag(1) = -2 / dr^2;
    upperDiag(1) = 2 / dr^2;
    mainDiag(idx) = -2 / dr^2;
    upperDiag(idx) = rIp ./ (rs(idx) * dr^2);
    lowerDiag(idx - 1) = rIm ./ (rs(idx) * dr^2);
    L = spdiags([[lowerDiag; 0], mainDiag, [0; upperDiag]], -1:1, Nr, Nr);
end

function local_validate_params(params)
    assert(params.R > 0 && params.SourceRadius >= 0 && params.SourceRadius <= params.R);
    assert(params.K > 0 && params.Da > 0 && params.Dm >= 0);
    assert(params.Alpha >= 0 && params.Alpha <= 1 && params.Beta >= 0);
    assert(params.Gamma >= 0 && params.SigmaM > 0);
    assert(params.LigandBoundaryAbsorption >= 0 && params.LigandBoundaryAbsorption <= 1);
    assert(params.MesodermBoundaryAbsorption >= 0 && params.MesodermBoundaryAbsorption <= 1);
    assert(params.InitialBasalFraction >= 0 && params.InitialBasalFraction <= 1);
    assert(params.T > 0 && params.Nr >= 3 && params.TimeStep > 0);
end
