function params = model_to_solver_parameters(p, numerics)
    %MODEL_TO_SOLVER_PARAMETERS Translate manuscript names to solver schema.

    params.SchemaVersion = 3;
    params.SolverVersion = 'publication_model';
    params.Condition = char(p.condition);
    params.R = p.domainRadius;
    params.SourceRadius = p.sourceRadius;
    params.K = p.clearanceRate;
    params.Da = p.ligandDiffusivity;
    params.Dm = p.mesodermDiffusivity;
    params.Alpha = p.upperResponseFraction;
    params.Beta = p.basalChemotacticResponse;
    params.Gamma = p.ingressionRate;
    params.MuM = p.initialCenterRadius;
    params.SigmaM = p.initialRadialWidth;
    params.LigandBoundaryAbsorption = local_boundary_absorption( p.ligandOuterBoundary);
    params.MesodermBoundaryAbsorption = local_boundary_absorption( p.mesodermOuterBoundary);
    params.T = numerics.durationHours;
    params.Nr = numerics.radialPoints;
    params.TimeStep = numerics.timeStepHours;
    params.MaxStoredSteps = numerics.maximumStoredSteps;
end

function absorption = local_boundary_absorption(boundary)
    boundary = lower(char(string(boundary)));
    switch boundary
        case 'dirichlet'
            absorption = 1;
        case {'neumann', 'no_flux', 'no-flux'}
            absorption = 0;
        otherwise
            error('manuscriptModel:BadBoundary', 'Boundary must be ''dirichlet'' or ''no_flux''.');
    end
end
