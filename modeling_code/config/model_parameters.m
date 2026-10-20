function p = model_parameters(condition)
    % Distances: um; time: h; diffusivities and chemotactic response: um^2/h.
    % Ligand production is one in arbitrary concentration units per hour.

    if nargin < 1 || isempty(condition), condition = 'control'; end
    condition = lower(char(string(condition)));
    if ~ismember(condition, {'control', 'shaplnr'})
        error('manuscriptModel:BadCondition', 'Condition must be ''control'' or ''shAPLNR''.');
    end
    p.condition = condition;
    p.domainRadius = 350;
    p.sourceRadius = 220;
    p.clearanceRate = 5;
    p.ligandDiffusivity = 1e4;
    p.mesodermDiffusivity = 30;
    p.upperResponseFraction = 0.5;
    p.basalChemotacticResponse = 5000;
    p.ingressionRate = 0.025;
    p.ligandOuterBoundary = "dirichlet";
    p.mesodermOuterBoundary = "dirichlet";
    if strcmp(condition, 'control')
        p.initialCenterRadius = 320;
        p.initialRadialWidth = 11;
    else
        p.initialCenterRadius = 305;
        p.initialRadialWidth = 12;
        p.basalChemotacticResponse = 0;
    end
end
