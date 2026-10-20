function result = run_model(condition, varargin)
    %RUN_MODEL Run control or shAPLNR with optional scientific-parameter overrides.

    p = inputParser;
    addRequired(p, 'condition', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Overrides', struct(), @isstruct);
    parse(p, condition, varargin{:});

    parameters = model_parameters(condition);
    overrides = p.Results.Overrides;
    names = fieldnames(overrides);
    unknown = setdiff(names, fieldnames(parameters));
    if ~isempty(unknown)
        error('manuscriptModel:UnknownOverride', 'Unknown parameter override(s): %s.', strjoin(unknown, ', '));
    end
    for i = 1:numel(names)
        parameters.(names{i}) = overrides.(names{i});
    end
    result = solve_chemotaxis_model( model_to_solver_parameters(parameters, numerical_settings()));
    result.manuscriptParameters = parameters;
    result.manuscriptNumerics = numerical_settings();
end
