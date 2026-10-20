function results = run_baseline_model(project, varargin)
    %RUN_BASELINE_MODEL Run or load control and shAPLNR manuscript baselines.

    p = inputParser;
    addRequired(p, 'project', @isstruct);
    addParameter(p, 'Recompute', false, @(x) islogical(x) || isscalar(x));
    parse(p, project, varargin{:});
    conditions = {'control', 'shAPLNR'};
    for i = 1:numel(conditions)
        path = fullfile(project.resultsRoot, [lower(conditions{i}) '_model_results.mat']);
        if logical(p.Results.Recompute) || ~isfile(path)
            result = run_model(conditions{i});
            save(path, 'result', '-v7.3');
        else
            loaded = load(path, 'result');
            result = loaded.result;
            % Normalize cache field aliases when loading.
            fieldAliases = {'m_apical', 'v_apical'};
            upperFields = {'m_upper', 'v_upper'};
            for j = 1:numel(fieldAliases)
                if isfield(result, fieldAliases{j})
                    result.(upperFields{j}) = result.(fieldAliases{j});
                    result = rmfield(result, fieldAliases{j});
                end
            end
            if isfield(result.manuscriptParameters, 'apicalResponseFraction')
                result.manuscriptParameters.upperResponseFraction = result.manuscriptParameters.apicalResponseFraction;
                result.manuscriptParameters = rmfield(result.manuscriptParameters, 'apicalResponseFraction');
            end
        end
        results.(conditions{i}) = result;
        results.paths.(conditions{i}) = path;
    end
end
