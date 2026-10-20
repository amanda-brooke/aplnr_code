function out = run_parameter_sweeps(varargin)
    %RUN_PARAMETER_SWEEPS Compute publication RMSE sweeps.

    cfg = sweep_settings();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    comparisonWindow = analysis.densityIntervalsHours(end, :) - analysis.experimentStartHour;
    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'OutputDir', fullfile(pwd, 'publication_model_sweeps'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'ControlReplicates', struct([]), @isstruct);
    addParameter(p, 'BaseParameters', struct(), @isstruct);
    addParameter(p, 'BaselineId', 'publication_model', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'FeatureFile', 'provided_in_memory', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'LambdaValues', cfg.lambda, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x > 0));
    addParameter(p, 'SourceRadiusValues', cfg.sourceRadius, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0));
    addParameter(p, 'BetaValues', cfg.beta, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0));
    addParameter(p, 'MesodermDiffusivityValues', cfg.mesodermDiffusivity, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0));
    addParameter(p, 'AlphaValues', cfg.alpha, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0) && all(x <= 1));
    addParameter(p, 'GammaValues', cfg.gamma, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0));
    addParameter(p, 'EdgeRadius', edgeRadius, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(p, 'DistanceFromEdgeLimit', [0 edgeRadius], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
    addParameter(p, 'ModelTimeWindowHours', comparisonWindow, @(x) isnumeric(x) && numel(x) == 2 && x(2) >= x(1));
    addParameter(p, 'Recompute', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'ShowProgress', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    parse(p, varargin{:});
    opts = p.Results;
    if isempty(fieldnames(opts.BaseParameters))
        error('run_parameter_sweeps:MissingBaseParameters', 'BaseParameters must contain a complete solver parameter structure.');
    end
    gridLengths = [numel(opts.SourceRadiusValues), numel(opts.LambdaValues), numel(opts.BetaValues), numel(opts.MesodermDiffusivityValues), numel(opts.AlphaValues), numel(opts.GammaValues)];
    if numel(unique(gridLengths)) ~= 1
        error('run_parameter_sweeps:GridSizeMismatch', 'All supplementary sweep axes must have the same number of values.');
    end

    outputDir = char(string(opts.OutputDir));
    baselineId = char(string(opts.BaselineId));
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    resultPath = fullfile(outputDir, [baselineId '_supplement_sweeps.mat']);
    config = local_config(opts);
    if isfile(resultPath) && ~logical(opts.Recompute)
        loaded = load(resultPath, 'out');
        if isequaln(loaded.out.config, config)
            out = loaded.out;
            fprintf('Loaded completed %s supplementary sweeps: %s\n', baselineId, resultPath);
            return;
        end
    end

    replicates = opts.ControlReplicates(:);
    if isempty(replicates)
        error('run_parameter_sweeps:MissingFeatures', 'Provide prepared ControlReplicates to compute sweeps.');
    end
    % FeatureFile records provenance only; data are supplied in memory.
    featureFile = char(string(opts.FeatureFile));
    if isempty(featureFile), featureFile = 'provided_in_memory'; end
    base = opts.BaseParameters;
    out = struct();
    out.schemaVersion = 4;
    out.baselineId = char(string(opts.BaselineId));
    out.base = base;
    out.config = config;
    out.outputDir = outputDir;
    out.featureFile = featureFile;
    out.sourceLambda = local_run_sweep( 'source_radius_lambda', base, replicates, double(opts.SourceRadiusValues(:).'), double(opts.LambdaValues(:).'), 'SourceRadius', 'Lambda', 'r_s (\mum)', '\lambda (\mum)', @local_source_lambda, base.SourceRadius, sqrt(base.Da / base.K), [baselineId '_density_rmse_sourceRadius_lambda'], outputDir, config, logical(opts.Recompute), logical(opts.ShowProgress));
    out.betaDiffusivity = local_run_sweep( 'beta_mesoderm_diffusivity', base, replicates, double(opts.BetaValues(:).'), double(opts.MesodermDiffusivityValues(:).'), 'Beta', 'Dm', '\beta (\mum^2/h)', 'D_m (\mum^2/h)', @local_beta_diffusivity, base.Beta, base.Dm, [baselineId '_density_rmse_beta_mesodermDiffusivity'], outputDir, config, logical(opts.Recompute), logical(opts.ShowProgress));
    out.alphaGamma = local_run_sweep( 'alpha_gamma', base, replicates, double(opts.AlphaValues(:).'), double(opts.GammaValues(:).'), 'Alpha', 'Gamma', '\alpha', '\gamma (h^{-1})', @local_alpha_gamma, base.Alpha, base.Gamma, [baselineId '_density_rmse_alpha_gamma'], outputDir, config, logical(opts.Recompute), logical(opts.ShowProgress));
    out.metadata = local_metadata(base, config, featureFile, opts);
    out.resultPath = resultPath;
    save(resultPath, 'out', '-v7.3');
    local_write_json(fullfile(outputDir, [baselineId '_supplement_sweeps_metadata.json']), out.metadata);
    fprintf('Saved completed %s supplementary sweeps: %s\n', baselineId, resultPath);
end

function sweep = local_run_sweep(id, base, replicates, xValues, yValues, xName, yName, xLabel, yLabel, mutateFunction, baselineX, baselineY, fileBase, outputDir, config, recompute, showProgress)
    checkpointPath = fullfile(outputDir, [fileBase '_checkpoint.mat']);
    nX = numel(xValues);
    nY = numel(yValues);
    rawRmse = nan(nY, nX);
    relativeRmse = nan(nY, nX);
    status = strings(nY, nX);
    errorMessage = strings(nY, nX);
    if isfile(checkpointPath) && ~recompute
        loaded = load(checkpointPath, 'checkpoint');
        checkpoint = loaded.checkpoint;
        if isequal(checkpoint.xValues, xValues) && isequal(checkpoint.yValues, yValues) && local_checkpoint_configs_match(checkpoint.config, config)
            rawRmse = checkpoint.rawRmse;
            relativeRmse = checkpoint.relativeRmse;
            status = checkpoint.status;
            errorMessage = checkpoint.errorMessage;
            fprintf('Resuming %s from %s\n', id, checkpointPath);
        end
    end

    started = tic;
    nTotal = nX * nY;
    for iy = 1:nY
        for ix = 1:nX
            if status(iy, ix) == "ok"
                continue;
            end
            completed = nnz(status == "ok") + nnz(status == "failed") + 1;
            if showProgress
                fprintf('%s %d/%d: %s=%.6g, %s=%.6g\n', id, completed, nTotal, xName, xValues(ix), yName, yValues(iy));
            end
            try
                params = mutateFunction(base, xValues(ix), yValues(iy));
                model = solve_chemotaxis_model(params);
                metrics = compare_model_experiment_density(model, replicates, 'EdgeRadius', config.edgeRadius, 'DistanceFromEdgeLimit', config.distanceFromEdgeLimit, 'ModelTimeWindowHours', config.modelTimeWindowHours);
                rawRmse(iy, ix) = metrics.rawRmse;
                relativeRmse(iy, ix) = metrics.relativeRmse;
                status(iy, ix) = "ok";
                errorMessage(iy, ix) = "";
            catch ME
                status(iy, ix) = "failed";
                errorMessage(iy, ix) = string(ME.message);
                warning('run_parameter_sweeps:RunFailed', '%s failed at %s=%.6g, %s=%.6g: %s', id, xName, xValues(ix), yName, yValues(iy), ME.message);
            end
        end
        checkpoint = struct('config', config, 'xValues', xValues, 'yValues', yValues, 'rawRmse', rawRmse, 'relativeRmse', relativeRmse, 'status', status, 'errorMessage', errorMessage);
        save(checkpointPath, 'checkpoint');
    end
    elapsedSeconds = toc(started);

    [bestRelative, bestLinearIndex] = min(relativeRmse(:), [], 'omitnan');
    if isempty(bestLinearIndex) || ~isfinite(bestRelative)
        best = struct('xValue', NaN, 'yValue', NaN, 'rawRmse', NaN, 'relativeRmse', NaN);
    else
        [bestY, bestX] = ind2sub(size(relativeRmse), bestLinearIndex);
        best = struct('xValue', xValues(bestX), 'yValue', yValues(bestY), 'rawRmse', rawRmse(bestY, bestX), 'relativeRmse', bestRelative);
    end

    resultTable = local_result_table(xValues, yValues, xName, yName, rawRmse, relativeRmse, status, errorMessage);
    csvPath = fullfile(outputDir, [fileBase '.csv']);
    writetable(resultTable, csvPath);
    sweep = struct('id', id, 'xName', xName, 'yName', yName, 'xLabel', xLabel, 'yLabel', yLabel, 'xValues', xValues, 'yValues', yValues, 'baselineX', baselineX, 'baselineY', baselineY, 'rawRmse', rawRmse, 'relativeRmse', relativeRmse, 'status', status, 'errorMessage', errorMessage, 'best', best, 'elapsedSeconds', elapsedSeconds, 'resultTable', resultTable, 'csvPath', csvPath, 'fileBase', fileBase);
    fprintf('%s best relative RMSE %.4g at %s=%.6g, %s=%.6g\n', id, best.relativeRmse, xName, best.xValue, yName, best.yValue);
end

function tableOut = local_result_table(xValues, yValues, xName, yName, rawRmse, relativeRmse, status, errorMessage)
    [xGrid, yGrid] = meshgrid(xValues, yValues);
    tableOut = table(repmat(string(xName), numel(xGrid), 1), repmat(string(yName), numel(xGrid), 1), xGrid(:), yGrid(:), rawRmse(:), relativeRmse(:), status(:), errorMessage(:), 'VariableNames', {'parameterX', 'parameterY', 'xValue', 'yValue', 'densityRmse', 'densityRelativeRmse', 'status', 'errorMessage'});
    tableOut = sortrows(tableOut, 'densityRelativeRmse', 'ascend', 'MissingPlacement', 'last');
end

function params = local_source_lambda(params, sourceRadius, lambda)
    params.SourceRadius = sourceRadius;
    params.Da = lambda.^2 .* params.K;
end

function params = local_beta_diffusivity(params, beta, diffusivity)
    params.Beta = beta;
    params.Dm = diffusivity;
end

function params = local_alpha_gamma(params, alpha, gamma)
    params.Alpha = alpha;
    params.Gamma = gamma;
end

function config = local_config(opts)
    % Retain the fixed convention fields so existing publication caches match.
    config = struct( 'baselineId', char(string(opts.BaselineId)), 'baseParameters', opts.BaseParameters, 'edgeRadius', double(opts.EdgeRadius), 'distanceFromEdgeLimit', double(opts.DistanceFromEdgeLimit(:).'), 'modelTimeWindowHours', double(opts.ModelTimeWindowHours(:).'), 'timeWindowPostBMPHours', double(opts.ModelTimeWindowHours(:).') + 24, 'areaWeighting', false, 'useReplicateEdgeRadius', false, 'lambdaValues', double(opts.LambdaValues(:).'), 'sourceRadiusValues', double(opts.SourceRadiusValues(:).'), 'betaValues', double(opts.BetaValues(:).'), 'mesodermDiffusivityValues', double(opts.MesodermDiffusivityValues(:).'), 'alphaValues', double(opts.AlphaValues(:).'), 'gammaValues', double(opts.GammaValues(:).'));
end

function tf = local_checkpoint_configs_match(first, second)
    comparisonFields = {'edgeRadius', 'distanceFromEdgeLimit', 'modelTimeWindowHours', 'areaWeighting', 'useReplicateEdgeRadius'};
    tf = all(isfield(first, comparisonFields)) && all(isfield(second, comparisonFields));
    for i = 1:numel(comparisonFields)
        fieldName = comparisonFields{i};
        tf = tf && isequaln(first.(fieldName), second.(fieldName));
    end
end

function metadata = local_metadata(base, config, featureFile, opts)
    metadata = struct( 'schemaVersion', 4, 'baselineId', char(string(opts.BaselineId)), 'condition', 'control', 'generatedAt', char(datetime('now', 'TimeZone', 'America/Los_Angeles', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')), 'solverVersion', base.SolverVersion, 'ligandConvention', 'unit Heaviside source', 'betaConvention', 'v_b=Beta*da/dr; v_a=Alpha*Beta*da/dr', 'betaRange', [min(opts.BetaValues) max(opts.BetaValues)], 'featureFile', featureFile, 'comparison', config, 'baseParameters', base);
end

function local_write_json(path, value)
    fid = fopen(path, 'w');
    if fid < 0
        error('run_parameter_sweeps:MetadataWriteFailed', 'Could not write %s.', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
end
