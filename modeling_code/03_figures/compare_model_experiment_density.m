function [metrics, modelProfile] = compare_model_experiment_density(model, replicates, varargin)

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'EdgeRadius', 350, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(p, 'DistanceFromEdgeLimit', [0 350], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
    addParameter(p, 'ModelTimeWindowHours', [24 30], @(x) isnumeric(x) && numel(x) == 2 && x(2) >= x(1));
    parse(p, varargin{:});
    opts = p.Results;

    [rModel, tModel, modelDensity] = model_areal_density_profile(model);
    replicates = replicates(:);
    repRaw = nan(numel(replicates), 1);
    repRelative = nan(numel(replicates), 1);

    for i = 1:numel(replicates)
        density = replicates(i).density;
        [rExp, tExp, expDensity] = experimental_areal_density_profile(density);
        tRel = tExp - min(tExp);

        timeWindow = double(opts.ModelTimeWindowHours(:).');
        okT = tRel >= max(min(tModel), timeWindow(1)) - 1e-9 & tRel <= min(max(tModel), timeWindow(2)) + 1e-9;
        xExp = double(opts.EdgeRadius) - rExp;
        xModel = double(opts.EdgeRadius) - rModel;
        [xModel, modelOrder] = sort(xModel);
        distanceLimit = double(opts.DistanceFromEdgeLimit(:).');
        okR = xExp >= distanceLimit(1) & xExp <= distanceLimit(2) & xExp >= min(xModel) & xExp <= max(xModel);
        if nnz(okT) < 2 || nnz(okR) < 2
            continue;
        end

        modelAtTimes = interp1(tModel(:), modelDensity, tRel(okT), 'linear', 'extrap');
        modelAtExp = interp1(xModel, modelAtTimes(:, modelOrder).', xExp(okR), 'linear', nan).';
        expAtExp = expDensity(okT, okR);
        [repRaw(i), repRelative(i)] = local_rmse(modelAtExp, expAtExp);
    end

    metrics = struct();
    metrics.rawRmse = mean(repRaw, 'omitnan');
    metrics.relativeRmse = mean(repRelative, 'omitnan');
    metrics.replicateRawRmse = repRaw;
    metrics.replicateRelativeRmse = repRelative;
    metrics.convention = local_convention(opts);
    modelProfile = struct('r', rModel, 't', tModel, 'density', modelDensity);
end

function [raw, relative] = local_rmse(modelValues, expValues)
    ok = isfinite(modelValues) & isfinite(expValues);
    if ~any(ok(:))
        raw = nan;
        relative = nan;
        return;
    end
    differences = modelValues - expValues;
    raw = sqrt(sum(differences(ok).^2) / nnz(ok));
    expScale = sqrt(sum(expValues(ok).^2) / nnz(ok));
    if isfinite(expScale) && expScale > 0
        relative = raw / expScale;
    else
        relative = nan;
    end
end

function convention = local_convention(opts)
    convention = struct( 'coordinate', 'fixed_distance_from_edge', 'edgeRadiusUm', double(opts.EdgeRadius), 'distanceFromEdgeLimitUm', double(opts.DistanceFromEdgeLimit(:).'), 'modelTimeWindowHours', double(opts.ModelTimeWindowHours(:).'), 'experimentalTimeReference', 'hours_since_first_replicate_frame', 'areaWeighting', false, 'useReplicateEdgeRadius', false, 'modelDensityNormalization', 'unit_annular_mass_per_stored_time', 'experimentalDensity', 'areal_density', 'aggregation', 'mean_of_replicate_rmse', 'relativeScale', 'experimental_rms_with_same_weights');
end
