function paths = prepare_experimental_data(project, varargin)

    p = inputParser;
    addRequired(p, 'project', @isstruct);
    addParameter(p, 'Recompute', false, @(x) islogical(x) || isscalar(x));
    parse(p, project, varargin{:});
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    conditions = {'control', 'treatment'};
    stems = {{'ctrl_1', 'ctrl_2', 'ctrl_3', 'ctrl_4'}, {'dox_1', 'dox_2', 'dox_3', 'dox_4'}};
    paths = struct();
    for i = 1:numel(conditions)
        expected = fullfile(project.preparedDataRoot, conditions{i}, 'target_features', 'experimental_target_features.mat');
        if logical(p.Results.Recompute) || ~isfile(expected)
            summaries = cell(numel(stems{i}), 1);
            for j = 1:numel(stems{i})
                fprintf('Preparing offset_24/%s/%s\n', conditions{i}, stems{i}{j});
                summaries{j} = local_prepare_replicate( project.rawDataRoot, conditions{i}, stems{i}{j}, analysis, edgeRadius);
            end
            replicateSummaries = vertcat(summaries{:});
            if ~isfolder(fileparts(expected)), mkdir(fileparts(expected)); end
            save(expected, 'replicateSummaries', '-v7.3');
        end
        paths.(conditions{i}) = expected;
    end
end

function summary = local_prepare_replicate(dataRoot, condition, stem, analysis, edgeRadius)
    files = resolve_dataset_files('offset_24', condition, 'DataRoot', dataRoot, 'FileStem', stem);
    edgeTable = read_data_table(files.edgeFile);
    spotTable = read_data_table(files.spotFile);
    [xPaths, yPaths, tPathsSeconds] = build_track_segments( edgeTable, spotTable, 'MinimumFrames', analysis.minimumFramesPerTrackSegment);

    [centerX, centerY, radius] = local_estimate_geometry({double(spotTable.POSITION_X(:))}, {double(spotTable.POSITION_Y(:))});
    geometry = struct('method', 'fixed_circle', 'centerX', centerX, 'centerY', centerY, 'radius', edgeRadius, 'edgeRadius', edgeRadius, 'measuredRadius', radius, 'edgeRadiusSource', 'override');

    timeOffsetHours = analysis.experimentStartHour;
    tPathsHours = cellfun(@(t) t ./ 3600 + timeOffsetHours, tPathsSeconds, 'UniformOutput', false);
    tMin = min(cellfun(@local_finite_min, tPathsHours));
    tMax = max(cellfun(@local_finite_max, tPathsHours));
    assert(isscalar(tMin) && isscalar(tMax) && isfinite(tMin) && isfinite(tMax) && tMax > tMin, 'modeling:InvalidTimeRange', 'No valid tracking time range for %s/%s.', condition, stem);
    rEdges = linspace(0, edgeRadius, analysis.radialBinCount + 1);
    tGrid = linspace(tMin, tMax, analysis.densityTimePointCount);

    spotTimeHours = double(spotTable.POSITION_T(:)) ./ 3600 + timeOffsetHours;
    spotRadius = hypot(double(spotTable.POSITION_X(:)) - centerX, double(spotTable.POSITION_Y(:)) - centerY);
    density = compute_density_timecourse( spotRadius, spotTimeHours, rEdges, tGrid, 'MinimumBinCount', analysis.minimumCellsPerFrameBin);
    radialVelocity = compute_segment_velocity_maps( xPaths, yPaths, tPathsSeconds, geometry, 'TimeOffsetHours', timeOffsetHours, 'DistanceEdges', rEdges, 'TimeEdgesHours', linspace(analysis.timeWindowHours(1), analysis.timeWindowHours(2), analysis.radialBinCount + 1), 'MinimumCount', analysis.minimumCellsPerFrameBin);

    settings = struct('nR', analysis.radialBinCount, 'nT', analysis.densityTimePointCount, 'rMax', edgeRadius, 'densityMode', 'area', 'timeOffsetHours', timeOffsetHours, 'densityConvention', 'all_unique_spots_with_minimum_bin_occupancy', 'kinematicConvention', 'nonbranching_segments_split_at_divisions', 'kinematicMinimumFrames', analysis.minimumFramesPerTrackSegment, 'geometryConvention', 'fixed_circle', 'edgeRadiusOverride', edgeRadius, 'densityMinimumBinCount', analysis.minimumCellsPerFrameBin, 'velocityMinimumBinCount', analysis.minimumCellsPerFrameBin, 'velocityMapConvention', 'segment_midpoint_radial_speed');
    summary = struct('fileStem', stem, 'files', files, 'density', density, 'radialVelocity', radialVelocity, 'geometry', geometry, 'settings', settings);
end

function [centerX, centerY, radius] = local_estimate_geometry(xPaths, yPaths)
    allX = vertcat(xPaths{:}); allY = vertcat(yPaths{:});
    valid = isfinite(allX) & isfinite(allY);
    allX = allX(valid); allY = allY(valid);
    assert(numel(allX) >= 10, 'modeling:InsufficientCirclePoints', 'Circle fitting requires at least 10 finite positions.');
    boundsCenter = [0.5 * (min(allX) + max(allX)), 0.5 * (min(allY) + max(allY))];
    boundsRadius = max([max(allX) - min(allX), max(allY) - min(allY)]) / 2;
    [centerX, centerY, radius] = local_estimate_ring_center(allX, allY, boundsCenter, boundsRadius);
end

function [centerX, centerY, radius] = local_estimate_ring_center(x, y, boundsCenter, boundsRadius)
    [cx, cy, r] = local_circle_fit(x, y);
    assert(local_circle_fit_is_reasonable(cx, cy, r, boundsCenter, boundsRadius), 'modeling:InvalidCircleFit', 'Circle fit is nonfinite or implausible relative to the observed positions.');
    ri = hypot(x - cx, y - cy);
    residual = abs(ri - r);
    medResidual = median(residual, 'omitnan');
    if isfinite(medResidual) && medResidual > 0
        keep = residual <= max(3 * medResidual, 0.02 * r);
        if nnz(keep) >= 10
            [cx2, cy2, r2] = local_circle_fit(x(keep), y(keep));
            assert(local_circle_fit_is_reasonable(cx2, cy2, r2, boundsCenter, boundsRadius), 'modeling:InvalidCircleFit', 'Refined circle fit is nonfinite or implausible relative to the observed positions.');
            cx = cx2; cy = cy2; r = r2;
        end
    end
    centerX = cx; centerY = cy; radius = r;
end

function [cx, cy, r] = local_circle_fit(x, y)
    A = [2*x(:), 2*y(:), ones(numel(x), 1)]; b = x(:).^2 + y(:).^2;
    assert(rank(A) == 3, 'modeling:InvalidCircleFit', 'Circle fitting requires noncollinear positions.');
    params = A \ b;
    cx = params(1); cy = params(2); r2 = params(3) + cx^2 + cy^2;
    assert(isfinite(r2) && r2 > 0, 'modeling:InvalidCircleFit', 'Circle fit produced an invalid radius.');
    r = sqrt(r2);
end

function tf = local_circle_fit_is_reasonable(cx, cy, r, boundsCenter, boundsRadius)
    tf = isfinite(cx) && isfinite(cy) && isfinite(r) && r > 0;
    if ~tf, return; end
    centerShift = hypot(cx - boundsCenter(1), cy - boundsCenter(2));
    tf = centerShift <= 2 * max(boundsRadius, eps) && r >= 0.1 * max(boundsRadius, eps) && r <= 5 * max(boundsRadius, eps);
end

function v = local_finite_min(x)
    x = x(isfinite(x));
    if isempty(x), v = nan; else, v = min(x); end
end

function v = local_finite_max(x)
    x = x(isfinite(x));
    if isempty(x), v = nan; else, v = max(x); end
end
