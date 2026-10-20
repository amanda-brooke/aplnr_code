function paths = plot_velocity_comparison( controlResult, treatmentResult, controlReplicates, treatmentReplicates, varargin)
    %PLOT_VELOCITY_COMPARISON Compare unsmoothed, equal-frame mean radial velocities.

    style = figure_style();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'controlResult', @isstruct);
    addRequired(p, 'treatmentResult', @isstruct);
    addRequired(p, 'controlReplicates', @isstruct);
    addRequired(p, 'treatmentReplicates', @isstruct);
    addParameter(p, 'TimeWindowHours', analysis.timeWindowHours, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'IntervalTimeWindowsHours', zeros(0, 2), @(x) isnumeric(x) && size(x, 2) == 2 && all(isfinite(x), 'all') && all(x(:, 2) > x(:, 1)) && size(x, 1) <= 5);
    addParameter(p, 'IntervalXLimits', style.supplementVelocityXLimits, @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1)));
    addParameter(p, 'IntervalModelLineStyle', style.intervalVelocityModelLineStyle, @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'IntervalModelColor', style.intervalVelocityModelColor, @(x) isnumeric(x) && numel(x) == 3 && all(isfinite(x)));
    addParameter(p, 'IntervalFigureSize', style.intervalVelocityFigureSize, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && all(x > 0));
    addParameter(p, 'IntervalConditionLabels', style.conditionLabels, @(x) iscellstr(x) && numel(x) == 2); %#ok<ISCLSTR>
    addParameter(p, 'OutputDir', fullfile(pwd, 'minimal_chemotaxis_figures'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'DataRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'MinimumFramesPerSegment', analysis.minimumFramesPerTrackSegment, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && mod(x, 1) == 0);
    addParameter(p, 'MinimumReplicates', analysis.velocity.minimumReplicates, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && mod(x, 1) == 0);
    addParameter(p, 'RadialBinWidth', edgeRadius / analysis.radialBinCount, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'EdgeRadius', edgeRadius, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'VelocityLimits', style.velocityLimits, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'MinimumDisplacementsPerFrameBin', analysis.velocity.minimumDisplacementsPerFrameBin, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && mod(x, 1) == 0);
    addParameter(p, 'MinimumValidFramesPerBin', analysis.velocity.minimumValidFramesPerReplicateBin, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && mod(x, 1) == 0);
    parse(p, controlResult, treatmentResult, controlReplicates, treatmentReplicates, varargin{:});
    opts = p.Results;

    if ~isempty(opts.IntervalTimeWindowsHours)
        paths = local_interval_overlay(controlResult, treatmentResult, controlReplicates, treatmentReplicates, opts);
        return;
    end

    timeWindow = double(opts.TimeWindowHours(:).');
    controlExperiment = local_read_experimental_profiles(controlReplicates, timeWindow, opts);
    treatmentExperiment = local_read_experimental_profiles(treatmentReplicates, timeWindow, opts);
    controlModel = local_model_profile(controlResult, timeWindow, controlExperiment.x);
    treatmentModel = local_model_profile(treatmentResult, timeWindow, treatmentExperiment.x);

    visibility = 'off';
    if logical(opts.Visible)
        visibility = 'on';
    end

    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.velocityFigureSize]);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    local_panel(nexttile(layout), controlExperiment, controlModel, 'Control', double(opts.MinimumReplicates), double(opts.EdgeRadius), double(opts.VelocityLimits), style.controlColor);
    local_panel(nexttile(layout), treatmentExperiment, treatmentModel, 'shAPLNR', double(opts.MinimumReplicates), double(opts.EdgeRadius), double(opts.VelocityLimits), style.treatmentColor);

    outputDir = char(string(opts.OutputDir));
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    baseName = sprintf('figure_model_experiment_inward_velocity_%g_%gh_equal_frame_mean_minFrames%d', timeWindow(1), timeWindow(2), round(double(opts.MinimumValidFramesPerBin)));
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    paths.fig = fullfile(outputDir, [baseName '.fig']);
    savefig(fig, paths.fig);
    exportgraphics(fig, paths.png, 'Resolution', double(opts.Resolution), 'BackgroundColor', 'white');
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    paths.timeWindowHours = timeWindow;
    paths.smoothingBins = 1;
    paths.minimumDisplacementsPerBin = 0;
    paths.minimumFractionOfPeakBinCount = 0;
    paths.minimumReplicates = double(opts.MinimumReplicates);
    paths.modelMinimumDensityFraction = 0;
    paths.edgeRadius = double(opts.EdgeRadius);
    paths.velocityLimits = double(opts.VelocityLimits);
    paths.controlExperimentalProfiles = controlExperiment.replicates;
    paths.treatmentExperimentalProfiles = treatmentExperiment.replicates;
    paths.controlModelProfile = controlModel;
    paths.treatmentModelProfile = treatmentModel;
    paths.distanceFromEdge = controlExperiment.x;
    paths.experimentalBinStatistic = 'mean';
    paths.experimentalTemporalWeighting = 'equal_frame';
    paths.minimumDisplacementsPerFrameBin = double(opts.MinimumDisplacementsPerFrameBin);
    paths.minimumValidFramesPerBin = double(opts.MinimumValidFramesPerBin);
    paths.minimumFramesPerSegment = double(opts.MinimumFramesPerSegment);

    if ~logical(opts.Visible)
        close(fig);
    else
        figure(fig);
    end
end

function paths = local_interval_overlay(controlResult, treatmentResult, controlReplicates, treatmentReplicates, opts)
    style = figure_style();
    intervals = double(opts.IntervalTimeWindowsHours);
    visibility = 'off';
    if logical(opts.Visible), visibility = 'on'; end
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 double(opts.IntervalFigureSize)], 'InvertHardcopy', 'off');
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    conditions = opts.IntervalConditionLabels;
    modelResults = {controlResult, treatmentResult};
    replicateSets = {controlReplicates, treatmentReplicates};
    for panel = 1:2
        intervalColors = local_time_colors(size(intervals, 1));
        ax = nexttile(layout);
        set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        hold(ax, 'on');
        for i = 1:size(intervals, 1)
            experiment = local_read_experimental_profiles( replicateSets{panel}, intervals(i, :), opts);
            meanProfile = mean(experiment.replicates, 2, 'omitnan');
            lowerProfile = min(experiment.replicates, [], 2, 'omitnan');
            upperProfile = max(experiment.replicates, [], 2, 'omitnan');
            support = sum(isfinite(experiment.replicates), 2) >= double(opts.MinimumReplicates);
            meanProfile(~support) = NaN;
            lowerProfile(~support) = NaN;
            upperProfile(~support) = NaN;
            local_band(ax, experiment.x, lowerProfile, upperProfile, intervalColors(i, :), style.intervalVelocityBandAlpha);
            label = sprintf('%g-%g h', intervals(i, 1), intervals(i, 2));
            plot(ax, experiment.x, meanProfile, '-', 'Color', intervalColors(i, :), 'LineWidth', style.intervalVelocityLineWidth, 'DisplayName', label);
        end
        modelWindow = [intervals(1, 1) intervals(end, 2)];
        model = local_model_profile(modelResults{panel}, modelWindow, experiment.x);
        plot(ax, experiment.x, model, char(string(opts.IntervalModelLineStyle)), 'Color', double(opts.IntervalModelColor), 'LineWidth', style.intervalVelocityLineWidth, 'DisplayName', sprintf('Model %g-%g h', modelWindow(1), modelWindow(2)));
        yline(ax, 0, ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        hold(ax, 'off');
        intervalXLimits = double(opts.IntervalXLimits);
        if isempty(intervalXLimits)
            intervalXLimits = [0 double(opts.EdgeRadius)];
        end
        xlim(ax, intervalXLimits);
        ylim(ax, double(opts.VelocityLimits));
        xticks(ax, [intervalXLimits(1) mean(intervalXLimits) intervalXLimits(2)]);
        velocityTicks = -5:5:10;
        yticks(ax, velocityTicks(velocityTicks >= opts.VelocityLimits(1) & velocityTicks <= opts.VelocityLimits(2)));
        xlabel(ax, 'Distance from edge (\mum)');
        ylabel(ax, 'Inward velocity (\mum/h)');
        ax.XLabel.FontSize = style.labelFontSize;
        ax.YLabel.FontSize = style.labelFontSize;
        title(ax, conditions{panel}, 'FontWeight', 'normal', 'Color', 'k', 'FontSize', style.labelFontSize);
        set(ax, 'FontName', style.fontName, 'FontSize', style.supplementTickFontSize, 'LineWidth', style.supplementAxesLineWidth, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'TickDir', 'in', 'Layer', 'top');
    end
    legendHandle = legend(layout.Children(end), 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'Box', 'off');
    legendHandle.Layout.Tile = 'south';
    legendHandle.TextColor = 'k';
    legendHandle.Color = 'w';
    legendHandle.FontSize = style.legendFontSize;

    outputDir = char(string(opts.OutputDir));
    if ~isfolder(outputDir), mkdir(outputDir); end
    baseName = sprintf('figure_model_experiment_inward_velocity_intervals_%g_%gh', intervals(1, 1), intervals(end, 2));
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    exportgraphics(fig, paths.png, 'Resolution', double(opts.Resolution), 'BackgroundColor', 'white');
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    paths.intervalsHours = intervals;
    paths.experimentalLineStyle = '-';
    paths.modelLineStyle = char(string(opts.IntervalModelLineStyle));
    paths.xLimits = double(opts.IntervalXLimits);
    paths.minimumFramesPerSegment = double(opts.MinimumFramesPerSegment);
    if ~logical(opts.Visible), close(fig); end
end

function experiment = local_read_experimental_profiles(replicates, timeWindow, opts)
    % Keep the profile calculation and its conventions shared across both figures.
    experiment = local_experimental_profiles( replicates, timeWindow, double(opts.RadialBinWidth), double(opts.EdgeRadius), double(opts.MinimumDisplacementsPerFrameBin), double(opts.MinimumValidFramesPerBin), char(string(opts.DataRoot)), double(opts.MinimumFramesPerSegment));
end

function out = local_experimental_profiles(replicates, timeWindow, radialBinWidth, edgeRadius, minimumPerFrameBin, minimumValidFrames, dataRoot, minimumFramesPerSegment)
    edges = (0:radialBinWidth:edgeRadius).';
    if edges(end) < edgeRadius, edges(end + 1) = edgeRadius; end
    x = 0.5 .* (edges(1:end-1) + edges(2:end));
    profiles = nan(numel(x), numel(replicates));
    binCounts = zeros(numel(x), numel(replicates));
    validFrameCounts = zeros(numel(x), numel(replicates));
    lineageCounts = zeros(numel(replicates), 1);
    uniqueSpotCounts = zeros(numel(replicates), 1);
    segmentCounts = zeros(numel(replicates), 1);
    contributingSegmentCounts = zeros(numel(replicates), 1);
    validDisplacementCounts = zeros(numel(replicates), 1);
    for j = 1:numel(replicates)
        spotFile = local_resolve_relocated_data_file( replicates(j).files.spotFile, 'spots', dataRoot);
        edgeFile = local_resolve_relocated_data_file( replicates(j).files.edgeFile, 'edges', dataRoot);
        spots = readtable(spotFile, 'VariableNamingRule', 'preserve');
        edgesTable = readtable(edgeFile, 'VariableNamingRule', 'preserve');
        [xTracks, yTracks, tTracks] = build_track_segments(edgesTable, spots, 'MinimumFrames', minimumFramesPerSegment);
        uniqueSpotCounts(j) = nnz(isfinite(double(spots.ID)));
        if ismember('TRACK_ID', spots.Properties.VariableNames)
            lineageCounts(j) = numel(unique(double(spots.TRACK_ID( isfinite(double(spots.TRACK_ID))))));
        end
        segmentCounts(j) = numel(xTracks);
        distanceFromEdge = [];
        inwardVelocity = [];
        displacementTime = [];
        for i = 1:numel(xTracks)
            dt = diff(tTracks{i}) ./ 3600;
            radius = hypot(xTracks{i} - replicates(j).geometry.centerX, yTracks{i} - replicates(j).geometry.centerY);
            radialMidpoint = 0.5 .* (radius(1:end-1) + radius(2:end));
            edgeDistance = edgeRadius - radialMidpoint;
            timeMidpoint = 0.5 .* (tTracks{i}(1:end-1) + tTracks{i}(2:end)) ./ 3600 + 24;
            dx = diff(xTracks{i}); dy = diff(yTracks{i});
            xMid = 0.5 .* (xTracks{i}(1:end-1) + xTracks{i}(2:end));
            yMid = 0.5 .* (yTracks{i}(1:end-1) + yTracks{i}(2:end));
            cx = repmat(replicates(j).geometry.centerX, size(xMid));
            cy = repmat(replicates(j).geometry.centerY, size(yMid));
            outwardX = (xMid - cx) ./ max(radialMidpoint, eps);
            outwardY = (yMid - cy) ./ max(radialMidpoint, eps);
            velocity = -(dx .* outwardX + dy .* outwardY) ./ dt;
            valid = isfinite(velocity) & isfinite(edgeDistance) & isfinite(dt) & dt > 0 & timeMidpoint >= timeWindow(1) & timeMidpoint <= timeWindow(2) & edgeDistance >= edges(1) & edgeDistance <= edges(end);
            if any(valid)
                contributingSegmentCounts(j) = contributingSegmentCounts(j) + 1;
                validDisplacementCounts(j) = validDisplacementCounts(j) + nnz(valid);
            end
            distanceFromEdge = [distanceFromEdge; edgeDistance(valid)]; %#ok<AGROW>
            inwardVelocity = [inwardVelocity; velocity(valid)]; %#ok<AGROW>
            displacementTime = [displacementTime; timeMidpoint(valid)]; %#ok<AGROW>
        end
        bins = discretize(distanceFromEdge, edges);
        valid = isfinite(bins) & isfinite(inwardVelocity);
        counts = accumarray(bins(valid), 1, [numel(x), 1], @sum, 0);
        binCounts(:, j) = counts;
        [profile, validFrameCount] = local_equal_frame_profile(bins(valid), displacementTime(valid), inwardVelocity(valid), numel(x), minimumPerFrameBin, minimumValidFrames);
        validFrameCounts(:, j) = validFrameCount;
        profiles(:, j) = profile;
    end

    function path = local_resolve_relocated_data_file(path, dataKind, dataRoot)
        path = char(string(path));
        if isfile(path), return; end
        [~, fileName, extension] = fileparts(path);
        condition = 'control';
        if contains(lower(path), [filesep 'treatment' filesep]) || contains(lower(fileName), 'dox')
            condition = 'treatment';
        end
        if ~isempty(dataRoot)
            datasetRoot = char(string(dataRoot));
        else
            project = model_pipeline_setup();
            datasetRoot = project.rawDataRoot;
        end
        candidate = fullfile(datasetRoot, condition, dataKind, [fileName extension]);
        if ~isfile(candidate)
            error('plot_velocity_comparison:MissingDataFile', 'Could not resolve relocated data file %s.', path);
        end
        path = candidate;
    end

    out = struct('x', x(:), 'replicates', profiles, 'binCounts', binCounts, 'validFrameCounts', validFrameCounts, 'lineageCounts', lineageCounts, 'uniqueSpotCounts', uniqueSpotCounts, 'segmentCounts', segmentCounts, 'contributingSegmentCounts', contributingSegmentCounts, 'validDisplacementCounts', validDisplacementCounts);
end

function [profile, validFrameCount] = local_equal_frame_profile( radialBin, timeHours, velocity, nRadialBins, minimumPerFrameBin, minimumValidFrames)
    [frameTimes, ~, frameBin] = unique(timeHours);
    nFrames = numel(frameTimes);
    counts = accumarray([radialBin, frameBin], 1, [nRadialBins, nFrames], @sum, 0);
    % Average within each frame/bin, then weight the retained frames equally.
    sums = accumarray([radialBin, frameBin], velocity, [nRadialBins, nFrames], @sum, 0);
    frameValues = sums ./ max(counts, 1);
    frameValues(counts < minimumPerFrameBin) = NaN;
    profile = mean(frameValues, 2, 'omitnan');
    validFrameCount = sum(isfinite(frameValues), 2);
    profile(validFrameCount < minimumValidFrames) = NaN;
end

function velocity = local_model_profile( model, timeWindow, xTarget)
    velocity = model_inward_velocity_profile(model, timeWindow, xTarget);
end

function local_panel(ax, experiment, model, conditionLabel, minimumReplicates, edgeRadius, velocityLimits, conditionColor)
    style = figure_style();
    profiles = experiment.replicates;
    meanProfile = mean(profiles, 2, 'omitnan');
    lowerProfile = min(profiles, [], 2, 'omitnan');
    upperProfile = max(profiles, [], 2, 'omitnan');
    support = sum(isfinite(profiles), 2) >= minimumReplicates;
    meanProfile(~support) = NaN;
    lowerProfile(~support) = NaN;
    upperProfile(~support) = NaN;

    hold(ax, 'on');
    yline(ax, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', style.axesLineWidth, 'HandleVisibility', 'off');
    valid = isfinite(lowerProfile) & isfinite(upperProfile);
    fill(ax, [experiment.x(valid); flipud(experiment.x(valid))], [lowerProfile(valid); flipud(upperProfile(valid))], conditionColor, 'FaceAlpha', style.velocityBandAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    for i = 1:size(profiles, 2)
        plot(ax, experiment.x, profiles(:, i), '--', 'Color', conditionColor + (1 - conditionColor) .* style.velocityReplicateFade, 'LineWidth', style.velocityReplicateLineWidth, 'HandleVisibility', 'off');
    end
    plot(ax, experiment.x, meanProfile, '--', 'Color', conditionColor, 'LineWidth', style.profileLineWidth, 'DisplayName', 'Experiment');
    plot(ax, experiment.x, model, '-', 'Color', conditionColor, 'LineWidth', style.profileLineWidth, 'DisplayName', 'Model');
    hold(ax, 'off');
    set(ax, 'XLim', [0 edgeRadius], 'YLim', velocityLimits, 'XTick', [0 edgeRadius / 2 edgeRadius], 'YTick', [-5 0 5 10], 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'LineWidth', style.axesLineWidth, 'TickDir', 'in', 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Inward velocity (\mum/h)');
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
    title(ax, conditionLabel, 'FontWeight', 'normal', 'FontSize', style.labelFontSize, 'Color', 'k');
    legendHandle = legend(ax, {'Experiment', 'Model'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2, 'Box', 'off');
    set(legendHandle, 'TextColor', 'k', 'Color', 'w', 'FontSize', style.legendFontSize);
end

function colors = local_time_colors(count)
    style = figure_style();
    palette = style.intervalColors;
    colors = palette(mod(0:count-1, size(palette, 1)) + 1, :);
end

function local_band(ax, x, lowerProfile, upperProfile, color, alphaValue)
    x = x(:);
    lowerProfile = lowerProfile(:);
    upperProfile = upperProfile(:);
    ok = isfinite(x) & isfinite(lowerProfile) & isfinite(upperProfile);
    if nnz(ok) < 2, return; end
    fill(ax, [x(ok); flipud(x(ok))], [lowerProfile(ok); flipud(upperProfile(ok))], color, 'FaceAlpha', alphaValue, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
