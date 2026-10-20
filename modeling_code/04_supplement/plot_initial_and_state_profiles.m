function paths = plot_initial_and_state_profiles( controlResult, treatmentResult, controlReplicates, treatmentReplicates, varargin)
    %PLOT_INITIAL_AND_STATE_PROFILES Initial-state and state-decomposition panels.

    % These panels use the same model results, experimental density records,
    % time windows, and annular normalization as the main density figures.

    style = figure_style();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'controlResult', @isstruct);
    addRequired(p, 'treatmentResult', @isstruct);
    addRequired(p, 'controlReplicates', @isstruct);
    addRequired(p, 'treatmentReplicates', @isstruct);
    addParameter(p, 'OutputDir', fullfile(pwd, 'supplementary_figure_panels'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'Panels', {'initial_conditions', 'state_decomposition'}, @(x) iscellstr(x) && all(ismember(x, {'initial_conditions', 'state_decomposition'}))); %#ok<ISCLSTR>
    addParameter(p, 'InitialTimeWindowHours', analysis.initialDensityWindowHours, @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
    addParameter(p, 'DensityIntervalsElapsedHours', analysis.densityIntervalsHours - analysis.experimentStartHour, @(x) isnumeric(x) && size(x, 2) == 2 && all(x(:, 2) > x(:, 1)));
    addParameter(p, 'EdgeRadius', edgeRadius, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    parse(p, controlResult, treatmentResult, controlReplicates, treatmentReplicates, varargin{:});
    opts = p.Results;

    visibility = 'off';
    if logical(opts.Visible), visibility = 'on'; end
    outputDir = char(string(opts.OutputDir));
    if ~isfolder(outputDir), mkdir(outputDir); end

    panels = p.Results.Panels;
    if ismember('initial_conditions', panels)
        fig = local_initial_conditions(controlResult, treatmentResult, controlReplicates, treatmentReplicates, visibility, double(opts.InitialTimeWindowHours), double(opts.EdgeRadius));
        paths.initialConditions = local_export(fig, outputDir, 'initial_condition_profiles_24_28h', double(opts.Resolution));
        if ~logical(opts.Visible), close(fig); end
    end
    if ismember('state_decomposition', panels)
        fig = local_state_decomposition(controlResult, visibility, double(opts.DensityIntervalsElapsedHours));
        paths.stateDecomposition = local_export(fig, outputDir, 'control_upper_basal_density_timecourse', double(opts.Resolution));
        if ~logical(opts.Visible), close(fig); end
    end
end

function fig = local_initial_conditions(control, treatment, controlReps, treatmentReps, visibility, timeWindow, edgeRadius)
    style = figure_style();
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.initialConditionFigureSize]);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    local_initial_panel(nexttile(layout), control, controlReps, 'Control', timeWindow, edgeRadius, true);
    local_initial_panel(nexttile(layout), treatment, treatmentReps, 'shAPLNR', timeWindow, edgeRadius, false);
    xlabel(layout, 'Distance from edge (\mum)', 'FontName', style.fontName, 'FontSize', style.supplementTickFontSize, 'Color', 'k');
end

function local_initial_panel(ax, model, replicates, titleText, timeWindow, edgeRadius, showYLabel)
    style = figure_style();
    [x, experimental] = local_experimental_initial(replicates, timeWindow, edgeRadius);
    r = edgeRadius - x;
    fitProfile = exp(-0.5 .* ((r - double(model.params.MuM)) ./ double(model.params.SigmaM)).^2);
    fitMass = trapz(r(end:-1:1), 2 * pi .* r(end:-1:1) .* fitProfile(end:-1:1));
    fitProfile = fitProfile ./ max(fitMass, eps);
    bar(ax, x, experimental, 1, 'FaceColor', style.initialHistogramColor, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x, fitProfile, '-', 'Color', style.initialConditionColor, 'LineWidth', style.initialConditionLineWidth);
    hold(ax, 'off');
    xlim(ax, [0 100]);
    ylim(ax, [0 2e-5]);
    pbaspect(ax, style.initialConditionAspectRatio);
    xticks(ax, [0 50 100]);
    yticks(ax, [0 1e-5 2e-5]);
    ax.YAxis.Exponent = 0;
    if showYLabel
        ax.YTickLabel = {'0', '1', '2'};
        ylabel(ax, 'Areal density (10^{-5}/\mum^2)');
    else
        ax.YTickLabel = {};
    end
    title(ax, titleText, 'FontWeight', 'normal', 'Color', 'k');
    local_style(ax);
end

function [x, profile] = local_experimental_initial(replicates, timeWindow, edgeRadius)
    profiles = [];
    x = [];
    for i = 1:numel(replicates)
        [r, t, density] = experimental_areal_density_profile(replicates(i).density);
        absoluteTime = t;
        if max(absoluteTime, [], 'omitnan') <= 31
            absoluteTime = absoluteTime + 24;
        end
        keep = absoluteTime >= timeWindow(1) & absoluteTime <= timeWindow(2);
        replicateProfile = mean(density(keep, :), 1, 'omitnan').';
        xRep = edgeRadius - r(:);
        [xRep, order] = sort(xRep, 'ascend');
        replicateProfile = replicateProfile(order);
        if isempty(x)
            x = xRep;
            profiles = nan(numel(x), numel(replicates));
        end
        profiles(:, i) = interp1(xRep, replicateProfile, x, 'linear', NaN); %#ok<AGROW>
    end
    profile = mean(profiles, 2, 'omitnan');
end

function fig = local_state_decomposition(model, visibility, intervals)
    style = figure_style();
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.stateFigureSize]);
    ax = axes(fig);
    [x, order] = sort(double(model.params.R) - model.rs(:), 'ascend');
    upper = local_orient(model.m_upper, model);
    basal = local_orient(model.m_basal, model);
    lineStyles = style.stateLineStyles;
    stateColors = style.stateColors;
    stateLabels = {'Total', 'Upper', 'Basal'};
    hold(ax, 'on');
    for i = 1:size(intervals, 1)
        keep = model.ts(:) >= intervals(i, 1) & model.ts(:) <= intervals(i, 2);
        normalizedUpper = nan(nnz(keep), numel(model.rs));
        normalizedBasal = normalizedUpper;
        rows = find(keep);
        for j = 1:numel(rows)
            total = max(upper(rows(j), :), 0) + max(basal(rows(j), :), 0);
            mass = trapz(model.rs(:), 2 * pi .* model.rs(:) .* total(:));
            normalizedUpper(j, :) = max(upper(rows(j), :), 0) ./ max(mass, eps);
            normalizedBasal(j, :) = max(basal(rows(j), :), 0) ./ max(mass, eps);
        end
        stateProfiles = [mean(normalizedUpper + normalizedBasal, 1, 'omitnan'); mean(normalizedUpper, 1, 'omitnan'); mean(normalizedBasal, 1, 'omitnan')];
        for state = 1:3
            displayName = '';
            if i == size(intervals, 1), displayName = stateLabels{state}; end
            plot(ax, x, stateProfiles(state, order), lineStyles{i}, 'Color', stateColors(state, :), 'LineWidth', style.stateLineWidth, 'DisplayName', displayName, 'HandleVisibility', local_handle_visibility(i == size(intervals, 1)));
        end
    end
    hold(ax, 'off');
    xlim(ax, [0 350]);
    ylim(ax, [0 1.2e-5]);
    xticks(ax, [0 175 350]);
    yticks(ax, [0 0.6e-5 1.2e-5]);
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Areal density (1/\mum^2)');
    local_style(ax);
    lgd = legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', style.legendFontSize);
    lgd.TextColor = 'k';
end

function value = local_handle_visibility(show)
    value = 'off';
    if show, value = 'on'; end
end

function values = local_orient(values, model)
    assert(isequal(size(values), [numel(model.ts), numel(model.rs)]), 'plot_initial_and_state_profiles:BadSeriesSize', 'Model state arrays must be time-by-radius.');
end

function local_style(ax)
    style = figure_style();
    set(ax, 'FontName', style.fontName, 'FontSize', style.supplementTickFontSize, 'LineWidth', style.supplementAxesLineWidth, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'TickDir', 'in', 'Layer', 'top', 'Color', 'w');
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
end

function paths = local_export(fig, outputDir, baseName, resolution)
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    paths.fig = fullfile(outputDir, [baseName '.fig']);
    exportgraphics(fig, paths.png, 'Resolution', resolution, 'BackgroundColor', 'white');
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, paths.fig);
end
