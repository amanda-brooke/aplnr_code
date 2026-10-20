function paths = plot_parameter_sweeps(sweeps, varargin)
    %PLOT_PARAMETER_SWEEPS Render sweep figures.

    style = figure_style();
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'sweeps', @(x) isstruct(x) || ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'DensityLimits', style.densityLimits, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'Panels', {'sourceLambda', 'betaDiffusivity', 'alphaGamma'}, @iscellstr);
    parse(p, sweeps, varargin{:});
    opts = p.Results;
    assert(all(ismember(opts.Panels, {'sourceLambda', 'betaDiffusivity', 'alphaGamma'})), 'Unknown supplementary sweep panel.');

    if ~isstruct(sweeps)
        loaded = load(char(string(sweeps)), 'out');
        sweeps = loaded.out;
    end
    outputDir = char(string(opts.OutputDir));
    if isempty(outputDir)
        outputDir = sweeps.outputDir;
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    visibility = 'off';
    if logical(opts.Visible)
        visibility = 'on';
    end

    paths = struct();
    modelTimeWindowHours = sweeps.config.modelTimeWindowHours;
    sourceExamples = local_corner_examples(sweeps.sourceLambda, sweeps.base);
    betaExamples = local_corner_examples(sweeps.betaDiffusivity, sweeps.base);
    alphaGammaExamples = local_corner_examples(sweeps.alphaGamma, sweeps.base);
    % At alpha=1, upper and basal dynamics are identical, so gamma does not
    % affect total density. Keep the low-gamma example and omit its duplicate.
    if abs(max(sweeps.alphaGamma.xValues) - 1) < 10 * eps
        alphaGammaExamples(4) = [];
    end
    isBaselineExample = arrayfun(@(example) abs(example.x - sweeps.alphaGamma.baselineX) < 10 * eps && abs(example.y - sweeps.alphaGamma.baselineY) < 10 * eps, alphaGammaExamples);
    alphaGammaExamples(isBaselineExample) = [];
    if ismember('sourceLambda', opts.Panels)
        sourceExamples = local_run_examples(sourceExamples);
    end
    if ismember('betaDiffusivity', opts.Panels)
        betaExamples = local_run_examples(betaExamples);
    end
    if ismember('alphaGamma', opts.Panels)
        alphaGammaExamples = local_run_examples(alphaGammaExamples);
    end
    baselineModel = solve_chemotaxis_model(sweeps.base);
    [~, ~, baselineLigand] = local_corner_profiles( baselineModel, modelTimeWindowHours);
    ligandNormalization = max(baselineLigand, [], 'omitnan');
    if ~(isfinite(ligandNormalization) && ligandNormalization > 0)
        ligandNormalization = 1;
    end
    commonColorLimits = local_common_color_limits( {sweeps.sourceLambda, sweeps.betaDiffusivity, sweeps.alphaGamma});
    paths.sharedColorLimits = commonColorLimits;
    names = {'sourceLambda', 'betaDiffusivity', 'alphaGamma'};
    exampleSets = {sourceExamples, betaExamples, alphaGammaExamples};
    for i = 1:numel(names)
        name = names{i};
        if ismember(name, opts.Panels)
            paths.(name) = local_render_sweep(sweeps.(name), exampleSets{i}, baselineModel, ligandNormalization, outputDir, visibility, double(opts.Resolution), commonColorLimits, double(opts.DensityLimits), modelTimeWindowHours);
        end
    end
end

function paths = local_render_sweep(sweep, examples, baselineModel, ligandNormalization, outputDir, visibility, resolution, colorLimits, densityLimits, modelTimeWindowHours)
    paths = struct();

    fig = local_composite_figure(sweep, examples, baselineModel, ligandNormalization, visibility, colorLimits, densityLimits, modelTimeWindowHours);
    paths.profiles = local_export(fig, outputDir, [sweep.fileBase '_profiles'], resolution);
    local_finish_figure(fig, visibility);
end

function fig = local_composite_figure(sweep, examples, baselineModel, ligandNormalization, visibility, colorLimits, densityLimits, modelTimeWindowHours)
    style = figure_style();
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.sweepFigureSize], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
    heatmapPosition = [0.10 0.30 0.22 0.47];
    ax = axes(fig, 'Units', 'normalized', 'Position', heatmapPosition);
    local_plot_heatmap(ax, sweep, colorLimits);
    hold(ax, 'on');
    [dx, dy] = local_cell_half_widths(sweep);
    for i = 1:numel(examples)
        boxHandle = rectangle(ax, 'Position', [examples(i).x-dx, examples(i).y-dy, 2*dx, 2*dy], 'EdgeColor', examples(i).color, 'LineWidth', style.sweepBoxLineWidth, 'FaceColor', 'none', 'Clipping', 'off');
        uistack(boxHandle, 'top');
    end
    set(ax, 'Layer', 'bottom');
    hold(ax, 'off');
    cb = colorbar(ax, 'southoutside');
    local_style_horizontal_colorbar(cb, [0.10 0.08 0.22 0.035]);
    set(ax, 'Position', heatmapPosition);

    axDensity = axes(fig, 'Units', 'normalized', 'Position', [0.41 0.30 0.22 0.47]);
    axLigand = axes(fig, 'Units', 'normalized', 'Position', [0.71 0.30 0.22 0.47]);
    local_plot_overlay_profiles(axDensity, axLigand, examples, baselineModel, ligandNormalization, densityLimits, modelTimeWindowHours);
end

function local_style_horizontal_colorbar(cb, position)
    style = figure_style();
    set(cb, 'Position', position, 'Color', 'k', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'LineWidth', style.colorbarLineWidth, 'AxisLocation', 'out');
    cb.Label.String = '';
end

function local_plot_overlay_profiles(axDensity, axLigand, examples, baselineModel, ligandNormalization, densityLimits, modelTimeWindowHours)
    style = figure_style();
    hold(axDensity, 'on');
    hold(axLigand, 'on');
    for i = 1:numel(examples)
        [x, densityProfile, ligandProfile] = local_corner_profiles(examples(i).model, modelTimeWindowHours);
        plot(axDensity, x, densityProfile, '-', 'Color', examples(i).color, 'LineWidth', style.sweepDensityLineWidth);
        plot(axLigand, x, ligandProfile ./ ligandNormalization, '-', 'Color', examples(i).color, 'LineWidth', style.sweepLigandLineWidth);
    end
    [x, densityProfile, ligandProfile] = local_corner_profiles( baselineModel, modelTimeWindowHours);
    plot(axDensity, x, densityProfile, '-', 'Color', 'k', 'LineWidth', style.sweepBaselineDensityLineWidth);
    local_style_density_overlay_axis(axDensity, densityLimits);
    plot(axLigand, x, ligandProfile ./ ligandNormalization, '-', 'Color', 'k', 'LineWidth', style.sweepBaselineLigandLineWidth);
    local_style_ligand_overlay_axis(axLigand);
    hold(axDensity, 'off');
    hold(axLigand, 'off');
end

function local_style_density_overlay_axis(ax, densityLimits)
    local_style_overlay_axis(ax, densityLimits);
    set(ax, 'XTick', [0 175 350], 'YTick', [densityLimits(1) mean(densityLimits) densityLimits(2)]);
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Areal density (1/\mum^2)');
end

function local_style_ligand_overlay_axis(ax)
    local_style_overlay_axis(ax, [0 1.2]);
    set(ax, 'XTick', [0 175 350], 'YTick', [0 0.4 0.8 1.2]);
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Relative ligand concentration');
end

function local_style_overlay_axis(ax, yLimits)
    style = figure_style();
    set(ax, 'XLim', [0 350], 'YLim', yLimits, 'XColor', 'k', 'YColor', 'k', 'Color', 'w', 'LineWidth', style.axesLineWidth, 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Layer', 'top', 'FontName', style.fontName, 'FontSize', style.tickFontSize);
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
    pbaspect(ax, [1 1 1]);
    grid(ax, 'off');
    box(ax, 'on');
end

function local_plot_heatmap(ax, sweep, colorLimits)
    style = figure_style();
    imagesc(ax, sweep.xValues, sweep.yValues, sweep.relativeRmse);
    set(ax, 'YDir', 'normal');
    axis(ax, 'tight');
    colormap(ax, flipud(gray(256)));
    if isempty(colorLimits)
        finiteValues = sweep.relativeRmse(isfinite(sweep.relativeRmse));
        if ~isempty(finiteValues)
            colorLimits = [min(finiteValues), local_percentile(finiteValues, 95)];
        end
    end
    if ~isempty(colorLimits)
        if colorLimits(2) <= colorLimits(1)
            colorLimits(2) = colorLimits(1) + eps(colorLimits(1));
        end
        clim(ax, colorLimits);
    end
    [xLabel, yLabel] = local_compact_axis_labels(sweep);
    xlabel(ax, xLabel, 'Color', 'k');
    ylabel(ax, yLabel, 'Color', 'k');
    set(ax, 'XTick', local_reduced_ticks(sweep.xValues), 'YTick', local_reduced_ticks(sweep.yValues), 'XColor', 'k', 'YColor', 'k', 'Color', 'w', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'LineWidth', style.axesLineWidth, 'TickDir', 'in', 'Layer', 'top');
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
    if strcmp(sweep.id, 'alpha_gamma')
        set(ax, 'YTick', [0 0.025 0.05], 'YTickLabel', {'0', '2.5', '5'});
        text(ax, -0.02, 1.03, '\times10^{-2}', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'Clipping', 'off');
    end
    xtickangle(ax, 0);
    pbaspect(ax, [1 1 1]);
    grid(ax, 'off');
    box(ax, 'on');
    hold(ax, 'on');
    [dx, dy] = local_cell_half_widths(sweep);
    [~, baselineX] = min(abs(sweep.xValues - sweep.baselineX));
    [~, baselineY] = min(abs(sweep.yValues - sweep.baselineY));
    rectangle(ax, 'Position', [sweep.xValues(baselineX)-dx, sweep.yValues(baselineY)-dy, 2*dx, 2*dy], 'EdgeColor', 'k', 'LineWidth', style.sweepBoxLineWidth, 'FaceColor', 'none');
    hold(ax, 'off');
end

function [xLabel, yLabel] = local_compact_axis_labels(sweep)
    switch sweep.id
        case 'source_radius_lambda'
            xLabel = 'r_s (\mum)';
            yLabel = '\lambda (\mum)';
        case 'beta_mesoderm_diffusivity'
            xLabel = '\beta (\mum^2/h)';
            yLabel = 'D_m (\mum^2/h)';
        case 'alpha_gamma'
            xLabel = '\alpha';
            yLabel = '\gamma (h^{-1})';
        otherwise
            xLabel = sweep.xLabel;
            yLabel = sweep.yLabel;
    end
end

function limits = local_common_color_limits(sweeps)
    finiteValues = [];
    for i = 1:numel(sweeps)
        values = sweeps{i}.relativeRmse;
        finiteValues = [finiteValues; values(isfinite(values))]; %#ok<AGROW>
    end
    if isempty(finiteValues)
        limits = [];
        return;
    end
    limits = [min(finiteValues), local_percentile(finiteValues, 95)];
    if limits(2) <= limits(1)
        limits(2) = limits(1) + eps(limits(1));
    end
end

function examples = local_corner_examples(sweep, base)
    style = figure_style();
    xLow = min(sweep.xValues);
    xHigh = max(sweep.xValues);
    yLow = min(sweep.yValues);
    yHigh = max(sweep.yValues);
    values = [xLow yLow; xHigh yLow; xLow yHigh; xHigh yHigh];
    colors = style.sweepCornerColors;
    examples = repmat(struct('x', NaN, 'y', NaN, 'color', [], 'params', struct(), 'model', struct()), 4, 1);
    for i = 1:4
        params = base;
        switch sweep.id
            case 'source_radius_lambda'
                params.SourceRadius = values(i, 1);
                params.Da = values(i, 2).^2 .* params.K;
            case 'beta_mesoderm_diffusivity'
                params.Beta = values(i, 1);
                params.Dm = values(i, 2);
            case 'alpha_gamma'
                params.Alpha = values(i, 1);
                params.Gamma = values(i, 2);
            otherwise
                error('plot_parameter_sweeps:UnknownSweep', 'Unknown sweep id "%s".', sweep.id);
        end
        examples(i) = struct('x', values(i, 1), 'y', values(i, 2), 'color', colors(i, :), 'params', params, 'model', struct());
    end
end

function examples = local_run_examples(examples)
    for i = 1:numel(examples)
        examples(i).model = solve_chemotaxis_model(examples(i).params);
    end
end

function [x, densityProfile, ligandProfile] = local_corner_profiles( model, modelTimeWindowHours)
    [~, t, density] = model_areal_density_profile(model);
    x = double(model.params.R) - model.rs(:);
    [x, order] = sort(x, 'ascend');
    keep = t >= modelTimeWindowHours(1) & t <= modelTimeWindowHours(2);
    densityProfile = mean(density(keep, order), 1, 'omitnan').';
    ligand = local_orient(model.a, numel(model.ts), numel(model.rs));
    ligandProfile = mean(ligand(keep, order), 1, 'omitnan').';
end

function values = local_orient(values, nT, nR)
    if isequal(size(values), [nT nR])
        return;
    end
    error('plot_parameter_sweeps:BadSeriesSize', 'Expected time-by-radius size [%d %d], got [%s].', nT, nR, num2str(size(values)));
end

function ticks = local_reduced_ticks(values)
    indices = unique(round(linspace(1, numel(values), min(5, numel(values)))));
    ticks = values(indices);
end

function [dx, dy] = local_cell_half_widths(sweep)
    dx = local_cell_half_width(sweep.xValues);
    dy = local_cell_half_width(sweep.yValues);
end

function halfWidth = local_cell_half_width(values)
    if numel(values) > 1
        halfWidth = min(diff(sort(unique(values)))) / 2;
    else
        halfWidth = max(0.5, abs(values(1)) * 0.05);
    end
end

function value = local_percentile(values, percentile)
    values = sort(values(:));
    position = 1 + (numel(values) - 1) * percentile / 100;
    lowerIndex = floor(position);
    upperIndex = ceil(position);
    if lowerIndex == upperIndex
        value = values(lowerIndex);
    else
        fraction = position - lowerIndex;
        value = values(lowerIndex) * (1 - fraction) + values(upperIndex) * fraction;
    end
end

function paths = local_export(fig, outputDir, baseName, resolution)
    paths = struct();
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    paths.fig = fullfile(outputDir, [baseName '.fig']);
    figureSize = fig.Position(3:4);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureSize], 'PaperSize', figureSize);
    print(fig, paths.png, '-dpng', sprintf('-r%d', round(resolution)));
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, paths.fig);
end

function local_finish_figure(fig, visibility)
    if strcmp(visibility, 'off')
        close(fig);
    else
        figure(fig);
    end
end
