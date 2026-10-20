function paths = plot_boundary_sensitivity(data, varargin)
    %PLOT_BOUNDARY_SENSITIVITY Plot boundary endpoints.

    analysis = analysis_settings();
    style = figure_style();
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'data', @(x) isstruct(x) || ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'DensityIntervalsElapsedHours', analysis.densityIntervalsHours - analysis.experimentStartHour, @(x) isnumeric(x) && size(x, 2) == 2 && all(isfinite(x), 'all') && all(x(:, 2) > x(:, 1)) && size(x, 1) <= 5);
    addParameter(p, 'DensityLimits', style.densityLimits, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    parse(p, data, varargin{:});
    opts = p.Results;

    if ~isstruct(data)
        loaded = load(char(string(data)), 'out');
        data = loaded.out;
    end
    outputDir = char(string(opts.OutputDir));
    if isempty(outputDir)
        outputDir = data.outputDir;
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    visibility = 'off';
    if logical(opts.Visible)
        visibility = 'on';
    end

    baselineParams = data.base;
    baselineParams.LigandBoundaryAbsorption = 1;
    baselineParams.MesodermBoundaryAbsorption = 1;
    ligandNoFluxParams = baselineParams;
    ligandNoFluxParams.LigandBoundaryAbsorption = 0;
    mesodermNoFluxParams = baselineParams;
    mesodermNoFluxParams.MesodermBoundaryAbsorption = 0;

    baseline = local_profiles(solve_chemotaxis_model(baselineParams));
    ligandNoFlux = local_profiles(solve_chemotaxis_model(ligandNoFluxParams));
    mesodermNoFlux = local_profiles(solve_chemotaxis_model(mesodermNoFluxParams));

    baselineId = 'publication_model';
    if isfield(data, 'baselineId')
        baselineId = char(string(data.baselineId));
    end
    intervals = double(opts.DensityIntervalsElapsedHours);
    densityLimits = double(opts.DensityLimits);
    ligandLimit = local_ligand_limit(baseline, ligandNoFlux, mesodermNoFlux);
    baseName = [baselineId '_boundary_endpoint_profiles'];
    fig = local_boundary_figure(baseline, ligandNoFlux, mesodermNoFlux, visibility, intervals, densityLimits, ligandLimit);
    paths = local_export(fig, outputDir, baseName, double(opts.Resolution));
    local_finish_figure(fig, visibility);

end

function fig = local_boundary_figure(baseline, ligandNoFlux, mesodermNoFlux, visibility, intervals, densityLimits, ligandLimit)
    style = figure_style();
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.boundaryFigureSize], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
    layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
    legendHandle = local_plot_overlay(nexttile(layout), baseline, ligandNoFlux, 'APELA outer boundary', true, true, true, false, true, intervals, densityLimits, ligandLimit);
    local_plot_overlay(nexttile(layout), baseline, mesodermNoFlux, 'Mesoderm outer boundary', false, true, true, true, true, intervals, densityLimits, ligandLimit);
    legendHandle.Layout.Tile = 'north';
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

function profiles = local_profiles(model)
    [~, t, density] = model_areal_density_profile(model);
    x = double(model.params.R) - model.rs(:);
    [x, order] = sort(x, 'ascend');
    ligand = local_orient(model.a, numel(model.ts), numel(model.rs));
    profiles = struct( 'x', x, 't', t(:), 'ligand', ligand(end, order).', 'density', density(:, order));
end

function legendHandle = local_plot_overlay(ax, baseline, noFlux, panelTitle, showLegend, showLeftAxis, showRightAxis, showXAxis, showTitle, intervals, densityLimits, ligandLimit)
    style = figure_style();
    legendHandle = gobjects(0);
    labels = arrayfun(@(i) sprintf('%g-%g h', intervals(i, 1) + 24, intervals(i, 2) + 24), 1:size(intervals, 1), 'UniformOutput', false);
    colors = local_time_colors(size(intervals, 1));

    yyaxis(ax, 'left');
    hold(ax, 'on');
    for i = 1:size(intervals, 1)
        baselineProfile = local_interval_mean( baseline.t, baseline.density, intervals(i, :));
        noFluxProfile = local_interval_mean( noFlux.t, noFlux.density, intervals(i, :));
        plot(ax, baseline.x, baselineProfile, '-', 'Color', colors(i, :), 'LineWidth', style.boundaryDensityLineWidth, 'HandleVisibility', 'off');
        plot(ax, noFlux.x, noFluxProfile, ':', 'Color', colors(i, :), 'LineWidth', style.boundaryDensityLineWidth, 'HandleVisibility', 'off');
    end
    if showLeftAxis
        ylabel(ax, 'Areal density (10^{-5}/\mum^2)', 'Color', 'k', 'FontSize', style.boundaryLabelFontSize);
    else
        ylabel(ax, '');
    end
    ylim(ax, densityLimits);

    yyaxis(ax, 'right');
    plot(ax, baseline.x, baseline.ligand, '-', 'Color', style.boundaryLigandColor, 'LineWidth', style.boundaryLigandLineWidth, 'HandleVisibility', 'off');
    plot(ax, noFlux.x, noFlux.ligand, ':', 'Color', style.boundaryLigandColor, 'LineWidth', style.boundaryLigandLineWidth, 'HandleVisibility', 'off');
    if showRightAxis
        ylabel(ax, 'Ligand concentration (a.u.)', 'Color', 'k', 'FontSize', style.boundaryLabelFontSize);
    else
        ylabel(ax, '');
    end
    ylim(ax, [0 ligandLimit]);

    if showXAxis
        xlabel(ax, 'Distance from edge (\mum)', 'Color', 'k', 'FontSize', style.boundaryLabelFontSize);
    else
        xlabel(ax, '');
    end
    if showTitle
        title(ax, panelTitle, 'FontWeight', 'normal', 'Color', 'k', 'FontSize', style.boundaryLabelFontSize);
    end
    xlim(ax, [0 350]);
    xticks(ax, [0 175 350]);
    local_style_axis(ax, showLeftAxis, showRightAxis, densityLimits, ligandLimit);
    if ~showXAxis
        set(ax, 'XTickLabel', {});
    end
    pbaspect(ax, [1 1 1]);

    if showLegend
        yyaxis(ax, 'left');
        legendHandles = gobjects(size(intervals, 1) + 3, 1);
        for i = 1:size(intervals, 1)
            legendHandles(i) = plot(ax, nan, nan, '-', 'Color', colors(i, :), 'LineWidth', style.boundaryDensityLineWidth);
        end
        nIntervals = size(intervals, 1);
        legendHandles(nIntervals + 1) = plot(ax, nan, nan, 'k-', 'LineWidth', style.boundaryDensityLineWidth);
        legendHandles(nIntervals + 2) = plot(ax, nan, nan, 'k:', 'LineWidth', style.boundaryDensityLineWidth);
        legendHandles(nIntervals + 3) = plot(ax, nan, nan, '-', 'Color', style.boundaryLigandLegendColor, 'LineWidth', style.boundaryLigandLineWidth);
        legendHandle = legend(ax, legendHandles, [labels {'Dirichlet', 'No flux', 'Ligand'}], 'NumColumns', 2, 'Box', 'off', 'TextColor', 'k', 'FontName', style.fontName, 'FontSize', style.boundaryLegendFontSize);
    end
    hold(ax, 'off');
end

function profile = local_interval_mean(t, values, interval)
    keep = t >= interval(1) & t <= interval(2);
    if any(keep)
        profile = mean(values(keep, :), 1, 'omitnan').';
    else
        profile = interp1(t, values, mean(interval), 'linear', 'extrap').';
    end
end

function colors = local_time_colors(nIntervals)
    fiveColors = [
        0.42 0.18 0.68
        0.20 0.42 0.78
        0.92 0.22 0.22
        0.95 0.56 0.20
        0.00 0.00 0.00
        ];
    if nIntervals == 3
        colors = fiveColors([1 3 4], :);
    else
        colors = fiveColors(1:nIntervals, :);
    end
end

function local_style_axis(ax, showLeftAxis, showRightAxis, densityLimits, ligandLimit)
    style = figure_style();
    set(ax, 'XColor', 'k', 'YColor', 'k', 'Color', 'w', 'FontName', style.fontName, 'FontSize', style.boundaryTickFontSize, 'LineWidth', style.supplementAxesLineWidth, 'TickDir', 'in', 'Layer', 'top');
    yyaxis(ax, 'left');
    if showLeftAxis
        ax.YColor = 'k';
        densityTicks = [densityLimits(1) mean(densityLimits) densityLimits(2)];
        ax.YAxis(1).Exponent = 0;
        set(ax, 'YTick', densityTicks, 'YTickLabel', compose('%g', densityTicks ./ 1e-5));
    else
        ax.YColor = 'k';
        set(ax, 'YTick', [], 'YTickLabel', {});
    end
    yyaxis(ax, 'right');
    if showRightAxis
        ax.YColor = 'k';
        set(ax, 'YTick', [0 ligandLimit / 2 ligandLimit]);
    else
        ax.YColor = 'k';
        set(ax, 'YTick', [], 'YTickLabel', {});
    end
    grid(ax, 'off');
    box(ax, 'on');
end

function limit = local_ligand_limit(varargin)
    maximum = 0;
    for i = 1:nargin
        maximum = max(maximum, max(varargin{i}.ligand, [], 'omitnan'));
    end
    if ~(isfinite(maximum) && maximum > 0)
        limit = 1;
        return;
    end
    limit = ceil(10 * 1.05 * maximum) / 10;
    if limit <= maximum
        limit = 1.1 * maximum;
    end
end

function values = local_orient(values, nT, nR)
    if isequal(size(values), [nT nR])
        return;
    end
    error('plot_boundary_sensitivity:BadSeriesSize', 'Expected time-by-radius size [%d %d], got [%s].', nT, nR, num2str(size(values)));
end
