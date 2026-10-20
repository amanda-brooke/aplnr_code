function paths = plot_experimental_speed_maps(controlReplicates, treatmentReplicates, varargin)
    %PLOT_EXPERIMENTAL_SPEED_MAPS Plot displacement-count-weighted radial speed maps.

    style = figure_style();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'controlReplicates', @isstruct);
    addRequired(p, 'treatmentReplicates', @isstruct);
    addParameter(p, 'OutputDir', fullfile(pwd, 'offset24_figures'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'AbsoluteTimeLimits', analysis.timeWindowHours, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'SpeedLimits', style.speedLimits, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'EdgeRadius', edgeRadius, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'RadialBinCount', analysis.radialBinCount, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && mod(x, 1) == 0);
    parse(p, controlReplicates, treatmentReplicates, varargin{:});
    opts = p.Results;

    visibility = 'off';
    if logical(opts.Visible)
        visibility = 'on';
    end
    [fig, baseName] = local_radial_speed_maps(controlReplicates, treatmentReplicates, visibility, double(opts.AbsoluteTimeLimits), double(opts.SpeedLimits), double(opts.EdgeRadius), double(opts.RadialBinCount));

    paths = local_export(fig, char(string(opts.OutputDir)), baseName, double(opts.Resolution));
    if ~logical(opts.Visible)
        close(fig);
    else
        figure(fig);
    end
end

function [fig, baseName] = local_radial_speed_maps( control, treatment, visibility, absoluteTimeLimits, speedLimits, edgeRadius, radialBinCount)
    style = figure_style();
    fig = local_figure(visibility, style.speedMapFigureSize);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(layout);
    local_plot_speed_map(ax, control, false, absoluteTimeLimits, speedLimits, edgeRadius, radialBinCount);
    ax = nexttile(layout);
    local_plot_speed_map(ax, treatment, true, absoluteTimeLimits, speedLimits, edgeRadius, radialBinCount);
    baseName = 'figure_experimental_radial_speed_maps';
end

function local_plot_speed_map(ax, replicates, addColorbar, absoluteTimeLimits, speedLimits, edgeRadius, radialBinCount)
    [x, elapsed, speed] = local_mean_radial_speed_map(replicates, edgeRadius, radialBinCount);
    imageHandle = imagesc(ax, x, elapsed + 24, speed.');
    imageHandle.AlphaData = isfinite(imageHandle.CData);
    set(ax, 'YDir', 'reverse');
    clim(ax, speedLimits);
    colormap(ax, flipud(magma_colormap(256)));
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Hours post BMP4');
    local_style_map_axis(ax, absoluteTimeLimits, edgeRadius);
    if addColorbar
        cb = colorbar(ax);
        cb.Ticks = [speedLimits(1) mean(speedLimits) speedLimits(2)];
        ylabel(cb, 'Radial speed (\mum/h)');
        local_style_colorbar(cb);
    end
end

function [xGrid, elapsedGrid, meanSpeed] = local_mean_radial_speed_map(replicates, edgeRadius, radialBinCount)
    xGrid = linspace(0, edgeRadius, radialBinCount).';
    first = replicates(1).radialVelocity;
    elapsedGrid = first.tCenters(:) - min(replicates(1).density.tHours(:));
    sumSpeed = zeros(numel(xGrid), numel(elapsedGrid));
    sumCounts = zeros(numel(xGrid), numel(elapsedGrid));
    for i = 1:numel(replicates)
        summary = replicates(i).radialVelocity;
        values = summary.avgV;
        counts = summary.counts;
        expectedSize = [numel(summary.distanceFromEdgeCenters), numel(summary.tCenters)];
        assert(isequal(size(values), expectedSize) && isequal(size(counts), expectedSize), 'plot_experimental_speed_maps:BadMapSize', 'Speed and count maps must be distance-by-time.');
        elapsed = summary.tCenters(:) - min(replicates(i).density.tHours(:));
        xReplicate = summary.distanceFromEdgeCenters(:);
        [xReplicate, order] = sort(xReplicate, 'ascend');
        values = values(order, :);
        counts = counts(order, :);
        [timeQuery, xQuery] = meshgrid(elapsedGrid, xGrid);
        valueQuery = interp2(elapsed, xReplicate, values, timeQuery, xQuery, 'linear', NaN);
        countQuery = interp2(elapsed, xReplicate, counts, timeQuery, xQuery, 'linear', NaN);
        ok = isfinite(valueQuery) & isfinite(countQuery) & countQuery > 0;
        weights = countQuery;
        weights(~ok) = 0;
        sumSpeed(ok) = sumSpeed(ok) + valueQuery(ok) .* weights(ok);
        sumCounts(ok) = sumCounts(ok) + weights(ok);
    end
    meanSpeed = sumSpeed ./ max(sumCounts, 1);
    meanSpeed(sumCounts == 0) = NaN;
end

function fig = local_figure(visibility, sizeInches)
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 sizeInches], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
end

function local_style_map_axis(ax, absoluteTimeLimits, edgeRadius)
    style = figure_style();
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'Color', 'w', 'LineWidth', style.axesLineWidth, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Layer', 'top', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'XLim', [0 edgeRadius], 'YLim', absoluteTimeLimits, 'XTick', [0 edgeRadius / 2 edgeRadius], 'YTick', [absoluteTimeLimits(1) mean(absoluteTimeLimits) absoluteTimeLimits(2)]);
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
end

function local_style_colorbar(cb)
    style = figure_style();
    set(cb, 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'LineWidth', style.colorbarLineWidth, 'Color', 'k');
    cb.Label.FontName = style.fontName;
    cb.Label.FontSize = style.colorbarLabelFontSize;
    cb.Label.Color = 'k';
end

function paths = local_export(fig, outputDir, baseName, resolution)
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    paths = struct();
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    paths.fig = fullfile(outputDir, [baseName '.fig']);
    exportgraphics(fig, paths.png, 'Resolution', resolution, 'BackgroundColor', 'white');
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, paths.fig);
end
