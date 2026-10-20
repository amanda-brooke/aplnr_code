function paths = plot_model_figure(kind, controlResult, treatmentResult, varargin)
    %PLOT_MODEL_FIGURE Build one canonical model-output figure.

    style = figure_style();
    analysis = analysis_settings();
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'kind', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addRequired(p, 'controlResult', @isstruct);
    addRequired(p, 'treatmentResult', @(x) isempty(x) || isstruct(x));
    addParameter(p, 'OutputDir', fullfile(pwd, 'minimal_chemotaxis_figures'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'DynamicsAbsoluteTimeLimits', analysis.timeWindowHours, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'SpeedLimits', style.speedLimits, @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(2) > x(1));
    addParameter(p, 'ApparentSpeedFrameIntervalHours', style.apparentSpeedFrameIntervalHours, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    parse(p, kind, controlResult, treatmentResult, varargin{:});
    opts = p.Results;

    kind = lower(strtrim(char(string(kind))));
    visibility = 'off';
    if logical(opts.Visible)
        visibility = 'on';
    end

    switch kind
        case 'ligand_drive'
            [fig, baseName] = local_ligand_drive(controlResult, visibility);
        case 'radial_speed_maps'
            local_require_treatment(treatmentResult, kind);
            [fig, baseName] = local_radial_speed_maps( controlResult, treatmentResult, visibility, double(opts.DynamicsAbsoluteTimeLimits), double(opts.SpeedLimits), double(opts.ApparentSpeedFrameIntervalHours));
        otherwise
            error('plot_model_figure:UnknownKind', 'Unknown figure kind "%s".', kind);
    end

    paths = local_export(fig, char(string(opts.OutputDir)), baseName, double(opts.Resolution));
    if ~logical(opts.Visible)
        close(fig);
    else
        figure(fig);
    end
end

function [fig, baseName] = local_ligand_drive(model, visibility)
    style = figure_style();
    fig = local_figure(visibility, style.ligandFigureSize);

    [x, order] = local_distance_from_edge(model);
    ligand = local_time_radius(model.a, model);
    ax = axes(fig);
    plot(ax, x, ligand(end, order), 'k-', 'LineWidth', style.profileLineWidth);
    hold(ax, 'on');
    sourceBoundary = double(model.params.R) - double(model.params.SourceRadius);
    xline(ax, sourceBoundary, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', style.referenceLineWidth);
    hold(ax, 'off');
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, '(a.u.)');
    ligandLimit = local_axis_ceiling(max(ligand(end, :), [], 'omitnan'));
    local_style_line_axis(ax, [0 350], [0 ligandLimit]);
    pbaspect(ax, [1 1 1]);
    baseName = 'figure_ligand_profile';
end

function [fig, baseName] = local_radial_speed_maps( control, treatment, visibility, absoluteTimeLimits, speedLimits, frameIntervalHours)
    style = figure_style();
    fig = local_figure(visibility, style.speedMapFigureSize);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(layout);
    local_plot_speed_map(ax, control, false, absoluteTimeLimits, speedLimits, frameIntervalHours);
    ax = nexttile(layout);
    local_plot_speed_map(ax, treatment, true, absoluteTimeLimits, speedLimits, frameIntervalHours);
    baseName = 'figure_model_radial_speed_maps';
end

function local_plot_speed_map(ax, model, addColorbar, absoluteTimeLimits, speedLimits, frameIntervalHours)
    [~, ~, density] = model_areal_density_profile(model);
    [x, order] = local_distance_from_edge(model);
    upper = max(local_time_radius(model.m_upper, model), 0);
    basal = max(local_time_radius(model.m_basal, model), 0);
    vUpper = local_time_radius(model.v_upper, model);
    vBasal = local_time_radius(model.v_basal, model);
    upperSpeed = model_apparent_radial_speed( vUpper, double(model.params.Dm), frameIntervalHours);
    basalSpeed = model_apparent_radial_speed( vBasal, double(model.params.Dm), frameIntervalHours);
    total = upper + basal;
    speed = (upper .* upperSpeed + basal .* basalSpeed) ./ max(total, eps);
    relativeDensity = density ./ max(density, [], 2, 'omitnan');
    speed(relativeDensity < 0.03) = NaN;

    imageHandle = imagesc(ax, x, model.ts(:) + 24, speed(:, order));
    imageHandle.AlphaData = isfinite(imageHandle.CData);
    set(ax, 'YDir', 'reverse');
    clim(ax, speedLimits);
    colormap(ax, flipud(magma_colormap(256)));
    xlabel(ax, 'Distance from edge (\mum)');
    ylabel(ax, 'Hours post BMP4');
    local_style_map_axis(ax, absoluteTimeLimits);
    if addColorbar
        cb = colorbar(ax);
        cb.Ticks = [speedLimits(1) mean(speedLimits) speedLimits(2)];
        ylabel(cb, 'Apparent radial speed (\mum/h)');
        local_style_colorbar(cb);
    end
end

function values = local_time_radius(values, model)
    nT = numel(model.ts);
    nR = numel(model.rs);
    if isequal(size(values), [nT nR])
        return;
    end
    error('plot_model_figure:BadSeriesSize', 'Expected time-by-radius size [%d %d], got [%s].', nT, nR, num2str(size(values)));
end

function [x, order] = local_distance_from_edge(model)
    x = double(model.params.R) - model.rs(:);
    [x, order] = sort(x, 'ascend');
end

function limit = local_axis_ceiling(value)
    if ~isfinite(value) || value <= 0
        limit = 1;
        return;
    end
    exponent = floor(log10(value));
    scale = 10 ^ exponent;
    mantissa = value ./ scale;
    candidates = [1 1.2 1.5 2 2.5 3 4 5 6 8 10];
    limit = candidates(find(candidates >= mantissa, 1, 'first')) .* scale;
end

function fig = local_figure(visibility, sizeInches)
    fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 sizeInches], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
end

function local_style_line_axis(ax, xLimits, yLimits)
    style = figure_style();
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'Color', 'w', 'LineWidth', style.axesLineWidth, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Layer', 'top', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'XLim', xLimits, 'YLim', yLimits, 'XTick', [xLimits(1) mean(xLimits) xLimits(2)], 'YTick', [yLimits(1) mean(yLimits) yLimits(2)]);
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
end

function local_style_map_axis(ax, absoluteTimeLimits)
    style = figure_style();
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'Color', 'w', 'LineWidth', style.axesLineWidth, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Layer', 'top', 'FontName', style.fontName, 'FontSize', style.tickFontSize, 'XLim', [0 350], 'YLim', absoluteTimeLimits, 'XTick', [0 175 350], 'YTick', [absoluteTimeLimits(1) mean(absoluteTimeLimits) absoluteTimeLimits(2)]);
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
    if strcmp(baseName, 'figure_ligand_profile')
        exportgraphics(fig, paths.png, 'Resolution', resolution, 'BackgroundColor', 'white', 'Padding', 'figure');
        exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Padding', 'figure');
    else
        exportgraphics(fig, paths.png, 'Resolution', resolution, 'BackgroundColor', 'white');
        exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
    savefig(fig, paths.fig);
end

function local_require_treatment(treatmentResult, kind)
    if isempty(treatmentResult)
        error('plot_model_figure:MissingTreatment', 'Figure kind "%s" requires a treatment result.', kind);
    end
end
