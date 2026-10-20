function paths = plot_density_comparison( controlResult, treatmentResult, controlReplicates, treatmentReplicates, varargin)
    %PLOT_DENSITY_COMPARISON Plot model and data by time interval.

    style = figure_style();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    p = inputParser;
    addParameter(p, 'OutputDir', fullfile(pwd, 'publication_model_figures'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'Visible', style.showFigures, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Resolution', style.exportResolution, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'IntervalsElapsedHours', analysis.densityIntervalsHours - analysis.experimentStartHour, @(x) isnumeric(x) && size(x, 2) == 2 && size(x, 1) == 3);
    addParameter(p, 'DensityLimits', style.densityLimits, @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
    addParameter(p, 'EdgeRadius', edgeRadius, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'RadialBinCount', analysis.radialBinCount, @(x) isnumeric(x) && isscalar(x) && x >= 2 && mod(x, 1) == 0);
    parse(p, varargin{:});
    opts = p.Results;

    conditions = { struct('name', 'Control', 'file', 'control', 'color', style.controlColor, 'model', controlResult, 'replicates', controlReplicates), struct('name', 'shAPLNR', 'file', 'shaplnr', 'color', style.treatmentColor, 'model', treatmentResult, 'replicates', treatmentReplicates)};
    outputDir = char(string(opts.OutputDir));
    if ~isfolder(outputDir), mkdir(outputDir); end
    visibility = 'off';
    if logical(opts.Visible), visibility = 'on'; end
    for conditionIndex = 1:numel(conditions)
        condition = conditions{conditionIndex};
        fig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 style.densityOverlayFigureSize], 'MenuBar', 'none', 'ToolBar', 'none');
        layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        for intervalIndex = 1:3
            ax = nexttile(layout);
            interval = double(opts.IntervalsElapsedHours(intervalIndex, :));
            [xExperiment, experimentMean, experimentLow, experimentHigh] = local_experiment_summary(condition.replicates, interval, double(opts.EdgeRadius), double(opts.RadialBinCount));
            [xModel, modelMean, modelLow, modelHigh] = local_model_summary(condition.model, interval);
            hold(ax, 'on');
            local_band(ax, xExperiment, experimentLow, experimentHigh, condition.color, style.densityExperimentBandAlpha);
            local_band(ax, xModel, modelLow, modelHigh, condition.color, style.densityModelBandAlpha);
            plot(ax, xExperiment, experimentMean, '--', 'Color', condition.color, 'LineWidth', style.densityExperimentLineWidth, 'DisplayName', 'Experiment');
            plot(ax, xModel, modelMean, '-', 'Color', condition.color, 'LineWidth', style.profileLineWidth, 'DisplayName', 'Model');
            hold(ax, 'off');
            xlim(ax, [0 double(opts.EdgeRadius)]);
            pbaspect(ax, [1 1 1]);
            ylim(ax, double(opts.DensityLimits));
            xticks(ax, [0 double(opts.EdgeRadius) / 2 double(opts.EdgeRadius)]);
            yticks(ax, [0 mean(opts.DensityLimits) opts.DensityLimits(2)]);
            if intervalIndex == 1
                ax.YAxis.Exponent = -5;
            else
                ax.YAxis.Exponent = 0;
            end
            title(ax, sprintf('%g-%g h', interval(1) + 24, interval(2) + 24), 'FontWeight', 'normal', 'FontSize', style.densityOverlayLabelFontSize, 'Color', 'k');
            if conditionIndex == 2, title(ax, ''); end
            if intervalIndex ~= 1
                ax.YTick = [];
            end
            set(ax, 'FontName', style.fontName, 'FontSize', style.densityOverlayTickFontSize, 'LineWidth', style.axesLineWidth, 'TickDir', 'in', 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        end
        if conditionIndex == 2
            xlabel(layout, 'Distance from edge (\mum)', 'FontName', style.fontName, 'FontSize', style.densityOverlayLabelFontSize, 'Color', 'k');
        end
        ylabel(layout, 'Areal density (1/\mum^2)', 'FontName', style.fontName, 'FontSize', style.densityOverlayLabelFontSize, 'Color', 'k');
        title(layout, condition.name, 'FontName', style.fontName, 'FontSize', style.densityOverlayLabelFontSize, 'FontWeight', 'normal', 'Color', 'k');
        legendHandle = legend(layout.Children(end), {'Experiment', 'Model'}, 'Location', 'northeast', 'Box', 'off');
        set(legendHandle, 'TextColor', 'k', 'Color', 'w', 'FontSize', style.densityOverlayLegendFontSize);
        baseName = ['figure_model_experiment_density_' condition.file];
        paths.(condition.file) = local_export( fig, outputDir, baseName, double(opts.Resolution));
        if ~logical(opts.Visible), close(fig); end
    end
end

function [x, profile, lowerProfile, upperProfile] = local_experiment_summary(replicates, interval, edgeRadius, nBins)
    x = linspace(0, edgeRadius, nBins).';
    profiles = nan(nBins, numel(replicates));
    for i = 1:numel(replicates)
        [radius, time, density] = experimental_areal_density_profile(replicates(i).density);
        elapsed = time - min(time);
        keep = elapsed >= interval(1) & elapsed <= interval(2);
        if any(keep)
            replicateProfile = mean(density(keep, :), 1, 'omitnan').';
        else
            replicateProfile = interp1(elapsed, density, mean(interval), 'linear', 'extrap').';
        end
        distance = edgeRadius - radius;
        [distance, order] = sort(distance(:));
        replicateProfile = replicateProfile(order);
        valid = isfinite(distance) & isfinite(replicateProfile);
        profiles(:, i) = interp1(distance(valid), replicateProfile(valid), x, 'linear', NaN);
    end
    profile = mean(profiles, 2, 'omitnan');
    lowerProfile = min(profiles, [], 2, 'omitnan');
    upperProfile = max(profiles, [], 2, 'omitnan');
    % Anchor the displayed experimental profile and band at the colony edge.
    profile(1) = 0;
    lowerProfile(1) = 0;
    upperProfile(1) = 0;
end

function [x, profile, lowerProfile, upperProfile] = local_model_summary(model, interval)
    [~, time, density] = model_areal_density_profile(model);
    x = double(model.params.R) - double(model.rs(:));
    [x, order] = sort(x);
    keep = time >= interval(1) & time <= interval(2);
    if any(keep)
        values = density(keep, order);
    else
        values = interp1(time, density(:, order), mean(interval), 'linear', 'extrap');
    end
    profile = mean(values, 1, 'omitnan').';
    lowerProfile = min(values, [], 1, 'omitnan').';
    upperProfile = max(values, [], 1, 'omitnan').';
end

function local_band(ax, x, lowerProfile, upperProfile, color, alphaValue)
    valid = isfinite(x) & isfinite(lowerProfile) & isfinite(upperProfile);
    if nnz(valid) < 2, return; end
    fill(ax, [x(valid); flipud(x(valid))], [lowerProfile(valid); flipud(upperProfile(valid))], color, 'FaceAlpha', alphaValue, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function paths = local_export(fig, outputDir, baseName, resolution)
    paths.png = fullfile(outputDir, [baseName '.png']);
    paths.pdf = fullfile(outputDir, [baseName '.pdf']);
    paths.fig = fullfile(outputDir, [baseName '.fig']);
    exportgraphics(fig, paths.png, 'Resolution', resolution, 'BackgroundColor', 'white');
    exportgraphics(fig, paths.pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, paths.fig);
end
