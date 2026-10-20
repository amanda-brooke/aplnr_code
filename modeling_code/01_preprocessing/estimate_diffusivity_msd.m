function out = estimate_diffusivity_msd(varargin)
    %ESTIMATE_DIFFUSIVITY_MSD Estimate D_m from treatment segments.
    %
    % Uses the two-dimensional relation MSD(tau) = 4 D tau + c. MSD values are
    % first averaged within each track segment and then across segments.
    % Tracking graphs are split at divisions so shared ancestors are not repeated.

    project = model_pipeline_setup();
    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'DataRoot', project.rawDataRoot, @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'OutputDir', fullfile(project.mainFigureRoot, 'treatment_diffusivity_msd'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'FileStems', {'dox_1', 'dox_2', 'dox_3', 'dox_4'}, @iscell);
    addParameter(p, 'LagStepHours', 1 / 6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MaxLagHours', 6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'FitLagHours', [0.5 3], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'BootstrapSamples', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'Visible', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    parse(p, varargin{:});

    dataRoot = char(string(p.Results.DataRoot));
    outputDir = char(string(p.Results.OutputDir));
    stems = p.Results.FileStems;
    lagStep = double(p.Results.LagStepHours);
    maxLag = double(p.Results.MaxLagHours);
    fitRange = sort(double(p.Results.FitLagHours(:).'));
    nBoot = round(double(p.Results.BootstrapSamples));
    if ~isfolder(outputDir), mkdir(outputDir); end

    edges = (-0.5 * lagStep):lagStep:(maxLag + 0.5 * lagStep);
    tau = (0:lagStep:maxLag).';
    nRep = numel(stems);
    replicates = struct([]);
    rng(double(p.Results.RandomSeed));

    for i = 1:nRep
        files = resolve_dataset_files('offset_24', 'treatment', 'DataRoot', dataRoot, 'FileStem', stems{i});
        spots = read_data_table(files.spotFile);
        edgesTable = read_data_table(files.edgeFile);
        [xPaths, yPaths, tPaths, segmentMetadata] = build_track_segments(edgesTable, spots, 'MinimumFrames', 10);
        trackMsd = local_track_msd(xPaths, yPaths, tPaths, edges);
        meanMsd = mean(trackMsd, 1, 'omitnan').';
        semMsd = std(trackMsd, 0, 1, 'omitnan').' ./ sqrt(sum(isfinite(trackMsd), 1).');
        fitMask = tau >= fitRange(1) & tau <= fitRange(2) & isfinite(meanMsd);
        [D, intercept, fitMsd] = local_fit(tau, meanMsd, fitMask);
        bootD = nan(nBoot, 1);
        for b = 1:nBoot
            sample = randi(size(trackMsd, 1), size(trackMsd, 1), 1);
            bootMean = mean(trackMsd(sample, :), 1, 'omitnan').';
            bootD(b) = local_fit(tau, bootMean, fitMask);
        end
        ci = prctile(bootD, [2.5 97.5]);
        replicate = struct('stem', stems{i}, 'segmentMetadata', segmentMetadata, 'nSegments', size(trackMsd, 1), 'tauHours', tau, 'trackMsd', trackMsd, 'meanMsd', meanMsd, 'semMsd', semMsd, 'fitMask', fitMask, 'fitMsd', fitMsd, 'D', D, 'intercept', intercept, 'bootstrapD', bootD, 'ci95', ci);
        if i == 1
            replicates = replicate;
        else
            replicates(i, 1) = replicate;
        end
    end

    D = vertcat(replicates.D);
    ci = vertcat(replicates.ci95);
    summary = table(string(stems(:)), vertcat(replicates.nSegments), D, ci(:, 1), ci(:, 2), 'VariableNames', {'replicate', 'nSegments', 'D_um2_per_h', 'ci95Low', 'ci95High'});
    writetable(summary, fullfile(outputDir, 'treatment_diffusivity_msd_summary.csv'));
    save(fullfile(outputDir, 'treatment_diffusivity_msd.mat'), 'replicates', 'summary', 'fitRange', 'lagStep', 'maxLag', '-v7.3');

    figures = local_plot(replicates, fitRange, outputDir, logical(p.Results.Visible));
    out = struct('replicates', replicates, 'summary', summary, 'fitLagHours', fitRange, 'figures', figures, 'outputDir', outputDir);
end

function trackMsd = local_track_msd(xPaths, yPaths, tPathsSeconds, edges)
    nTracks = numel(xPaths);
    nBins = numel(edges) - 1;
    trackMsd = nan(nTracks, nBins);
    for i = 1:nTracks
        x = double(xPaths{i}(:)); y = double(yPaths{i}(:));
        t = double(tPathsSeconds{i}(:)) / 3600;
        sums = zeros(1, nBins); counts = zeros(1, nBins);
        for j = 1:(numel(t) - 1)
            dt = t((j + 1):end) - t(j);
            keep = dt <= edges(end);
            if ~any(keep), continue; end
            dr2 = (x((j + 1):end) - x(j)).^2 + (y((j + 1):end) - y(j)).^2;
            bins = discretize(dt(keep), edges);
            values = dr2(keep);
            valid = isfinite(bins) & isfinite(values);
            sums = sums + accumarray(bins(valid), values(valid), [nBins 1], @sum, 0).';
            counts = counts + accumarray(bins(valid), 1, [nBins 1], @sum, 0).';
        end
        trackMsd(i, counts > 0) = sums(counts > 0) ./ counts(counts > 0);
    end
end

function [D, intercept, fitMsd] = local_fit(tau, msd, mask)
    D = nan; intercept = nan; fitMsd = nan(size(tau));
    if sum(mask) < 2, return; end
    p = polyfit(tau(mask), msd(mask), 1);
    D = p(1) / 4;
    intercept = p(2);
    fitMsd(mask) = polyval(p, tau(mask));
end

function files = local_plot(reps, fitRange, outputDir, visible)
    fig = figure('Color', 'w', 'Visible', local_on_off(visible), 'Units', 'inches', 'Position', [1 1 3.8 3.8]);
    ax = axes(fig, 'Position', [0.21 0.19 0.74 0.76]); hold(ax, 'on');
    tau = reps(1).tauHours;
    replicateMsd = horzcat(reps.meanMsd);
    meanMsd = mean(replicateMsd, 2, 'omitnan');
    sdMsd = std(replicateMsd, 0, 2, 'omitnan');
    validBand = isfinite(meanMsd) & isfinite(sdMsd);
    fill(ax, [tau(validBand); flipud(tau(validBand))], [meanMsd(validBand) - sdMsd(validBand); flipud(meanMsd(validBand) + sdMsd(validBand))], [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    for i = 1:numel(reps)
        plot(ax, tau, reps(i).meanMsd, '-', 'Color', [0.68 0.68 0.68], 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    plot(ax, tau, meanMsd, 'k-', 'LineWidth', 2.8);
    fitMask = tau >= fitRange(1) & tau <= fitRange(2) & isfinite(meanMsd);
    [~, ~, fitMsd] = local_fit(tau, meanMsd, fitMask);
    mesodermRed = [0.78 0.28 0.22];
    plot(ax, tau, fitMsd, '--', 'Color', mesodermRed, 'LineWidth', 3.6);

    xlim(ax, [0 tau(end)]); ylim(ax, [0 inf]);
    xlabel(ax, 'Lag time (h)'); ylabel(ax, 'MSD (\mum^2)');
    local_style(ax);

    png = fullfile(outputDir, 'treatment_diffusivity_msd.png');
    pdf = fullfile(outputDir, 'treatment_diffusivity_msd.pdf');
    figPath = fullfile(outputDir, 'treatment_diffusivity_msd.fig');
    exportgraphics(fig, png, 'Resolution', 300);
    exportgraphics(fig, pdf, 'ContentType', 'vector');
    savefig(fig, figPath);
    files = struct('png', png, 'pdf', pdf, 'fig', figPath);
end

function local_style(ax)
    set(ax, 'FontName', 'Arial', 'FontSize', 9, 'LineWidth', 1, 'Box', 'on', 'TickDir', 'in', 'XColor', 'k', 'YColor', 'k', 'Color', 'w', 'Layer', 'top');
    grid(ax, 'off');
    ax.FontSize = 13; ax.LineWidth = 1.2;
    ax.XLabel.FontSize = 16; ax.YLabel.FontSize = 16;
    ax.XLabel.Color = 'k'; ax.YLabel.Color = 'k';
    pbaspect(ax, [1 1 1]);
end

function value = local_on_off(tf)
    if tf, value = 'on'; else, value = 'off'; end
end
