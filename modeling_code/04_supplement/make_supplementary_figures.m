function outputs = make_supplementary_figures( results, experiment, project, varargin)
    %MAKE_SUPPLEMENTARY_FIGURES Plot supplied sweeps, boundary tests, and panels E-G.

    p = inputParser;
    addRequired(p, 'results', @isstruct);
    addRequired(p, 'experiment', @isstruct);
    addRequired(p, 'project', @isstruct);
    addParameter(p, 'Panels', {'S4A', 'S4B', 'S4C', 'S4D', 'S4E', 'S4F', 'S4G'}, @iscellstr);
    addParameter(p, 'Sweeps', struct(), @isstruct);
    parse(p, results, experiment, project, varargin{:});
    panels = p.Results.Panels;
    assert(all(ismember(panels, {'S4A', 'S4B', 'S4C', 'S4D', 'S4E', 'S4F', 'S4G'})), 'Unknown Supplement S4 panel.');
    outputs = struct();
    analysis = analysis_settings();
    style = figure_style();
    edgeRadius = model_parameters('control').domainRadius;
    intervalsElapsed = analysis.densityIntervalsHours - analysis.experimentStartHour;
    sweepDir = fullfile(project.supplementaryRoot, 'parameter_sweeps');
    boundaryDir = fullfile(project.supplementaryRoot, 'boundary_sensitivity');
    panelDir = fullfile(project.supplementaryRoot, 'figure_panels');
    for path = {sweepDir, boundaryDir, panelDir}
        if ~isfolder(path{1}), mkdir(path{1}); end
    end

    base = model_to_solver_parameters(model_parameters('control'), numerical_settings());
    if any(ismember(panels, {'S4A', 'S4B', 'S4C'}))
        sweeps = p.Results.Sweeps;
        if isempty(fieldnames(sweeps))
            error('make_supplementary_figures:MissingSweeps', 'Compute or load sweeps with compute_supplementary_sweeps and pass them using Sweeps for S4A-C.');
        end
        names = {'sourceLambda', 'betaDiffusivity', 'alphaGamma'};
        selected = names(ismember({'S4A', 'S4B', 'S4C'}, panels));
        outputs.sweeps = plot_parameter_sweeps( sweeps, 'OutputDir', sweepDir, 'Visible', style.showFigures, 'DensityLimits', style.densityLimits, 'Panels', selected);
    end

    if ismember('S4D', panels)
        boundary = struct('baselineId', 'publication_model', 'base', base, 'outputDir', boundaryDir);
        outputs.boundary = plot_boundary_sensitivity( boundary, 'OutputDir', boundaryDir, 'Visible', style.showFigures, 'DensityIntervalsElapsedHours', intervalsElapsed, 'DensityLimits', style.densityLimits);
    end

    if ismember('S4E', panels)
        velocityOptions = publication_velocity_options(project);
        outputs.velocityIntervals = plot_velocity_comparison( results.control, results.shAPLNR, experiment.control, experiment.shAPLNR, 'IntervalTimeWindowsHours', analysis.densityIntervalsHours, 'VelocityLimits', style.supplementVelocityLimits, velocityOptions{:}, 'OutputDir', panelDir, 'Visible', style.showFigures);
    end
    if ismember('S4F', panels)
        outputs.initialConditions = plot_initial_and_state_profiles( results.control, results.shAPLNR, experiment.control, experiment.shAPLNR, 'Panels', {'initial_conditions'}, 'OutputDir', panelDir, 'Visible', style.showFigures, 'DensityIntervalsElapsedHours', intervalsElapsed, 'EdgeRadius', edgeRadius);
    end
    if ismember('S4G', panels)
        outputs.stateDecomposition = plot_initial_and_state_profiles( results.control, results.shAPLNR, experiment.control, experiment.shAPLNR, 'Panels', {'state_decomposition'}, 'OutputDir', panelDir, 'Visible', style.showFigures, 'DensityIntervalsElapsedHours', intervalsElapsed, 'EdgeRadius', edgeRadius);
    end
    outputs.directories = struct('sweeps', sweepDir, 'boundarySensitivity', boundaryDir, 'figurePanels', panelDir);
end
