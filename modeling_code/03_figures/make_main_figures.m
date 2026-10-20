function outputs = make_main_figures(results, experiment, project, varargin)
    %MAKE_MAIN_FIGURES Generate all publication-model main-text components.

    p = inputParser;
    addParameter(p, 'Panels', {'4B', '4C', '4D', '4E', '4F'}, @iscellstr);
    parse(p, varargin{:});
    panels = p.Results.Panels;
    assert(all(ismember(panels, {'4B', '4C', '4D', '4E', '4F'})), 'Unknown Figure 4 panel.');
    outputs = struct();
    analysis = analysis_settings();
    style = figure_style();
    edgeRadius = model_parameters('control').domainRadius;
    intervalsElapsed = analysis.densityIntervalsHours - analysis.experimentStartHour;

    if ismember('4B', panels)
        outputs.ligand = plot_model_figure( 'ligand_drive', results.control, [], 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures);
    end
    if ismember('4E', panels)
        outputs.densityOverlays = plot_density_comparison( results.control, results.shAPLNR, experiment.control, experiment.shAPLNR, 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures, 'IntervalsElapsedHours', intervalsElapsed, 'DensityLimits', style.densityLimits, 'EdgeRadius', edgeRadius, 'RadialBinCount', analysis.radialBinCount);
    end
    if ismember('4F', panels)
        outputs.modelMaps = plot_model_figure( 'radial_speed_maps', results.control, results.shAPLNR, 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures, 'DynamicsAbsoluteTimeLimits', analysis.timeWindowHours, 'SpeedLimits', style.speedLimits, 'ApparentSpeedFrameIntervalHours', style.apparentSpeedFrameIntervalHours);
        outputs.experimentalMaps = plot_experimental_speed_maps(experiment.control, experiment.shAPLNR, 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures, 'AbsoluteTimeLimits', analysis.timeWindowHours, 'SpeedLimits', style.speedLimits, 'EdgeRadius', edgeRadius, 'RadialBinCount', analysis.radialBinCount);
    end
    if ismember('4D', panels)
        velocityOptions = publication_velocity_options(project);
        outputs.inwardVelocity = plot_velocity_comparison( results.control, results.shAPLNR, experiment.control, experiment.shAPLNR, 'TimeWindowHours', analysis.timeWindowHours, 'VelocityLimits', style.velocityLimits, velocityOptions{:}, 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures);
    end
    if ismember('4C', panels)
        outputs.diffusivity = estimate_diffusivity_msd( 'DataRoot', project.rawDataRoot, 'OutputDir', project.mainFigureRoot, 'Visible', style.showFigures);
    end
end
