function sweeps = compute_supplementary_sweeps(experiment, project, varargin)
    %COMPUTE_SUPPLEMENTARY_SWEEPS Compute or load the shared S4A-C cache.
    p = inputParser;
    addParameter(p, 'Recompute', false, @islogical);
    parse(p, varargin{:});
    cfg = sweep_settings();
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    interval = analysis.densityIntervalsHours(end, :) - analysis.experimentStartHour;
    base = model_to_solver_parameters(model_parameters('control'), numerical_settings());
    sweeps = run_parameter_sweeps( 'OutputDir', fullfile(project.supplementaryRoot, 'parameter_sweeps'), 'BaseParameters', base, 'BaselineId', 'publication_model', 'ControlReplicates', experiment.control, 'FeatureFile', experiment.paths.control, 'ModelTimeWindowHours', interval, 'EdgeRadius', edgeRadius, 'DistanceFromEdgeLimit', [0 edgeRadius], 'LambdaValues', cfg.lambda, 'SourceRadiusValues', cfg.sourceRadius, 'BetaValues', cfg.beta, 'MesodermDiffusivityValues', cfg.mesodermDiffusivity, 'AlphaValues', cfg.alpha, 'GammaValues', cfg.gamma, 'Recompute', p.Results.Recompute, 'ShowProgress', true);
end
