function options = publication_velocity_options(project)
    %PUBLICATION_VELOCITY_OPTIONS Shared analysis conventions for Figure 4D and S4E.
    analysis = analysis_settings();
    edgeRadius = model_parameters('control').domainRadius;
    options = { 'RadialBinWidth', edgeRadius / analysis.radialBinCount, 'EdgeRadius', edgeRadius, 'MinimumFramesPerSegment', analysis.minimumFramesPerTrackSegment, 'MinimumDisplacementsPerFrameBin', analysis.velocity.minimumDisplacementsPerFrameBin, 'MinimumValidFramesPerBin', analysis.velocity.minimumValidFramesPerReplicateBin, 'MinimumReplicates', analysis.velocity.minimumReplicates, 'DataRoot', project.rawDataRoot};
end
