function style = figure_style()

    style.densityLimits = [0 1.6e-5];
    style.speedLimits = [0 30];
    style.velocityLimits = [-5 10];
    style.supplementVelocityLimits = [-5 15];
    style.supplementVelocityXLimits = [0 200];
    style.apparentSpeedFrameIntervalHours = 1 / 6;
    style.showFigures = true;
    style.exportResolution = 300;
    style.conditionLabels = {'Control', 'shAPLNR'};
    style.intervalVelocityFigureSize = [5.2 3.7];
    style.intervalColors = [0.42 0.18 0.68; 0.92 0.22 0.22; 0.95 0.56 0.20];
    style.initialConditionFigureSize = [4.0 3.4];
    style.initialConditionAspectRatio = [0.5 1 1];
    style.initialConditionLineWidth = 2.0;
    style.initialConditionColor = [0.68 0.25 0.18];
    style.sweepFigureSize = [7.2 3.4];

    % Shared typography and line weights.
    style.fontName = 'Arial';
    style.tickFontSize = 10;
    style.labelFontSize = 12;
    style.legendFontSize = 9;
    style.axesLineWidth = 0.9;
    style.colorbarLineWidth = 0.8;
    style.colorbarLabelFontSize = 11;
    style.profileLineWidth = 2.7;
    style.referenceLineWidth = 1.4;
    style.controlColor = [0.68 0.25 0.18];
    style.treatmentColor = [0 0 0];
    style.ligandFigureSize = [3.5 3.5];
    style.speedMapFigureSize = [6.6 3.0];

    % S4A-C: corner colors are lower-left, lower-right, upper-left, upper-right.
    style.sweepCornerColors = [0.000 0.447 0.698; 0.902 0.624 0.000; 0.800 0.475 0.655; 0.000 0.620 0.451];
    style.sweepDensityLineWidth = 1.8;
    style.sweepLigandLineWidth = 1.6;
    style.sweepBaselineDensityLineWidth = 2.4;
    style.sweepBaselineLigandLineWidth = 2.2;
    style.sweepBoxLineWidth = 1.8;

    % Figure 4D and S4E: pooled and interval velocity comparisons.
    style.velocityFigureSize = [6.8 3.6];
    style.velocityBandAlpha = 0.12;
    style.velocityReplicateFade = 0.60;
    style.velocityReplicateLineWidth = 1;
    style.intervalVelocityLineWidth = 2.4;
    style.intervalVelocityModelLineStyle = '-';
    style.intervalVelocityModelColor = [0 0 0];
    style.intervalVelocityBandAlpha = 0.10;
    style.supplementTickFontSize = 13;
    style.supplementAxesLineWidth = 1.4;

    % S4D: boundary-condition comparison.
    style.boundaryFigureSize = [4.6 7.6];
    style.boundaryDensityLineWidth = 3.0;
    style.boundaryLigandLineWidth = 2.6;
    style.boundaryLigandColor = [0.35 0.35 0.35];
    style.boundaryLigandLegendColor = [0.45 0.45 0.45];
    style.boundaryTickFontSize = 15;
    style.boundaryLabelFontSize = 18;
    style.boundaryLegendFontSize = 14;

    % S4F-G: initial conditions and state decomposition.
    style.initialHistogramColor = [0.72 0.72 0.72];
    style.stateFigureSize = [4.0 3.4];
    style.stateColors = [0 0 0; 0.82 0.55 0.50; 0.56 0.15 0.10];
    style.stateLineStyles = {':', '--', '-'};
    style.stateLineWidth = 2.5;

    % Figure 4E: model-experiment density overlays.
    style.densityOverlayFigureSize = [9.2 3.9];
    style.densityOverlayTickFontSize = 16;
    style.densityOverlayLabelFontSize = 18;
    style.densityOverlayLegendFontSize = 14;
    style.densityExperimentLineWidth = 2.4;
    style.densityExperimentBandAlpha = 0.12;
    style.densityModelBandAlpha = 0.05;
end
