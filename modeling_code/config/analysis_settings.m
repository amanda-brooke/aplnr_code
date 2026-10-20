function analysis = analysis_settings()
    analysis.experimentStartHour = 24;
    analysis.timeWindowHours = [24 54];
    analysis.densityIntervalsHours = [24 34; 34 44; 44 54];
    analysis.initialDensityWindowHours = [24 28];
    analysis.radialBinCount = 60;
    analysis.densityTimePointCount = 60;
    analysis.minimumCellsPerFrameBin = 2;
    analysis.minimumFramesPerTrackSegment = 10;
    analysis.velocity.minimumDisplacementsPerFrameBin = 2;
    analysis.velocity.minimumValidFramesPerReplicateBin = 10;
    analysis.velocity.minimumReplicates = 4;
end
