function out = compute_density_timecourse(radius, timeHours, rEdges, tGridHours, varargin)
    %COMPUTE_DENSITY_TIMECOURSE Areal density from unique spot detections.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'MinimumBinCount', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && mod(x, 1) == 0);
    parse(p, varargin{:});

    radius = double(radius(:));
    timeHours = double(timeHours(:));
    rEdges = double(rEdges(:));
    tGridHours = double(tGridHours(:));
    minimumBinCount = double(p.Results.MinimumBinCount);
    assert(numel(radius) == numel(timeHours), 'compute_density_timecourse:SizeMismatch', 'radius and timeHours must have the same number of elements.');

    ok = isfinite(radius) & isfinite(timeHours);
    radius = radius(ok);
    timeHours = timeHours(ok);
    [frameTimes, ~, frameIndex] = unique(timeHours);
    nFrames = numel(frameTimes);
    nBins = numel(rEdges) - 1;
    frameDensity = nan(nFrames, nBins);
    rawFrameCounts = zeros(nFrames, 1);
    retainedFrameCounts = zeros(nFrames, 1);
    rawBinCounts = zeros(nFrames, nBins);
    retainedBinCounts = zeros(nFrames, nBins);
    annulusAreas = pi .* (rEdges(2:end).^2 - rEdges(1:end-1).^2);

    for i = 1:nFrames
        values = radius(frameIndex == i);
        counts = histcounts(values, rEdges);
        rawBinCounts(i, :) = counts;
        rawFrameCounts(i) = numel(values);
        counts(counts < minimumBinCount) = 0;
        retainedBinCounts(i, :) = counts;
        retainedFrameCounts(i) = sum(counts);
        denominator = max(retainedFrameCounts(i), 1);
        frameDensity(i, :) = counts ./ (denominator .* annulusAreas.');
    end

    if nFrames == 1
        densityOver = repmat(frameDensity, numel(tGridHours), 1);
        counts = repmat(retainedFrameCounts, numel(tGridHours), 1);
        rawCounts = repmat(rawFrameCounts, numel(tGridHours), 1);
        rawCountsByBin = repmat(rawBinCounts, numel(tGridHours), 1);
        retainedCountsByBin = repmat(retainedBinCounts, numel(tGridHours), 1);
    else
        densityOver = interp1(frameTimes, frameDensity, tGridHours, 'linear', NaN);
        counts = interp1(frameTimes, retainedFrameCounts, tGridHours, 'linear', NaN);
        rawCounts = interp1(frameTimes, rawFrameCounts, tGridHours, 'linear', NaN);
        rawCountsByBin = interp1(frameTimes, rawBinCounts, tGridHours, 'nearest', NaN);
        retainedCountsByBin = interp1( frameTimes, retainedBinCounts, tGridHours, 'nearest', NaN);
    end

    out = struct( 'rEdges', rEdges, 'rCenters', 0.5 .* (rEdges(1:end-1) + rEdges(2:end)), 'tHours', tGridHours, 'densityOver', densityOver, 'counts', counts, 'rawCounts', rawCounts, 'rawCountsByBin', rawCountsByBin, 'retainedCountsByBin', retainedCountsByBin, 'densityMode', 'area', 'source', 'unique_spots', 'minimumBinCount', minimumBinCount, 'nUniqueSpots', numel(radius));
end
