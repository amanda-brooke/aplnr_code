function summary = compute_segment_velocity_maps(xPaths, yPaths, tPaths, geometry, varargin)

    p = inputParser;
    addParameter(p, 'TimeOffsetHours', 0, @isscalar);
    addParameter(p, 'DistanceEdges', 0:25:350, @isnumeric);
    addParameter(p, 'TimeEdgesHours', 24:4:48, @isnumeric);
    addParameter(p, 'MinimumCount', 10, @isscalar);
    parse(p, varargin{:});
    dEdges = double(p.Results.DistanceEdges(:));
    tEdges = double(p.Results.TimeEdgesHours(:));

    distance = []; time = []; speed = []; inward = []; angle = [];
    for i = 1:numel(xPaths)
        x = double(xPaths{i}(:)); y = double(yPaths{i}(:));
        tSeconds = double(tPaths{i}(:));
        dtHours = diff(tSeconds) ./ 3600;
        xMid = 0.5 .* (x(1:end-1) + x(2:end));
        yMid = 0.5 .* (y(1:end-1) + y(2:end));
        tMidSeconds = 0.5 .* (tSeconds(1:end-1) + tSeconds(2:end));
        [rMid, dMid, cx, cy] = radial_coordinates(xMid, yMid, geometry);
        dx = diff(x); dy = diff(y);
        outwardX = (xMid - cx) ./ max(rMid, eps);
        outwardY = (yMid - cy) ./ max(rMid, eps);
        edgeInward = -(dx .* outwardX + dy .* outwardY) ./ dtHours;
        edgeSpeed = abs(edgeInward);
        inwardX = -outwardX; inwardY = -outwardY;
        edgeAngle = atan2(dx .* inwardY - dy .* inwardX, dx .* inwardX + dy .* inwardY);
        tMidHours = tMidSeconds ./ 3600 + double(p.Results.TimeOffsetHours);
        ok = isfinite(dMid) & isfinite(tMidHours) & isfinite(edgeSpeed) & isfinite(edgeInward) & dtHours > 0;
        distance = [distance; dMid(ok)]; 
        time = [time; tMidHours(ok)]; 
        speed = [speed; edgeSpeed(ok)]; 
        inward = [inward; edgeInward(ok)]; 
        angle = [angle; edgeAngle(ok)]; 
    end

    dBin = discretize(distance, dEdges);
    tBin = discretize(time, tEdges);
    valid = isfinite(dBin) & isfinite(tBin);
    shape = [numel(dEdges) - 1, numel(tEdges) - 1];
    counts = accumarray([dBin(valid), tBin(valid)], 1, shape, @sum, 0);
    meanSpeed = local_mean_map(dBin(valid), tBin(valid), speed(valid), shape);
    meanInward = local_mean_map(dBin(valid), tBin(valid), inward(valid), shape);
    minimumCount = double(p.Results.MinimumCount);
    meanSpeed(counts < minimumCount) = NaN;
    meanInward(counts < minimumCount) = NaN;
    summary = struct('distanceFromEdgeEdges', dEdges, 'distanceFromEdgeCenters', 0.5 .* (dEdges(1:end-1) + dEdges(2:end)), 'rCenters', double(geometry.edgeRadius) - 0.5 .* (dEdges(1:end-1) + dEdges(2:end)), 'tEdges', tEdges, 'tCenters', 0.5 .* (tEdges(1:end-1) + tEdges(2:end)), 'avgV', meanSpeed, 'avgSpeed', meanSpeed, 'avgInwardVelocity', meanInward, 'counts', counts, 'angles', angle, 'minimumCount', minimumCount, 'meta', struct('version', 4, 'quantity', 'radial_speed', 'projection', 'midpoint_radial_unit_vector'));
end

function values = local_mean_map(dBin, tBin, data, shape)
    sumValues = accumarray([dBin, tBin], data, shape, @sum, 0);
    counts = accumarray([dBin, tBin], 1, shape, @sum, 0);
    values = sumValues ./ max(counts, 1);
    values(counts == 0) = NaN;
end
