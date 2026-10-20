function [xPaths, yPaths, tPaths, metadata] = build_track_segments(edgeTable, spotTable, varargin)

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'MinimumFrames', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && mod(x, 1) == 0);
    parse(p, varargin{:});
    minimumFrames = double(p.Results.MinimumFrames);

    requiredEdges = {'SPOT_SOURCE_ID', 'SPOT_TARGET_ID'};
    requiredSpots = {'ID', 'POSITION_X', 'POSITION_Y', 'POSITION_T'};
    assert(all(ismember(requiredEdges, edgeTable.Properties.VariableNames)), 'build_track_segments:MissingEdgeColumns', 'Edge table is missing required source/target columns.');
    assert(all(ismember(requiredSpots, spotTable.Properties.VariableNames)), 'build_track_segments:MissingSpotColumns', 'Spot table is missing required ID/position/time columns.');

    source = double(edgeTable.SPOT_SOURCE_ID);
    target = double(edgeTable.SPOT_TARGET_ID);
    validEdge = isfinite(source) & isfinite(target);
    source = source(validEdge);
    target = target(validEdge);
    allNodes = unique([source; target]);
    [~, sourceNode] = ismember(source, allNodes);
    [~, targetNode] = ismember(target, allNodes);
    inDegree = accumarray(targetNode, 1, [numel(allNodes), 1]);
    outDegree = accumarray(sourceNode, 1, [numel(allNodes), 1]);
    outgoingEdges = accumarray(sourceNode, (1:numel(source)).', [numel(allNodes), 1], @(v) {v}, {[]});

    spotId = double(spotTable.ID);
    validSpot = isfinite(spotId);
    idToRow = containers.Map(spotId(validSpot), find(validSpot));
    starts = find(inDegree(sourceNode) ~= 1 | outDegree(sourceNode) ~= 1);
    xPaths = cell(numel(source), 1);
    yPaths = cell(numel(source), 1);
    tPaths = cell(numel(source), 1);
    segmentNodes = cell(numel(source), 1);
    used = false(numel(source), 1);
    nSegments = 0;

    for startIndex = starts(:).'
        if used(startIndex), continue; end
        nodes = [source(startIndex), target(startIndex)];
        used(startIndex) = true;
        currentNode = targetNode(startIndex);
        while inDegree(currentNode) == 1 && outDegree(currentNode) == 1
            candidates = outgoingEdges{currentNode};
            edgeIndex = candidates(find(~used(candidates), 1, 'first'));
            if isempty(edgeIndex), break; end
            nodes(end + 1) = target(edgeIndex); 
            used(edgeIndex) = true;
            currentNode = targetNode(edgeIndex);
        end
        nSegments = nSegments + 1;
        [xPaths{nSegments}, yPaths{nSegments}, tPaths{nSegments}] = local_coordinates(nodes, spotTable, idToRow);
        segmentNodes{nSegments} = nodes;
    end

    % Handles cycles or malformed components without silently dropping edges.
    for edgeIndex = find(~used).'
        nodes = [source(edgeIndex), target(edgeIndex)];
        used(edgeIndex) = true;
        nSegments = nSegments + 1;
        [xPaths{nSegments}, yPaths{nSegments}, tPaths{nSegments}] = local_coordinates(nodes, spotTable, idToRow);
        segmentNodes{nSegments} = nodes;
    end

    xPaths = xPaths(1:nSegments);
    yPaths = yPaths(1:nSegments);
    tPaths = tPaths(1:nSegments);
    segmentNodes = segmentNodes(1:nSegments);
    allSegmentCount = nSegments;
    keep = cellfun(@numel, tPaths) >= minimumFrames;
    xPaths = xPaths(keep);
    yPaths = yPaths(keep);
    tPaths = tPaths(keep);
    segmentNodes = segmentNodes(keep);
    metadata = struct('segmentNodes', {segmentNodes}, 'nSegments', nnz(keep), 'nSegmentsBeforeLengthFilter', allSegmentCount, 'minimumFrames', minimumFrames, 'nEdges', numel(source), 'nEdgesRetained', sum(cellfun(@(t) numel(t) - 1, tPaths)), 'nDivisionNodes', nnz(outDegree > 1));
end

function [x, y, t] = local_coordinates(nodes, spots, idToRow)
    present = arrayfun(@(id) isKey(idToRow, id), nodes);
    rows = cell2mat(values(idToRow, num2cell(nodes(present))));
    x = double(spots.POSITION_X(rows));
    y = double(spots.POSITION_Y(rows));
    t = double(spots.POSITION_T(rows));
    [t, order] = sort(t);
    x = x(order);
    y = y(order);
end
