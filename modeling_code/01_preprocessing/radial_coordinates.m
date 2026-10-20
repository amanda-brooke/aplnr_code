function [radius, distanceFromEdge, centerX, centerY] = radial_coordinates(x, y, geometry)
    %RADIAL_COORDINATES Coordinates relative to the fitted colony center.

    x = double(x(:)); y = double(y(:));
    assert(~isfield(geometry, 'centerTimes'), 'modeling:UnsupportedGeometry', 'Only fixed-center geometry is supported.');
    validateattributes(geometry.centerX, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(geometry.centerY, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(geometry.edgeRadius, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    assert(numel(x) == numel(y), 'modeling:CoordinateSizeMismatch', 'X and Y coordinates must have equal lengths.');
    centerX = repmat(double(geometry.centerX), size(x));
    centerY = repmat(double(geometry.centerY), size(x));
    radius = hypot(x - centerX, y - centerY);
    distanceFromEdge = double(geometry.edgeRadius) - radius;
end
