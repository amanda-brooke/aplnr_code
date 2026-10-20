function velocity = model_inward_velocity_profile( model, timeWindowHoursPostBMP4, distanceFromEdge, varargin)
    %MODEL_INWARD_VELOCITY_PROFILE Density-weighted drift.

    p = inputParser;
    addRequired(p, 'model', @isstruct);
    addRequired(p, 'timeWindowHoursPostBMP4', @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
    addRequired(p, 'distanceFromEdge', @isnumeric);
    addParameter(p, 'ExperimentStartHour', 24, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    parse(p, model, timeWindowHoursPostBMP4, distanceFromEdge, varargin{:});

    xModel = double(model.params.R) - model.rs(:);
    [xModel, order] = sort(xModel);
    upperDensity = local_time_radius(model.m_upper, model);
    basalDensity = local_time_radius(model.m_basal, model);
    upperVelocity = local_time_radius(model.v_upper, model);
    basalVelocity = local_time_radius(model.v_basal, model);
    elapsedWindow = double(timeWindowHoursPostBMP4) - double(p.Results.ExperimentStartHour);
    use = model.ts(:) >= elapsedWindow(1) & model.ts(:) <= elapsedWindow(2);
    if ~any(use)
        error('model_inward_velocity_profile:EmptyVelocityTimeWindow', 'The requested velocity interval does not overlap model output.');
    end
    velocityModel = (sum(upperDensity(use, :) .* (-upperVelocity(use, :)) + basalDensity(use, :) .* (-basalVelocity(use, :)), 1) ./ max(sum(upperDensity(use, :) + basalDensity(use, :), 1), eps)).';
    velocity = interp1(xModel, velocityModel(order), double(distanceFromEdge), 'linear', NaN);
end

function values = local_time_radius(values, model)
    if ~isequal(size(values), [numel(model.ts), numel(model.rs)])
        error('model_inward_velocity_profile:BadVelocitySeriesSize', 'Model state and velocity arrays must be time-by-radius.');
    end
end
