function [r, t, densityProfile] = model_areal_density_profile(model)
    %MODEL_AREAL_DENSITY_PROFILE Return unit-mass model density profiles.
    %
    % Each stored upper-plus-basal profile is normalized by its annular
    % mass, integral 2*pi*r*density dr, matching the offset-24 figure workflow.

    r = model.rs(:);
    t = model.ts(:);
    expectedSize = [numel(t), numel(r)];
    assert(isequal(size(model.m_upper), expectedSize) && isequal(size(model.m_basal), expectedSize), 'model_areal_density_profile:BadSeriesSize', 'Upper and basal density arrays must be time-by-radius.');
    totalDensity = max(model.m_upper, 0) + max(model.m_basal, 0);

    densityProfile = nan(numel(t), numel(r));
    for i = 1:numel(t)
        profile = max(totalDensity(i, :).', 0);
        totalMass = trapz(r, 2 * pi * r .* profile);
        if isfinite(totalMass) && totalMass > 0
            densityProfile(i, :) = (profile ./ totalMass).';
        end
    end
end

