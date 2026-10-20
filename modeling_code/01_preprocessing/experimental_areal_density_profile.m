function [r, t, areaDensity] = experimental_areal_density_profile(density)

    arguments
        density (1, 1) struct
    end

    r = density.rCenters(:); t = density.tHours(:);
    assert(strcmp(density.densityMode, 'area'), 'experimental_areal_density_profile:UnsupportedDensityMode', 'Expected densityMode=''area'' from prepare_experimental_data.');
    areaDensity = density.densityOver;
    if ~isequal(size(areaDensity), [numel(t), numel(r)])
        error('experimental_areal_density_profile:BadDensitySize', 'Expected time-by-radius size [%d %d], got [%s].', numel(t), numel(r), num2str(size(areaDensity)));
    end
end
