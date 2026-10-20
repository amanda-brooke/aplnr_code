function apparentSpeed = model_apparent_radial_speed(radialDrift, diffusivity, frameIntervalHours)

    if ~(isnumeric(diffusivity) && all(isfinite(diffusivity(:))) && all(diffusivity(:) >= 0) && (isscalar(diffusivity) || isequal(size(diffusivity), size(radialDrift))))
        error('model_apparent_radial_speed:BadDiffusivity', 'Diffusivity must be finite, nonnegative, and scalar or the same size as radialDrift.');
    end
    if ~(isnumeric(frameIntervalHours) && isscalar(frameIntervalHours) && isfinite(frameIntervalHours) && frameIntervalHours > 0)
        error('model_apparent_radial_speed:BadFrameInterval', 'Frame interval must be a finite positive scalar in hours.');
    end

    mu = radialDrift .* frameIntervalHours;
    sigma = sqrt(2 .* diffusivity .* frameIntervalHours);
    if all(sigma(:) == 0)
        apparentSpeed = abs(radialDrift);
        return;
    end

    z = abs(mu) ./ sigma;
    expectedAbsStep = sigma .* sqrt(2 / pi) .* exp(-0.5 .* z.^2) + abs(mu) .* erf(z ./ sqrt(2));
    apparentSpeed = expectedAbsStep ./ frameIntervalHours;
    zeroSigma = sigma == 0;
    if any(zeroSigma(:))
        apparentSpeed(zeroSigma) = abs(radialDrift(zeroSigma));
    end
end
