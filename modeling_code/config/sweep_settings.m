function cfg = sweep_settings()
    cfg.sourceRadius = linspace(180, 260, 21);
    cfg.lambda = linspace(5, 95, 21);
    cfg.beta = 2500:250:7500;
    cfg.mesodermDiffusivity = 10:2:50;
    cfg.alpha = 0:0.05:1;
    cfg.gamma = 0:0.0025:0.05;
end
