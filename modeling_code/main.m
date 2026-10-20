close all; clc;
project = model_pipeline_setup();
fprintf('Outputs: %s\n', project.outputRoot);

%% Prepare experimental data
prepare_experimental_data(project, 'Recompute', false);
experiment = load_experimental_data(project);

%% Baseline models
results = run_baseline_model(project, 'Recompute', false); 

%% Figure 4B: Steady APELA profile
make_main_figures(results, experiment, project,'Panels', {'4B'});

%% Figure 4C: shAPLNR diffusivity from MSD
make_main_figures(results, experiment, project, 'Panels', {'4C'});

%% Figure 4D: Model-experiment inward velocity
make_main_figures(results, experiment, project, 'Panels', {'4D'});

%% Figure 4E: Model-experiment density by interval
make_main_figures(results, experiment, project, 'Panels', {'4E'});

%% Figure 4F: Experimental and model radial-speed maps
make_main_figures(results, experiment, project, 'Panels', {'4F'});

%% S4A-C: Parameter sweeps
sweeps = compute_supplementary_sweeps(experiment, project); 
make_supplementary_figures(results, experiment, project, 'Panels', {'S4A', 'S4B', 'S4C'}, 'Sweeps', sweeps);

%% S4D: Ligand and mesoderm boundary sensitivity
make_supplementary_figures(results, experiment, project, 'Panels', {'S4D'});

%% S4E: Inward velocity by experimental interval
make_supplementary_figures(results, experiment, project, 'Panels', {'S4E'});

%% S4F: Initial-condition histograms and Gaussian fits
make_supplementary_figures(results, experiment, project, 'Panels', {'S4F'});

%% S4G: Total, upper, and basal density
make_supplementary_figures(results, experiment, project, 'Panels', {'S4G'});

%% Export Figure 4 assembly panels
export_main_assembly_panels(project.mainFigureRoot);
