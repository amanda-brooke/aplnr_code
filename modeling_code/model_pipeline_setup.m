function project = model_pipeline_setup()
    root = fileparts(mfilename('fullpath'));
    folders = {'config', '01_preprocessing', '02_model', '03_figures', '04_supplement'};
    addpath(root);
    for i = 1:numel(folders)
        addpath(fullfile(root, folders{i}));
    end
    project = project_paths(root);
    for folder = {project.outputRoot, project.resultsRoot, project.mainFigureRoot, project.supplementaryRoot}
        if ~isfolder(folder{1}), mkdir(folder{1}); end
    end
end
