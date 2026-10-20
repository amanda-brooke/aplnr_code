function project = project_paths(root)
    if nargin == 0
        root = fileparts(fileparts(mfilename('fullpath')));
    end
    project.root = root;
    project.rawDataRoot = getenv('MODELING_DATA_ROOT');
    if isempty(project.rawDataRoot)
        project.rawDataRoot = fullfile(root, 'data', 'offset_24');
    end
    project.outputRoot = getenv('MODELING_OUTPUT_ROOT');
    if isempty(project.outputRoot)
        project.outputRoot = fullfile(root, 'outputs');
    end
    project.preparedDataRoot = fullfile(project.outputRoot, 'prepared_data', 'offset_24');
    project.resultsRoot = fullfile(project.outputRoot, 'results');
    project.mainFigureRoot = fullfile(project.outputRoot, 'main_figures');
    project.supplementaryRoot = fullfile(project.outputRoot, 'supplementary');
end
