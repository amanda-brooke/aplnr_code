function data = load_experimental_data(project)
    %LOAD_EXPERIMENTAL_DATA Load summaries using the configured raw-data location.

    paths = prepare_experimental_data(project, 'Recompute', false);
    data.control = local_load(paths.control, project, 'control');
    data.shAPLNR = local_load(paths.treatment, project, 'treatment');
    data.paths = paths;
end

function replicates = local_load(path, project, condition)
    loaded = load(path, 'replicateSummaries');
    replicates = loaded.replicateSummaries;
    fields = {'spotFile', 'edgeFile', 'trackFile'};
    folders = {'spots', 'edges', 'tracks'};
    for i = 1:numel(replicates)
        for j = 1:numel(fields)
            if ~isfield(replicates(i).files, fields{j}), continue; end
            oldPath = char(string(replicates(i).files.(fields{j})));
            if isempty(oldPath), continue; end
            oldPath = strrep(oldPath, '\', filesep);
            [~, name, extension] = fileparts(oldPath);
            candidate = fullfile(project.rawDataRoot, condition, folders{j}, [name extension]);
            % Cached absolute paths must not select another copy of the data.
            replicates(i).files.(fields{j}) = candidate;
        end
    end
end
