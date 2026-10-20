function files = resolve_dataset_files(datasetGroup, datasetName, varargin)
    % Resolve edge/spot/track CSVs directly under the configured data root.
    %
    % Expected layout:
    %   DataRoot/<datasetName>/edges/*.csv
    %   DataRoot/<datasetName>/spots/*.csv
    %   DataRoot/<datasetName>/tracks/*.csv (optional)
    % datasetGroup identifies the dataset in returned metadata.

    parser = inputParser;
    parser.FunctionName = mfilename;
    addRequired(parser, 'datasetGroup', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addRequired(parser, 'datasetName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'DataRoot', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data'), @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'FileStem', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parse(parser, datasetGroup, datasetName, varargin{:});

    datasetGroup = char(string(parser.Results.datasetGroup));
    datasetName = char(string(parser.Results.datasetName));
    dataRoot = char(string(parser.Results.DataRoot));
    fileStem = char(string(parser.Results.FileStem));

    if ~isfolder(dataRoot)
        error('resolve_dataset_files:MissingDataRoot', 'Data root does not exist: %s', dataRoot);
    end

    datasetDir = fullfile(dataRoot, datasetName);
    if ~isfolder(datasetDir)
        error('resolve_dataset_files:NotFound', 'Dataset directory does not exist: %s', datasetDir);
    end
    trackDir = fullfile(datasetDir, 'tracks');

    files = struct( 'datasetGroup', datasetGroup, 'datasetName', datasetName, 'fileStem', fileStem, 'dataRoot', dataRoot, 'datasetDir', datasetDir, 'edgeDir', fullfile(datasetDir, 'edges'), 'spotDir', fullfile(datasetDir, 'spots'), 'trackDir', trackDir, 'edgeFile', local_single_csv(fullfile(datasetDir, 'edges'), 'edges', fileStem), 'spotFile', local_single_csv(fullfile(datasetDir, 'spots'), 'spots', fileStem), 'trackFile', '' );

    if isfolder(trackDir)
        files.trackFile = local_optional_single_csv(trackDir, 'tracks', fileStem);
    end
end

function filePath = local_single_csv(folderPath, label, fileStem)
    listing = dir(fullfile(folderPath, '*.csv'));
    listing = listing(~[listing.isdir]);
    listing = local_filter_listing(listing, fileStem);

    if isempty(listing)
        allListing = dir(fullfile(folderPath, '*.csv'));
        allListing = allListing(~[allListing.isdir]);
        names = {allListing.name};
        if isempty(names), names = {'<none>'}; end
        if isempty(fileStem)
            error(local_folder_error_id(label), 'No CSV files found in %s.', folderPath);
        end
        error(local_folder_error_id(label), 'No CSV in %s matched FileStem "%s". Available files:\n%s', folderPath, fileStem, strjoin(names, newline));
    end

    if numel(listing) ~= 1
        names = {listing.name};
        if isempty(fileStem)
            error(local_folder_error_id(label), ['Expected exactly one CSV in %s. Found %d:\n%s\n' 'Pass FileStem to select one, e.g. Analysis(''offset_24'', ''control'', ''FileStem'', ''ctrl_3'').'], folderPath, numel(listing), strjoin(names, newline));
        end
        error(local_folder_error_id(label), 'FileStem "%s" matched multiple CSVs in %s:\n%s', fileStem, folderPath, strjoin(names, newline));
    end
    filePath = fullfile(listing(1).folder, listing(1).name);
end

function filePath = local_optional_single_csv(folderPath, label, fileStem)
    listing = dir(fullfile(folderPath, '*.csv'));
    listing = listing(~[listing.isdir]);
    listing = local_filter_listing(listing, fileStem);
    if isempty(listing)
        filePath = '';
        return;
    end
    if numel(listing) > 1
        names = {listing.name};
        error(local_folder_error_id(label), 'Expected at most one CSV in %s. Found %d:\n%s', folderPath, numel(listing), strjoin(names, newline));
    end
    filePath = fullfile(listing(1).folder, listing(1).name);
end

function listing = local_filter_listing(listing, fileStem)
    if isempty(fileStem)
        return;
    end

    stem = lower(strtrim(fileStem));
    keep = false(size(listing));
    for i = 1:numel(listing)
        [~, base, ~] = fileparts(listing(i).name);
        [~, label] = fileparts(char(string(listing(i).folder)));
        normalizedStem = normalize_dataset_file_stem(base, label);
        keep(i) = strcmpi(base, stem) || strcmpi(normalizedStem, stem) || startsWith(lower(normalizedStem), [stem '_']);
    end
    listing = listing(keep);
end

function errId = local_folder_error_id(label)
    errId = ['resolve_dataset_files:Bad' upper(label(1)) label(2:end) 'Folder'];
end
