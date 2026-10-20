function stem = normalize_dataset_file_stem(baseName, label)
    %NORMALIZE_DATASET_FILE_STEM Normalize dataset CSV basenames to a shared stem.
    %
    % Examples:
    %   ctrl_1_edges   -> ctrl_1
    %   c1_spots_1000  -> c1_1000
    %   c1_tracks_1000 -> c1_1000

    baseName = char(string(baseName));
    label = lower(char(string(label)));
    baseLower = lower(baseName);

    if strcmp(baseLower, label)
        stem = '';
        return;
    end

    suffix = ['_' label];
    middleToken = ['_' label '_'];
    prefix = [label '_'];

    if endsWith(baseLower, suffix)
        stem = baseName(1:(end - strlength(suffix)));
    elseif contains(baseLower, middleToken)
        stem = regexprep(baseName, ['_' label '_'], '_', 'ignorecase', 'once');
    elseif startsWith(baseLower, prefix)
        stem = baseName((strlength(prefix) + 1):end);
    else
        stem = baseName;
    end

    stem = regexprep(stem, '_+', '_');
    stem = regexprep(stem, '^_', '');
    stem = regexprep(stem, '_$', '');
end
