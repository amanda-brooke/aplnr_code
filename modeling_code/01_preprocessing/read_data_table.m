function T = read_data_table(fileName, varargin)
    % Read a CSV at the path returned by resolve_dataset_files.
    if ~isfile(fileName)
        error('read_data_table:FileNotFound', 'CSV file does not exist: %s', fileName);
    end
    T = readtable(fileName, varargin{:});
end
