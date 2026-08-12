function tbl = otr_read_numeric_table(filename, delimiter)
%OTR_READ_NUMERIC_TABLE Read a small numeric text/CSV table robustly.
%   TBL = OTR_READ_NUMERIC_TABLE(FILENAME) skips blank/metadata lines until
%   it finds a row containing at least two numeric fields.  A preceding
%   nonnumeric row of matching width is retained in TBL.headers.  UTF-8 BOM
%   bytes are removed.  No import or signal-processing toolbox is required.
%
%   The returned struct contains:
%       data        numeric matrix (short rows are padded with NaN)
%       headers     cell array of column labels (possibly empty strings)
%       metadata    raw text lines preceding the numeric data
%       filename    absolute or supplied input path
%       delimiter   delimiter used to split records

%   This reader targets MATLAB R2019b and newer and intentionally avoids
%   readtable/detectImportOptions, whose inferred rules vary by release.

%   See also OTR_READ_SCOPE_CSV.

    if nargin < 1 || isempty(filename)
        error('otr:MissingFilename', 'A filename is required.');
    end
    if nargin < 2 || isempty(delimiter)
        delimiter = local_detect_delimiter(filename);
    end
    if exist(filename, 'file') ~= 2
        error('otr:FileNotFound', 'Cannot find numeric table: %s', filename);
    end

    fid = fopen(filename, 'rt');
    if fid < 0
        error('otr:FileOpenFailed', 'Cannot open numeric table: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = cell(0, 1);
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if ~isempty(line)
            % Remove a UTF-8 BOM.  Depending on the locale MATLAB exposes
            % it either as Unicode U+FEFF or as its three raw UTF-8 bytes.
            if double(line(1)) == 65279
                line(1) = [];
            elseif numel(line) >= 3 && isequal(double(line(1:3)), [239 187 191])
                line(1:3) = [];
            end
        end
        lines{end + 1, 1} = line; %#ok<AGROW>
    end

    numericStart = [];
    numericRows = cell(0, 1);
    rowWidths = zeros(0, 1);
    for ii = 1:numel(lines)
        fields = local_split(lines{ii}, delimiter);
        [values, validCount] = local_numeric_fields(fields);
        if validCount >= 2
            if isempty(numericStart)
                numericStart = ii;
            end
            numericRows{end + 1, 1} = values; %#ok<AGROW>
            rowWidths(end + 1, 1) = numel(values); %#ok<AGROW>
        elseif ~isempty(numericStart)
            % Ignore blank/footer records after the numeric body.  A
            % nonnumeric record in the middle must not shift later rows.
            continue;
        end
    end

    if isempty(numericStart)
        error('otr:NoNumericData', 'No row with two numeric fields in %s.', filename);
    end

    ncol = max(rowWidths);
    data = NaN(numel(numericRows), ncol);
    for ii = 1:numel(numericRows)
        values = numericRows{ii};
        data(ii, 1:numel(values)) = values;
    end

    % Remove columns that are completely empty.  This also handles a
    % trailing delimiter without changing intentional all-NaN gaps inside.
    keep = any(isfinite(data), 1);
    data = data(:, keep);

    headers = repmat({''}, 1, ncol);
    headerLine = numericStart - 1;
    while headerLine >= 1 && isempty(strtrim(lines{headerLine}))
        headerLine = headerLine - 1;
    end
    if headerLine >= 1
        candidate = local_split(lines{headerLine}, delimiter);
        [~, nNumeric] = local_numeric_fields(candidate);
        if nNumeric < 2
            m = min(numel(candidate), ncol);
            headers(1:m) = candidate(1:m);
        end
    end
    headers = headers(keep);
    for ii = 1:numel(headers)
        headers{ii} = strtrim(strrep(headers{ii}, '"', ''));
    end

    tbl = struct();
    tbl.data = data;
    tbl.headers = headers;
    tbl.metadata = lines(1:max(0, numericStart - 1));
    tbl.filename = filename;
    tbl.delimiter = delimiter;
end

function delimiter = local_detect_delimiter(filename)
    fid = fopen(filename, 'rt');
    if fid < 0
        error('otr:FileOpenFailed', 'Cannot open numeric table: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    sample = '';
    for ii = 1:25
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        sample = [sample line]; %#ok<AGROW>
    end
    if numel(strfind(sample, sprintf('\t'))) > numel(strfind(sample, ',')) %#ok<STREMP>
        delimiter = sprintf('\t');
    elseif ~isempty(strfind(sample, ',')) %#ok<STREMP>
        delimiter = ',';
    elseif ~isempty(strfind(sample, ';')) %#ok<STREMP>
        delimiter = ';';
    else
        delimiter = 'whitespace';
    end
end

function fields = local_split(line, delimiter)
    line = strtrim(line);
    if isempty(line)
        fields = {''};
    elseif strcmp(delimiter, 'whitespace')
        fields = regexp(line, '\s+', 'split');
    else
        fields = regexp(line, regexptranslate('escape', delimiter), 'split');
    end
end

function [values, validCount] = local_numeric_fields(fields)
    values = NaN(1, numel(fields));
    validCount = 0;
    for jj = 1:numel(fields)
        token = strtrim(strrep(fields{jj}, '"', ''));
        if isempty(token)
            continue;
        end
        value = sscanf(token, '%f', 1);
        if ~isempty(value)
            values(jj) = value;
            validCount = validCount + 1;
        end
    end
end
